#!/usr/bin/env python

# Code for classifying AuxTel mount errors
# 
# Craig Lage - 28Dec23

import nest_asyncio
nest_asyncio.apply()
import sys, time, os, asyncio
from datetime import datetime
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta

from lsst_efd_client import EfdClient
from lsst.daf.butler import Butler

import lsst.summit.utils.butlerUtils as butlerUtils
from lsst.summit.utils.efdUtils import calcNextDay
from lsst.summit.utils.utils import dayObsIntToString
from astro_metadata_translator import ObservationInfo
from lsst_efd_client import merge_packed_time_series as mpts

from scipy.interpolate import UnivariateSpline
from scipy.fft import fft, fftfreq


###############SUBROUTINES################

NON_TRACKING_IMAGE_TYPES = ['BIAS',
                            'FLAT',
                            'DARK']

AUXTEL_ANGLE_TO_EDGE_OF_FIELD_ARCSEC = 280.0
MOUNT_IMAGE_WARNING_LEVEL = .25  # this determines the colouring of the cells in the table, yellow for this
MOUNT_IMAGE_BAD_LEVEL = .4


def _getEfdData(client, dataSeries, startTime, endTime):
    """A synchronous warpper for geting the data from the EFD.

    This exists so that the top level functions don't all have to be async def.
    """
    loop = asyncio.get_event_loop()
    return loop.run_until_complete(client.select_time_series(dataSeries, ['*'], startTime.utc, endTime.utc))

def findSpline(xf, yf, numKnots=18):
    # Finds a spline with a specified number of knots
    # by binary search through smoothing values
    s1 = 1; s2 = 1.0E7
    spline = UnivariateSpline(xf, yf, s=s1)
    knots = spline.get_knots()
    if len(knots) < numKnots:
        return s1, knots
    count = 0
    while count < 50:
        count += 1
        s = np.sqrt(s1 * s2)
        spline = UnivariateSpline(xf, yf, s=s)
        knots = spline.get_knots()
        if len(knots) > numKnots:
            s1 = s
        else:
            s2 = s
        if abs(len(knots) - numKnots) < 2:
            return s, knots
    return s,knots

def calculateFFTKnots(dataId, butler, client, limit=0.25, fig=None):
    expRecord = butlerUtils.getExpRecordFromDataId(butler, dataId)
    dayString = dayObsIntToString(expRecord.day_obs)
    seqNumString = str(expRecord.seq_num)
    dataIdString = f"{dayString} - seqNum {seqNumString}"

    imgType = expRecord.observation_type.upper()
    if imgType in NON_TRACKING_IMAGE_TYPES:
        return False

    exptime = expRecord.exposure_time
    if exptime < 4.99:
        return False

    tStart = expRecord.timespan.begin.tai.to_value("isot")
    tEnd = expRecord.timespan.end.tai.to_value("isot")
    elevation = 90.0 - expRecord.zenith_angle

    # TODO: DM-33859 remove this once it can be got from the expRecord
    md = butler.get('raw.metadata', dataId, detector=0)
    obsInfo = ObservationInfo(md)
    azimuth = obsInfo.altaz_begin.az.value

    t_start = Time(tStart, scale='tai')
    t_end = Time(tEnd, scale='tai')

    mount_position = _getEfdData(client, "lsst.sal.ATMCS.mount_AzEl_Encoders", t_start, t_end)

    az = mpts(mount_position, 'azimuthCalculatedAngle', stride=1)
    el = mpts(mount_position, 'elevationCalculatedAngle', stride=1)

    # Calculate the tracking errors
    az_vals = np.array(az.values[:, 0])
    el_vals = np.array(el.values[:, 0])
    times = np.array(az.values[:, 1])
    # The fits are much better if the time variable
    # is centered in the interval
    fit_times = times - times[int(len(az.values[:, 1]) / 2)]

    # Fit with a polynomial
    az_fit = np.polyfit(fit_times, az_vals, 4)
    el_fit = np.polyfit(fit_times, el_vals, 4)
    az_model = np.polyval(az_fit, fit_times)
    el_model = np.polyval(el_fit, fit_times)

    # Errors in arcseconds
    az_error = (az_vals - az_model) * 3600
    el_error = (el_vals - el_model) * 3600

    # Calculate RMS
    az_rms = np.sqrt(np.mean(az_error * az_error))
    el_rms = np.sqrt(np.mean(el_error * el_error))

    # Calculate Image impact RMS
    image_az_rms = az_rms * np.cos(el_vals[0] * np.pi / 180.0)
    image_el_rms = el_rms
    tot_rms = np.sqrt(image_az_rms**2 + image_el_rms**2)

    if tot_rms < limit:
        return [tot_rms, None, True]
    else:
        # Check the timebase errors
        # The classifier has difficulty distinguishing these, so I added a dedicated test
        cRIO_ts = mount_position["cRIO_timestamp"]
        timestamps = cRIO_ts.values
        deltaTs = []
        for n in range(1, len(timestamps)):
            deltaTs.append(timestamps[n] - timestamps[n-1])
        if np.max(deltaTs) > 1.05:
            timebaseFlag = True
        else:
            timebaseFlag = False

        if fig:
            axs = fig.subplots(1,2)
            plt.subplots_adjust(wspace=0.3)
            names = ["Azimuth", "Elevation"]
        # Calculate the FFT knots
        knotList = []
        numKnots = 18
        for i, error in enumerate([az_error, el_error]):
            # Number of samples in time series
            N = len(error)
            SAMPLE_RATE = 100 # Samples/sec
            yf = fft(error)
            xf = fftfreq(N, 1 / SAMPLE_RATE)
            yf = np.abs(fft(error))
            # Truncate the FFT above 5 Hz
            count = 0
            for n in range(len(xf)):
                if xf[n] < 5.0:
                    count += 1
                else:
                    break
            xf = xf[0:count]
            yf = yf[0:count]

            s, knots = findSpline(xf, yf, numKnots=numKnots)
            spline = UnivariateSpline(xf, yf, s=s)
            knots = spline.get_knots()
            values = spline(knots)
            # Truncate or add to get the desired number of knots
            if len(knots) < numKnots:
                numToAdd = int(numKnots - len(knots))
                lenKnots = knots[-1] - knots[0]
                for n in range(numToAdd):
                    newKnot = knots[-1] + 0.01 * (n+1)
                    newValue = values[-1]
                    knots = np.append(knots, newKnot)
                    values = np.append(values, newValue)
            if len(knots) > numKnots:
                knots = np.delete(knots, [-2])
                values = np.delete(values, [-2])
            for n in range(numKnots):
                knotList.append(knots[n])
                knotList.append(values[n])
            if fig:
                axs[i].set_title(names[i])
                axs[i].plot(xf, yf, color='blue')
                axs[i].plot(xf, spline(xf), ls='--', color='red')
                axs[i].scatter(knots, values, marker='x', color='red')
                axs[i].set_xlabel("Frequency(Hz)")
                axs[i].set_ylabel("Amplitude")
        return [tot_rms, knotList, timebaseFlag]

