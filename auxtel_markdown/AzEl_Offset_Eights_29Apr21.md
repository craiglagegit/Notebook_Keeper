---
jupyter:
  jupytext:
    formats: ipynb,markdown//md
    text_representation:
      extension: .md
      format_name: markdown
      format_version: '1.3'
      jupytext_version: 1.14.0
  kernelspec:
    display_name: LSST
    language: python
    name: lsst
---

## AuxTel AzEl offsets - 28-Apr-21

In this notebook, investigate az-el offsets from 11-Mar-21\
This is going to analyze a pair of exposures

```python
import sys, time, os, asyncio, glob

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pickle as pkl
import pandas as pd
import astropy.io.fits as pf
from astropy.time import Time, TimeDelta
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, AltAz, ICRS, EarthLocation, Angle, FK5
import astropy.units as u

from lsst.daf.butler import Butler as gen3Butler
from lsst.daf.persistence import Butler as gen2Butler
from lsst_efd_client import EfdClient
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python jupyter={"source_hidden": true}
# Get EFD client
client = EfdClient('summit_efd')

def merge_packed_time_series(packed_dataframe, base_field, stride=1, 
                             ref_timestamp_col="cRIO_timestamp", internal_time_scale="tai"):
    """Select fields that are time samples and unpack them into a dataframe.
            Parameters
            ----------
            packedDF : `pandas.DataFrame`
                packed data frame containing the desired data
            base_field :  `str`
                Base field name that will be expanded to query all
                vector entries.
            stride : `int`, optional
                Only use every stride value when unpacking.  Must be a factor
                of the number of packed values.
                (1 by default)
            ref_timestamp_col : `str`, optional
                Name of the field name to use to assign timestamps to unpacked
                vector fields (default is 'cRIO_timestamp').
            internal_time_scale : `str`, optional
                Time scale to use when converting times to internal formats
                ('tai' by default). Equivalent to EfdClient.internal_scale
        Returns
            -------
            result : `pandas.DataFrame`
                A `pandas.DataFrame` containing the results of the query.
            """
    
    packed_fields = [k for k in packed_dataframe.keys() if k.startswith(base_field)]
    packed_fields = sorted(packed_fields, key=lambda k: int(k[len(base_field):]))  # sort by pack ID
    npack = len(packed_fields)
    if npack%stride != 0:
        raise RuntimeError(f"Stride must be a factor of the number of packed fields: {stride} v. {npack}")
    packed_len = len(packed_dataframe)
    n_used = npack//stride   # number of raw fields being used
    output = np.empty(n_used*packed_len)
    times = np.empty_like(output, dtype=packed_dataframe[ref_timestamp_col][0])
    
    if packed_len == 1:
        dt = 0
    else:
        dt = (packed_dataframe[ref_timestamp_col][1] - packed_dataframe[ref_timestamp_col][0])/npack
    for i in range(0, npack, stride):
        i0 = i//stride
        output[i0::n_used] = packed_dataframe[f"{base_field}{i}"]
        times[i0::n_used] = packed_dataframe[ref_timestamp_col] + i*dt
     
    timestamps = Time(times, format='unix', scale=internal_time_scale).datetime64
    return pd.DataFrame({base_field:output, "times":times}, index=timestamps)

# This makes the plots more readable
def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar
```

```python
# Gen3 butler
dayObs = 20210311
REPO_DIR = '/repo/main'
butler_g3 = gen3Butler(REPO_DIR, collections="LATISS/raw/all")
```

```python jupyter={"outputs_hidden": true} tags=[]
# First get all of the tracking data

firstExpId = 2021031100418
lastExpId = 2021031100424
charVisits = {}
for expId in range(firstExpId, lastExpId+1):
    charVisit = {}
    expId1 = expId-1
    expId2 = expId
    mData1 = butler_g3.get('raw.metadata', detector=0, exposure=expId1)
    exp1Start = Time(mData1['DATE-BEG'],scale='tai')
    exp1End = Time(mData1['DATE-END'],scale='tai')
    mData2 = butler_g3.get('raw.metadata', detector=0, exposure=expId2)
    exp2Start = Time(mData2['DATE-BEG'],scale='tai')
    exp2End = Time(mData2['DATE-END'],scale='tai')
    charVisit['mData'] = mData2
    # These are for finding the timestamps of the offset event
    offsets = await client.select_time_series("lsst.sal.ATPtg.command_offsetAzEl", ['*'],
                                          exp1Start, exp2End)
    for i, offset in enumerate(offsets.values):
        offsetTime = Time(offsets.index[0],scale='tai')
        commandedAzShift = offsets.values[0][0]
        commandedElShift = offsets.values[0][1]
        charVisit['commandedAzShift'] = commandedAzShift
        charVisit['commandedElShift'] = commandedElShift
    # Now get the mount data
    mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                          exp1Start, exp2End)
    az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1)
    el = merge_packed_time_series(mount_position, 'elevationCalculatedAngle', stride=1)

    # Calculate the tracking shift
    az_vals = np.array(az.values.tolist())[:,0]
    el_vals = np.array(el.values.tolist())[:,0]
    times = np.array(az.values.tolist())[:,1]

    shiftIndex = np.where(times>offsetTime.tai.unix)[0][0]
    valsPerSecond = 100
    n1 = np.where(times>exp1Start.tai.unix)[0][0] # Start of fist exposure
    n2 = np.where(times>exp1End.tai.unix)[0][0] # End of first exposure
    n3 = np.where(times>exp2Start.tai.unix)[0][0] # Start of second exposure
    try:
        n4 = np.where(times>exp2End.tai.unix)[0][0] # End of second exposure
    except:
        n4 = len(times)

    # Fit the tracking before the offset with a quadratic
    az_vals_fit = az_vals[n1:n2]
    el_vals_fit = el_vals[n1:n2]
    times_fit = times[n1:n2]

    az_fit = np.polyfit(times_fit, az_vals_fit, 2)
    el_fit = np.polyfit(times_fit, el_vals_fit, 2)

    az_model = az_fit[0] * times * times + az_fit[1] * times + az_fit[2]
    el_model = el_fit[0] * times * times + el_fit[1] * times + el_fit[2]

    # Apply this model and calculate the departure from the model after the shift
    az_error = (az_vals - az_model) * 3600
    el_error = (el_vals - el_model) * 3600

    az_shift = np.mean(az_error[n3:n4])
    el_shift = np.mean(el_error[n3:n4])
    charVisit['az_shift'] = az_shift
    charVisit['el_shift'] = el_shift
    charVisits[f"{expId}"] = charVisit
    
    fig = plt.figure(figsize=(8,8))
    plt.subplots_adjust(hspace=0.4)
    plt.subplot(3,1,1)
    plt.title(f"Commands - 11-Mar-21", fontsize = 18)
    plt.plot(times, times*0.0-1.0, color='green')
    plt.plot([exp1Start.tai.unix, exp1End.tai.unix], [0.1,0.1], color='red', lw=2.0, label=f"exp:{expId1}")
    plt.plot([exp2Start.tai.unix, exp2End.tai.unix], [0.1,0.1], color='red', lw=2.0, label=f"exp:{expId2}")
    plt.plot([offsetTime.tai.unix,offsetTime.tai.unix], [0.1, 0.4], color='blue', lw=2.0, label="Commanded Offset")
    plt.ylim(0.0, 1.0)
    plt.legend()
    plt.subplot(3,1,2)
    plt.title(f"Azimuth change", fontsize = 18)
    plt.plot(times, az_error, color='green')
    plt.text(np.mean(times), np.mean(az_error), f"az shift = {az_shift:.2f} arcsec")
    plt.text(np.mean(times), np.mean(az_error)*0.5, f"commanded az shift = {commandedAzShift:.1f} arcsec")
    plt.ylabel("Arcsec")
    plt.subplot(3,1,3)
    plt.title(f"Elevation change", fontsize = 18)
    plt.plot(times, el_error, color='green')
    plt.text(np.mean(times), np.mean(el_error), f"el shift = {el_shift:.2f} arcsec")
    plt.text(np.mean(times), np.mean(el_error)*0.5, f"commanded el shift = {commandedElShift:.1f} arcsec")
    plt.ylabel("Arcsec")
    plt.xlabel("Time (Unix - seconds)")
    plt.savefig(f"/project/cslage/AuxTel/offsets/Offsets_{expId1}_{expId2}_29Apr21.pdf")
```

```python

```

```python jupyter={"outputs_hidden": true}
# Now get the image data
charConfig = CharacterizeImageConfig()
charConfig.doMeasurePsf = False#True
charConfig.doApCorr = False
charConfig.doDeblend = False
charConfig.repair.doCosmicRay = True
charConfig.repair.doInterpolate = True   
charConfig.detection.minPixels = 500
charTask = CharacterizeImageTask(config=charConfig)

for expId in range(firstExpId, lastExpId+1):
    charVisit = charVisits[f"{expId}"]
    exp = gen2_butler.get('quickLookExp', detector=0, expId=expId)
    charResult = charTask.run(exp)
    sourceCatalog = charResult.sourceCat
    maxFlux = np.nanmax(sourceCatalog['base_CircularApertureFlux_3_0_instFlux'])
    selectBrightestSource = sourceCatalog['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.99
    brightestSource = sourceCatalog.subset(selectBrightestSource)
    brightestCentroid = (brightestSource['base_SdssCentroid_x'][0], \
                         brightestSource['base_SdssCentroid_y'][0])
    brightCatalog = sourceCatalog.subset(sourceCatalog['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.001)
    print(f"expId:{expId}. Found {len(sourceCatalog)} sources, {len(brightCatalog)} bright sources")
    print(f"Brightest centroid at {brightestCentroid}")
    charVisit['exp'] = exp
    charVisit['brightestCentroid'] = brightestCentroid
    charVisit['brightCatalog'] = brightCatalog
    charVisits[f"{expId}"] = charVisit

```

```python
# Plot the images
Ncenter = (3000,3000)
Nlength = 500.0
yShift = 150.0
shiftCounter = 0
plt.figure(figsize=(8,8))

for expId in range(firstExpId, lastExpId+1):
    seqNo = expId - 2021031100000
    charVisit = charVisits[f"{expId}"]
    plt.subplot(1,1,1)
    plt.title(f"Figure Eight - {firstExpId} - {lastExpId}",fontsize=18)
    plt.scatter([charVisit['brightestCentroid'][0]],[charVisit['brightestCentroid'][1]] \
                ,color='yellow', marker='+', s=400)
    if seqNo == 418:
        arr = charVisit['exp'].image.array
        arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
        img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')

        cat = charVisit['brightCatalog']
        rotpa = Angle(charVisit['mData']['ROTPA']*u.deg)
        plt.arrow(Ncenter[0],Ncenter[1], -Nlength*np.sin(rotpa), Nlength*np.cos(rotpa),\
            color='green', width = 20)
        plt.text(Ncenter[0]-Nlength*np.sin(rotpa)+100,Ncenter[1]+Nlength*np.cos(rotpa)+100, 'N', \
            color='green', fontsize=12, weight='bold')
        plt.arrow(Ncenter[0],Ncenter[1], Nlength*np.cos(rotpa), Nlength*np.sin(rotpa),\
            color='green', width = 20)
        plt.text(Ncenter[0]+Nlength*np.cos(rotpa)+100,Ncenter[1]+Nlength*np.sin(rotpa)+100, 'E', \
            color='green', fontsize=12, weight='bold')

    if seqNo in [418, 421, 424]:
        plt.text((charVisit['brightestCentroid'][0]+100), \
                 (charVisit['brightestCentroid'][1] + shiftCounter*yShift), \
                 f"{seqNo} az={charVisit['commandedAzShift']} el={charVisit['commandedElShift']}", \
                 color='yellow', fontsize=12, weight='bold')
        
        shiftCounter += 1
    else:
        plt.text((charVisit['brightestCentroid'][0]+100), \
                 (charVisit['brightestCentroid'][1]), \
                 f"{seqNo} az={charVisit['commandedAzShift']}\n       el={charVisit['commandedElShift']}", \
                 color='yellow', fontsize=12, weight='bold')
    #plt.xlim(0,4000)
    #plt.ylim(4000,0)
plt.tight_layout(h_pad=1)
plt.savefig(f"/project/cslage/AuxTel/offsets/Offsets_Figure_Eight_{firstExpId}_{lastExpId}_29Apr21.pdf")
```

```python
outfile = open('/project/cslage/AuxTel/offsets/offsets_figure_eight_418_424_29apr21.pkl','wb')

pkl.dump(charVisits,outfile)
outfile.close()
```

```python

```
