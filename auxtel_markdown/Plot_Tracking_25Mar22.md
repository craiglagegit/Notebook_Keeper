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

# AuxTel Plot tracking - 25-Mar-22

In this notebook, investigate again mount tracking events in Feb-Mar, 2022\
This is after the EFD was converted to UTC.\
Thanks to Simon Krughoff for contributions.

```python
import sys, time, os, asyncio
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from astropy.time import Time, TimeDelta
from lsst.daf.butler import Butler
from lsst_efd_client import EfdClient
```

```python
# Get EFD client and the butler
client = EfdClient('ldf_stable_efd')
butler = Butler('/repo/main', collections="LATISS/raw/all")
```

```python
expIdList = [2022021600091,2022021600097,2022021600101,2022021600107,2022021600109,\
            2022021600721,2022021700299,2022021700311,2022031600536,2022031600795,\
            2022031700454,2022031700461]

doOffset = True

pdf = PdfPages("/project/cslage/AuxTel/offsets/In_Position_Errors_29Mar22.pdf")
for expId in expIdList:
    try:
        mData = butler.get('raw.metadata', detector=0, exposure=expId)
        # Need to convert DATE_BEG and DATE_END to UTC to sync up with the EFD
        date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
        date_end = Time(mData['DATE-END'], format='isot', scale='tai')
        # Use these for finding the "allAxesInPosition" timestamp
        # The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
        before = 10.0
        after = 30.0
        if expId < 2021101300000:
            # EFD was switched to UTC on 20211013
            tai_offset = 37.0
        else:
            tai_offset = 0.0

        start = date_beg.utc - TimeDelta(before, format='sec') + TimeDelta(tai_offset, format='sec')
        end = date_end.utc + TimeDelta(after, format='sec') + TimeDelta(tai_offset, format='sec')
        inPosition = await client.select_time_series("lsst.sal.ATMCS.logevent_allAxesInPosition", "inPosition", start, end)
        inPosition = inPosition[inPosition['inPosition']==True] 

        # Use these for finding the shutter status timestamp
        # The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
        # They agree within a few milliseconds.
        before = 5.0
        after = 30.0
        inPos = Time(inPosition.index[0])
        inPosM1 = inPos - TimeDelta(1.0, format='sec')
        tstart = inPos - TimeDelta(before, format='sec')
        tend = inPos + TimeDelta(after, format='sec')
        shutter = await client.select_time_series("lsst.sal.ATCamera.logevent_shutterDetailedState", "substate", tstart, tend)
        # These match within msec with the DATE-BEG and DATE-END timestamps in the header,
        # after we have converted DATE_END and DATE_BEG to UTC
        shutter_open = shutter.index[0]
        try:
            shutter_close = shutter.index[1]
        except:
            shutter_close = shutter.index[0] + pd.Timedelta(seconds=mData['EXPTIME'])

        # Now get the mount tracking info for a time before the inPosition timestamp,
        # and after the shutter closed time stamp
        before = 10.0
        after = 10.0
        inPos = Time(inPosition.index[0])
        tstart = inPos - TimeDelta(before, format='sec')
        tend = Time(shutter_close) + TimeDelta(after, format='sec')
        az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "azimuthCalculatedAngle",  tstart, tend)
        az_target = await client.select_time_series("lsst.sal.ATMCS.logevent_target", "azimuth",  tstart, tend)
        el = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "elevationCalculatedAngle",  tstart, tend)
        el_target = await client.select_time_series("lsst.sal.ATMCS.logevent_target", "elevation",  tstart, tend)

        if doOffset:
            offset = (tstart.jd - az.index[0].to_julian_date()) * 86400.0
            initial_offset = offset
            az.index += pd.DateOffset(seconds=offset)
            el.index += pd.DateOffset(seconds=offset)
            az_target.index += pd.DateOffset(seconds=offset)
            el_target.index += pd.DateOffset(seconds=offset)
        else:
            initial_offset = 0.0
    except:
        continue
    # Now plot it
    az_target_vals = np.array(az_target.values.tolist())[:,0]
    az_target_plus = az_target_vals + 0.0005
    az_target_minus = az_target_vals - 0.0005
    az_target_times = np.array(az_target.index.tolist())
    el_target_vals = np.array(el_target.values.tolist())[:,0]
    el_target_plus = el_target_vals + 0.0005
    el_target_minus = el_target_vals - 0.0005
    el_target_times = np.array(el_target.index.tolist())

    plotstart = (inPos - TimeDelta(5.0, format='sec')).to_datetime()
    plotend = (inPos + TimeDelta(15.0, format='sec')).to_datetime()

    
    fig = plt.figure(figsize=(8,8))
    plt.subplots_adjust(wspace=0.5)
    plt.suptitle(f"Mount Tracking - ExpId {expId}", fontsize = 16)
    # Azimuth axis
    plt.subplot(2,2,1)
    ax1 = az['azimuthCalculatedAngle'].plot(color='red')
    ax1.plot(az_target_times, az_target_vals, label='azimuth target', color='blue')
    ax1.plot(az_target_times, az_target_plus, ls = '--', color='blue')
    ax1.plot(az_target_times, az_target_minus, ls = '--', color='blue')
    ax1.set_title("Azimuth axis\n", fontsize=12)
    ax1.axvline(inPos.to_datetime(), color="green", linestyle="--", label="All Axes In Position")
    ax1.axvline(inPosM1.to_datetime(), color="green", linestyle="-.", label="In Position - 1.0")
    ax1.axvline(shutter_open, color='cyan', linestyle="--", label="Shutter_Open")
    ax1.axvline(shutter_close, color='cyan', linestyle="--", label="Shutter_Closed")
    ax1.set_ylabel("Degrees")
    #ax1.legend(loc='lower left')
    # Elevation axis
    plt.subplot(2,2,2)
    ax2 = el['elevationCalculatedAngle'].plot(label="Calculated Angle", color='red')
    ax2.plot(el_target_times, el_target_vals, label='Target', color='blue')
    ax2.plot(el_target_times, el_target_plus, ls='--', color='blue')
    ax2.plot(el_target_times, el_target_minus, ls='--', color='blue')
    ax2.set_title("Elevation axis\n", fontsize=12)
    ax2.axvline(inPos.to_datetime(), color="green", linestyle="--", label="All Axes In Position")
    ax2.axvline(inPosM1.to_datetime(), color="green", linestyle="-.", label="In Position - 1.0")
    ax2.axvline(shutter_open, color='cyan', linestyle="--", label="Shutter_Open")
    ax2.axvline(shutter_close, color='cyan', linestyle="--", label="Shutter_Closed")
    ax2.set_ylabel("Degrees")
    ax2.legend(bbox_to_anchor=(-0.5, -0.2))
    ax2.text(0.0,-0.5, f"Time offset = {initial_offset:.2f} seconds", transform=ax2.transAxes)

    ax1.set_xlim(plotstart, plotend)
    ax2.set_xlim(plotstart, plotend)


    
    pdf.savefig(fig)  # saves the current figure into a pdf page
    plt.close()
pdf.close()

```

```python

```
