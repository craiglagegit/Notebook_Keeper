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

# AuxTel Plot tracking - 29-Oct-21

In this notebook, investigate again mount tracking on 29-Oct-21\
This is after the EFD was converted to UTC.\
Thanks to Simon Krughoff for contributions.

```python
import sys, time, os, asyncio

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
%matplotlib widget
import pandas as pd
from astropy.time import Time, TimeDelta
from lsst.daf.butler import Butler
```

```python
from lsst_efd_client import EfdClient
from lsst_efd_client import  __version__ as efdVersion
print(efdVersion)
```

```python
# Get EFD client and the butler
client = EfdClient('ldf_stable_efd')
butler = Butler('/repo/main', collections="LATISS/raw/all")
```

```python
# Get one header data using Gen3 butler
# This confirms that the DATE_BEG and DATE_END timestamps remain in TAI, as specified.

expId = 2021113000005
mData = butler.get('raw.metadata', detector=0, exposure=expId)
print(f"{expId} \t {mData['TIMESYS']} \t {mData['DATE']} \t {mData['DATE-BEG']} \t {mData['DATE-END']}")
```

```python
# Need to convert DATE_BEG and DATE_END to UTC to sync up with the EFD
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
date_end = Time(mData['DATE-END'], format='isot', scale='tai')
print(date_beg.utc, date_end.utc)
print(date_beg.tai, date_end.tai)
```

```python
# Use these for finding the "allAxesInPosition" timestamp
# The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
before = 150.0
after = 30.0
if expId < 2021101300000:
    # EFD was switched to UTC on 20211013
    tai_offset = 37.0
else:
    tai_offset = 0.0

start = date_beg.utc - TimeDelta(before, format='sec') + TimeDelta(tai_offset, format='sec')
end = date_end.utc + TimeDelta(after, format='sec') + TimeDelta(tai_offset, format='sec')
print(start, end)

inPosition = await client.select_time_series("lsst.sal.ATMCS.logevent_allAxesInPosition", "inPosition", start, end)
inPosition = inPosition[inPosition['inPosition']==True] 
print(inPosition)
```

```python
# Use these for finding the shutter status timestamp
# The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
# They agree within a few milliseconds.
before = 5.0
after = 150.0
inPos = Time(inPosition.index[0])
tstart = inPos - TimeDelta(before, format='sec')
tend = inPos + TimeDelta(after, format='sec')

shutter = await client.select_time_series("lsst.sal.ATCamera.logevent_shutterDetailedState", "substate", tstart, tend)

# These match within msec with the DATE-BEG and DATE-END timestamps in the header,
# after we have converted DATE_END and DATE_BEG to UTC
print(shutter.index[0], date_beg.utc)
print(shutter.index[1], date_end.utc)
```

```python
# Now get the mount tracking info for a time before and after the inPosition timestamp.

before = 150.0
after = 30.0
inPos = Time(inPosition.index[0])
tstart = inPos - TimeDelta(before, format='sec')
tend = inPos + TimeDelta(after, format='sec')
packed_tstart = tstart + TimeDelta(tai_offset, format='sec')
packed_tend = tend + TimeDelta(tai_offset, format='sec')

print(f"{inPos} \t {tstart} \t {tend}")
az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "azimuthCalculatedAngle",  packed_tstart, packed_tend)
az_target = await client.select_time_series("lsst.sal.ATMCS.logevent_target", "azimuth",  tstart, tend)
print(f"Tstart={tstart}, Start of dataFrame = {az.index[0]}, {az_target.index[0]}")
```

```python
# Plot it
az_target_vals = np.array(az_target.values.tolist())[:,0]
az_target_times = np.array(az_target.index.tolist())
fig = plt.figure(figsize = (8,6))
plt.suptitle(f"Mount Tracking - ExpId {expId}", fontsize = 18)
# Azimuth axis
plt.subplot(1,1,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.plot(az_target_times, az_target_vals, label='azimuth target', color='blue')
ax1.set_title("Azimuth axis", fontsize=16)
ax1.axvline(inPosition.index[0], color="green", linestyle="--", label="All Axes In Position")
ax1.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
ax1.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax1.set_ylabel("Degrees")
ax1.legend()
#plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{expId}_29Oct21.pdf")

```

```python tags=[]
# Add a time offset to the packed series.
offset = 24.75# seconds
az_vals = np.array(az.values.tolist())[:,0]
az_times = np.array(az.index.tolist())
for i, t in enumerate(az_times):
    az_times[i] = t + pd.Timedelta(value=offset*1.0E9) + pd.Timedelta(value=tai_offset*1.0E9)# pandas Timedelta is in nsec
az_target_vals = np.array(az_target.values.tolist())[:,0]
az_target_times = np.array(az_target.index.tolist())
# Interpolate target values to match packed values
az_target_interp = np.zeros([len(az_times)])
lastj = 0
for i,azTime in enumerate(az_times):
    if (Time(azTime) - Time(az_target_times[0])) < 0:
        az_target_interp[i] = -1.0
    else:
        for j in range(lastj, len(az_target_times)):
            dT = Time(azTime) - Time(az_target_times[j])
            dT = dT.sec
            #print(j, azTime, azTTime, dT)
            if dT < 0 :
                break
        if j == 0:
            az_target_interp[i] = -1.0
            continue
        dT1 = Time(az_target_times[j]) - Time(az_target_times[j-1])
        dT1 = dT1.sec
        slope = (az_target_vals[j] - az_target_vals[j-1]) / dT1
        dT2 = Time(azTime) - Time(az_target_times[j-1])
        dT2 = dT2.sec
        az_target_interp[i] = az_target_vals[j-1] + slope * dT2
        lastj = j - 1
        if i%500 == 0:
            print(i, j, az_target_vals[j-1], az_target_vals[j], slope, az_target_interp[i])

```

```python
# Plot it with the offset
fig = plt.figure(figsize = (8,6))
plt.subplots_adjust(hspace=0.5)
plt.suptitle(f"Mount Tracking - ExpId {expId} - Offset={offset} sec", fontsize = 16)
# Azimuth axis
plt.subplot(2,1,1)
ax1 = az_target['azimuth'].plot(legend=True, label='azimuth target', color='blue')
ax1.plot(az_times, az_vals, color='red', label='azimuthCalculatedAngle')
ax1.set_title("Azimuth axis", fontsize=16)
ax1.axvline(inPosition.index[0], color="green", linestyle="--", label="All Axes In Position")
ax1.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
ax1.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax1.set_ylabel("Degrees")
ax1.legend()
plt.subplot(2,1,2)
ax2 = az_target['azimuth'].plot(legend=False)
ax2.plot(az_times, abs(az_vals - az_target_interp), color='magenta', label='azimuthError')
ax2.set_title("Azimuth Error", fontsize=16)
ax2.axvline(inPosition.index[0], color="green", linestyle="--", label="All Axes In Position")
ax2.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
ax2.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax2.set_ylabel("Degrees")
ax2.set_ylim(0.0, 0.4)
ax2.legend()

plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{expId}_Offset_24p7_09Dec21.pdf")

```

```python
for i in range(15000, 15500):
    print(az_times[i], az_vals[i], az_target_interp[i])
```

```python

```
