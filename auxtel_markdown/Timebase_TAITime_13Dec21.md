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
expIds = [2021090900534, 2021101400007, 2021101400015, 2021110200563, 2021110400123, 2021113000005, 2021120200004]
print("expId\t\t\t Index Time  - taiTime (seconds)")
for expId in expIds:
    mData = butler.get('raw.metadata', detector=0, exposure=expId)
    #print(f"{expId} \t {mData['TIMESYS']} \t {mData['DATE']} \t {mData['DATE-BEG']} \t {mData['DATE-END']}")
    # Need to convert DATE_BEG and DATE_END to UTC to sync up with the EFD
    date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
    date_end = Time(mData['DATE-END'], format='isot', scale='tai')
    before = 30.0
    after = 30.0
    if expId < 2021101300000:
        # EFD was switched to UTC on 20211013
        # tai_offset compensates for this
        tai_offset = 37.0
    else:
        tai_offset = 0.0

    start = date_beg.utc - TimeDelta(before, format='sec') + TimeDelta(tai_offset, format='sec')
    end = date_end.utc + TimeDelta(after, format='sec') + TimeDelta(tai_offset, format='sec')
    #print(start, end)

    az_target = await client.select_time_series("lsst.sal.ATMCS.logevent_target",["azimuth", "taiTime"],  start, end)
    tIndex = Time(az_target.index[0]) - TimeDelta(tai_offset, format='sec')
    tTai = Time(az_target['taiTime'][0], format='unix_tai')
    dt = tIndex - tTai
    # The az shift and el shift are the same.
    #print(tIndex_az.isot, (tTai - TimeDelta(tai_shift, format='sec')).isot, dt_az.sec)
    print(f"{expId}\t\t {dt.sec}")
```

```python
az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ["azimuthCalculatedAngle", "taiTime"],  tstart, tend)
```

```python
print(date_beg.utc, date_end.utc)
print(date_beg.tai, date_end.tai)
```

```python
# Use these for finding the "allAxesInPosition" timestamp
# The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
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
az_target = await client.select_time_series("lsst.sal.ATMCS.logevent_target",["azimuth", "taiTime"],  tstart, tend)
```

```python
az_target
```

```python
test = pd.to_datetime(az_target['taiTime'][0])
print(test)
```

```python
test
```

```python
test = Time(az_target['taiTime'][0], format='unix_tai')
print(test.format)
print(test)
print(test.isot)
```

```python
tIndex = az_target.index[0]
tTai = Time(az_target['taiTime'][0], format='unix_tai').isot
print(tIndex, tTai)
```

```python
tIndex = Time(az_target.index[0])
tTai = Time(az_target['taiTime'][0], format='unix_tai')
print(tIndex.format, tTai.format)
dt = tIndex - tTai
print(tIndex.isot, tTai.isot, dt.sec)
```

```python
az_target
```

```python
az_target_vals = Time(np.array(az_target.values.tolist())[:,1], format='unix_tai')
print(az_target_vals[0:3].isot)
```

```python
pd.Timestamp(az_target_times[0].isot)
```

```python

```

```python
el_target_times[0]
```

```python
type(az_target_times[0])
```

```python
az_target_vals[0]
```

```python
az
```

```python
# Get the header data
#expId = 2021101400007
#expId = 2021110200563
expId = 2021120200004
mData = butler.get('raw.metadata', detector=0, exposure=expId)
#print(f"{expId} \t {mData['TIMESYS']} \t {mData['DATE']} \t {mData['DATE-BEG']} \t {mData['DATE-END']}")
# Need to convert DATE_BEG and DATE_END to UTC to sync up with the EFD
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
date_end = Time(mData['DATE-END'], format='isot', scale='tai')

# Get the inPosition timestamp

tai_shift = 37.0

shift = 37.0

before = 150.0
after = 30.0
if expId < 2021101300000:
    # EFD was switched to UTC on 20211013
    # tai_offset compensates for this
    tai_offset = 37.0
else:
    tai_offset = 0.0

start = date_beg.utc - TimeDelta(before, format='sec') + TimeDelta(tai_offset, format='sec')
end = date_end.utc + TimeDelta(after, format='sec') + TimeDelta(tai_offset, format='sec')

inPosition = await client.select_time_series("lsst.sal.ATMCS.logevent_allAxesInPosition", "inPosition", start, end)
inPosition = inPosition[inPosition['inPosition']==True] 

# Now get the mount tracking info for a time before and after the inPosition timestamp.

before = 120.0
after = 30.0
inPos = Time(inPosition.index[0])
tstart = inPos - TimeDelta(before, format='sec')
tend = inPos + TimeDelta(after, format='sec')

print(f"{inPos} \t {tstart} \t {tend}")
az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "azimuthCalculatedAngle",  tstart, tend)
az_target = await client.select_time_series("lsst.sal.ATMCS.logevent_target", ["azimuth", "taiTime"],  tstart, tend)
el = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "elevationCalculatedAngle",  tstart, tend)
el_target = await client.select_time_series("lsst.sal.ATMCS.logevent_target", ["elevation", "taiTime"],  tstart, tend)

print(f"Tstart={tstart}, Start of dataFrame = {az.index[0]}, {az_target.index[0]}")

# Plot it
az_target_vals = np.array(az_target.values.tolist())[:,0]
az_target_taiTimes = Time(np.array(az_target.values.tolist())[:,1]-tai_shift-shift, format='unix_tai')
el_target_vals = np.array(el_target.values.tolist())[:,0]
el_target_taiTimes = Time(np.array(el_target.values.tolist())[:,1]-tai_shift-shift, format='unix_tai')

az_target_times = []
el_target_times = []

for i in range(len(az_target_taiTimes)):
    az_target_times.append(pd.Timestamp(az_target_taiTimes[i].isot))
for i in range(len(el_target_taiTimes)):    
    el_target_times.append(pd.Timestamp(el_target_taiTimes[i].isot))

fig = plt.figure(figsize = (12,6))
plt.suptitle(f"Mount Tracking - ExpId {expId}", fontsize = 18)
# Azimuth axis
plt.subplot(1,2,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.plot(az_target_times, az_target_vals, label='azimuth target', color='blue')
ax1.set_title("Azimuth axis", fontsize=16)
ax1.axvline(inPosition.index[0], color="green", linestyle="--", label="All Axes In Position")
#ax1.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
#ax1.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax1.set_ylabel("Degrees")
ax1.legend()
# Elevation axis
plt.subplot(1,2,2)
ax2 = el['elevationCalculatedAngle'].plot(legend=True, color='red')
ax2.plot(el_target_times, el_target_vals, label='elevation target', color='blue')
ax2.set_title("Elevation axis", fontsize=16)
ax2.axvline(inPosition.index[0], color="green", linestyle="--", label="All Axes In Position")
#ax2.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
#ax2.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax2.set_ylabel("Degrees")
ax2.legend()

plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_taiTime_{expId}_13Dec21.pdf")

```

```python tags=[]
# Add a time offset to the packed series.
offset = 0.0# seconds

az_vals = np.array(az.values.tolist())[:,0]
az_times = np.array(az.index.tolist())
el_vals = np.array(el.values.tolist())[:,0]
el_times = np.array(el.index.tolist())

for i, t in enumerate(az_times):
    az_times[i] = az_times[i] + pd.Timedelta(value=offset*1.0E9) # pandas Timedelta is in nsec
    el_times[i] = el_times[i] + pd.Timedelta(value=offset*1.0E9) # pandas Timedelta is in nsec
az_target_vals = np.array(az_target.values.tolist())[:,0]
az_target_times = np.array(az_target.index.tolist())
el_target_vals = np.array(el_target.values.tolist())[:,0]
el_target_times = np.array(el_target.index.tolist())

# Interpolate target values to match packed values
az_target_interp = np.zeros([len(az_times)])
el_target_interp = np.zeros([len(el_times)])
lastj = 0
for i,azTime in enumerate(az_times):
    if (Time(azTime) - Time(az_target_times[0])) < 0:
        az_target_interp[i] = -1.0
        el_target_interp[i] = -1.0
    else:
        for j in range(lastj, len(az_target_times)):
            dT = Time(azTime) - Time(az_target_times[j])
            dT = dT.sec
            #print(j, azTime, azTTime, dT)
            if dT < 0 :
                break
        if j == 0:
            az_target_interp[i] = -1.0
            el_target_interp[i] = -1.0
            continue
        dT1 = Time(az_target_times[j]) - Time(az_target_times[j-1])
        dT1 = dT1.sec
        slope_az = (az_target_vals[j] - az_target_vals[j-1]) / dT1
        slope_el = (el_target_vals[j] - el_target_vals[j-1]) / dT1
        dT2 = Time(azTime) - Time(az_target_times[j-1])
        dT2 = dT2.sec
        az_target_interp[i] = az_target_vals[j-1] + slope_az * dT2
        el_target_interp[i] = el_target_vals[j-1] + slope_el * dT2
        lastj = j - 1
        if i%500 == 0:
            print(i, j, az_target_vals[j-1], az_target_vals[j], slope, az_target_interp[i])

```

```python
# Plot it
az_target_vals = np.array(az_target.values.tolist())[:,0]
az_target_times = np.array(az_target.index.tolist())
el_target_vals = np.array(el_target.values.tolist())[:,0]
el_target_times = np.array(el_target.index.tolist())

fig = plt.figure(figsize = (8,8))
plt.subplots_adjust(hspace=0.5, wspace=0.5)
plt.suptitle(f"Mount Tracking - ExpId {expId}", fontsize = 18)
# Azimuth axis
plt.subplot(2,2,1)
ax1 = az_target['azimuth'].plot(legend=True, label='Azimuth target', color='blue')
ax1.plot(az_times, az_vals, label='azimuth', color='red')
ax1.set_title("Azimuth axis", fontsize=16)
ax1.axvline(inPosition.index[0], color="green", linestyle="--", label="All Axes In Position")
#ax1.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
#ax1.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax1.set_ylabel("Degrees")
ax1.legend()
# Elevation axis
plt.subplot(2,2,2)
ax2 = el_target['elevation'].plot(legend=True, label='Elevation target', color='blue')
ax2.plot(el_times, el_vals, label='elevation', color='red')
ax2.set_title("Elevation axis", fontsize=16)
ax2.axvline(inPosition.index[0], color="green", linestyle="--", label="All Axes In Position")
#ax2.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
#ax2.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax2.set_ylabel("Degrees")
ax2.legend()
plt.subplot(2,2,3)
ax3 = az_target['azimuth'].plot(legend=False)
ax3.plot(az_times, abs(az_vals - az_target_interp), color='magenta', label='azimuthError')
ax3.set_title("Azimuth Error", fontsize=16)
ax3.axvline(inPosition.index[0], color="green", linestyle="--", label="All Axes In Position")
#ax3.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
#ax3.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax3.set_ylabel("Degrees")
ax3.set_ylim(0.0, 0.4)
ax3.legend()
plt.subplot(2,2,4)
ax4 = el_target['elevation'].plot(legend=False)
ax4.plot(el_times, abs(el_vals - el_target_interp), color='magenta', label='elevationError')
ax4.set_title("Elevation Error", fontsize=16)
ax4.axvline(inPosition.index[0], color="green", linestyle="--", label="All Axes In Position")
#ax4.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
#ax4.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax4.set_ylabel("Degrees")
ax4.set_ylim(0.0, 0.2)
ax4.legend()

plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{expId}_Offset_0p0_09Dec21.pdf")

```

```python tags=[]
for i in range(15000, 15500):
    print(az_times[i], az_vals[i], az_target_interp[i])
```

```python

```
