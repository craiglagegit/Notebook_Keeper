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
# Get EFD client
client = EfdClient('ldf_stable_efd')
```

```python
# Get one header data using Gen3 butler
# This confirms that the DATE_BEG and DATE_END timestamps remain in TAI, as specified.

expId = 2021101400007
butler = Butler('/repo/main', collections="LATISS/raw/all")

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
before = 5.0
after = 5.0
start = date_beg.utc - TimeDelta(before, format='sec')
end = date_end.utc + TimeDelta(after, format='sec')
print(start, end)

inPosition = await client.select_time_series("lsst.sal.ATMCS.logevent_allAxesInPosition", "inPosition", start, end)
inPosition = inPosition[inPosition['inPosition']==True] 
print(inPosition)
```

```python
# The result says that this timestamp is in UTC.
# This is CORRECT
inPosition.index[0]
```

```python
# Use these for finding the shutter status timestamp
# The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
# They agree within a few milliseconds.

shutter = await client.select_time_series("lsst.sal.ATCamera.logevent_shutterDetailedState", "substate", start, end)

# These match within msec with the DATE-BEG and DATE-END timestamps in the header,
# after we have converted DATE_END and DATE_BEG to UTC
print(shutter.index[0], date_beg.utc)
print(shutter.index[1], date_end.utc)
```

```python
# Now get the mount tracking info for a time before and after the inPosition timestamp.
before = 5.0
after = 5.0
inPos = Time(inPosition.index[0])
tstart = inPos - TimeDelta(before, format='sec')
tend = inPos + TimeDelta(after, format='sec')
print(f"{inPos} \t {tstart} \t {tend}")
az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "azimuthCalculatedAngle",  tstart, tend)
print(f"Tstart={tstart}, Start of dataFrame = {az.index[0]}")
```

```python
# Plot it
fig = plt.figure(figsize = (12,6))
plt.suptitle(f"Mount Tracking - ExpId {expId}", fontsize = 18)
# Azimuth axis
plt.subplot(1,1,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.set_title("Azimuth axis", fontsize=16)
ax1.axvline(inPos.isot, color="green", linestyle="--", label="All Axes In Position")
ax1.axvline(shutter.index[0], color='cyan', linestyle="--", label="Exp_Start")
ax1.axvline(shutter.index[1], color='magenta', linestyle="--", label="Exp_End")
ax1.set_ylabel("Degrees")
ax1.legend()
#plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{expId}_29Oct21.pdf")

```

```python
# Now get and plot similar data for the 10 exposures taken

firstExpId = 2021101400006
mData = butler.get('raw.metadata', detector=0, exposure=firstExpId)
date_beg = Time(mData['DATE-BEG'],format='isot', scale='tai')
lastExpId = 2021101400017
mData = butler.get('raw.metadata', detector=0, exposure=lastExpId)
date_end = Time(mData['DATE-END'],format='isot', scale='tai')

before = 5.0
after = 5.0
tstart = date_beg.utc - TimeDelta(before, format='sec')
tend = date_end.utc + TimeDelta(after, format='sec')

# Get the inPosition timestamps
inPosition = await client.select_time_series("lsst.sal.ATMCS.logevent_allAxesInPosition", "inPosition", tstart, tend)
inPosition = inPosition[inPosition['inPosition']==True] 
#print(inPosition)

# Get the shutter open and close timestamps

shutter = await client.select_time_series("lsst.sal.ATCamera.logevent_shutterDetailedState", "substate", tstart, tend)

open_shutter = shutter[shutter['substate']==2] 
#print(open_shutter)

close_shutter = shutter[shutter['substate']==1] 
#print(close_shutter)

# Get the mount tracking data

az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "azimuthCalculatedAngle",  tstart, tend)
print(f"Tstart={tstart}, Start of dataFrame = {az.index[0]}")

# Plot it
fig = plt.figure(figsize = (8,6))
plt.suptitle(f"Mount Tracking - ExpIds {firstExpId} - {lastExpId}", fontsize = 18)
# Azimuth axis
plt.subplot(1,1,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.set_title("Azimuth axis", fontsize=16)
for n in range(len(inPosition)):
    if n==0:
        ax1.axvline(inPosition.index[n], color="green", linestyle="--", label="All Axes In Position")
    else:
        ax1.axvline(inPosition.index[n], color="green", linestyle="--", label="")
for n in range(len(open_shutter)):
    if n==0:
        ax1.axvline(open_shutter.index[n], color='cyan', linestyle="--", label="Exp_Start")
    else:
        ax1.axvline(open_shutter.index[n], color='cyan', linestyle="--", label="")
for n in range(len(close_shutter)):
    if n==0:
        ax1.axvline(close_shutter.index[n], color='magenta', linestyle="--", label="Exp_End")
    else:
        ax1.axvline(close_shutter.index[n], color='magenta', linestyle="--", label="")
        
ax1.set_ylabel("Degrees")
ax1.legend()
#plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{firstExpId}_{lastExpId}_29Oct21.pdf")

```

```python

```
