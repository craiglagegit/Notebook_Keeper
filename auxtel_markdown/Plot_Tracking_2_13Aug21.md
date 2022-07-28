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

## AuxTel Mount fails - 13-Aug-21

In this notebook, investigate again mount tracking on 25-May-21\
I can get things to line up, by repeatedly telling the code that UTC times are really TAI. \
This is all a big mess, but this seems to work for now.\
Modifying to use a single image.

```python
import sys, time, os, asyncio

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.time import Time, TimeDelta
from lsst_efd_client import EfdClient, merge_packed_time_series
from lsst.daf.butler import Butler
```

```python
# Get EFD client
client = EfdClient('ldf_stable_efd')
```

```python
# Get the header data using Gen3 butler
# These two exposures are of two different objects, so we know we slewed between them
# Note that the DATE keyword is ~ 30 seconds before DATE-BEG and DATE-END
# This doesn't seem right

expId = 2021052500184
butler = Butler('/repo/main', collections="LATISS/raw/all")

mData = butler.get('raw.metadata', detector=0, exposure=expId)
print(f"{expId} \t {mData['OBJECT']} \t {mData['DATE']} \t {mData['DATE-BEG']} \t {mData['DATE-END']}")
```

```python
# Use these for finding the "allAxesInPosition" timestamp
# The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
# But are these times in UTC, or TAI?
before = 60.0
after = 5.0
start = Time(mData['DATE-BEG'],format='isot', scale='utc') - TimeDelta(before, format='sec')
end = Time(mData['DATE-END'],format='isot', scale='utc') + TimeDelta(after, format='sec')
print(start, end)
timestamp = f"time >= '{start}+00:00' AND time <= '{end}+00:00'"
query = f'SELECT "inPosition" FROM "efd"."autogen"."lsst.sal.ATMCS.logevent_allAxesInPosition"\
    WHERE {timestamp} and inPosition = true'

inPosition = await client.influx_client.query(query)
print(inPosition)
```

```python
# Why is the timestamp in UTC?
inPosition.index[0]
```

```python
# Use these for finding the shutter status timestamp
# The inPosition timestamp makes sense with the DATE-BEG and DATE-END times

timestamp = f"time >= '{start}+00:00' AND time <= '{end}+00:00'"
query = f'SELECT "substate" FROM "efd"."autogen"."lsst.sal.ATCamera.logevent_shutterDetailedState"\
    WHERE {timestamp}'

shutter = await client.influx_client.query(query)
#print(inPosition)
```

```python
# These match perfectly with the DATE-BEG and DATE-END timestamps
print(shutter)
```

```python
# Now get the mount tracking info for a time before and after the inPosition timestamp.
before = 60.0
after = 60.0
inPos = Time(inPosition.index[0], scale='tai') # We lie to it and tell it it is TAI.
tstart = inPos - TimeDelta(before, format='sec')
tend = inPos + TimeDelta(after, format='sec')
print(f"{inPos} \t {tstart} \t {tend}")
```

```python
# Get and plot the data
# Note that when it gets the data, it adds another 37 seconds to tstart and tend!!!
# If I change merge_packed_time_series internal_time_scale to 'utc', then it doesn't do this.
mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                          tstart, tend)
nasmyth_position = await client.select_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", ['*'],
                                          tstart, tend)

az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1, internal_time_scale="utc")
el = merge_packed_time_series(mount_position, 'elevationCalculatedAngle', stride=1, internal_time_scale="utc")
rot = merge_packed_time_series(nasmyth_position, 'nasmyth2CalculatedAngle', stride=1, internal_time_scale="utc")

# Plot it
fig = plt.figure(figsize = (16,6))
plt.suptitle(f"Mount Tracking - ExpId {expId}", fontsize = 18)
# Azimuth axis
plt.subplot(1,1,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.set_title("Azimuth axis", fontsize=16)
ax1.axvline(inPos.isot, color="green", linestyle="--", label="In Position")
#ax1.axvline(Time(mDatas[expIds[0]]['DATE-BEG'], scale='utc').isot, color='blue', linestyle="--", label="Exp1_Start")
#ax1.axvline(Time(mDatas[expIds[0]]['DATE-END'], scale='utc').isot, color='green', linestyle="--", label="Exp1_End")
ax1.axvline(Time(mData['DATE-BEG'], scale='utc').isot, color='cyan', linestyle="--", label="Exp2_Start")
ax1.axvline(Time(mData['DATE-END'], scale='utc').isot, color='magenta', linestyle="--", label="Exp2_End")
ax1.set_ylabel("Degrees")
ax1.legend()
#plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{expIds[1]}_25May21.pdf")

```

```python

```
