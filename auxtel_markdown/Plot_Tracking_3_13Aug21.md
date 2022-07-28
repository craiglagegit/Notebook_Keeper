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

# AuxTel Mount fails - 13-Aug-21

In this notebook, investigate again mount tracking on 04-Aug-21\
With the steps outlined below, this seems to work for now.\

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
# Note that the DATE keyword is ~ 32 seconds before DATE-BEG and DATE-END
# This is because the DATE keyword is in UTC, as specified by the FITS spec.
# The DATE-BEG and DATE-END keywords are in TAI, as specified by the TIMESYS keyword

expId = 2021080400011
butler = Butler('/repo/main', collections="LATISS/raw/all")

mData = butler.get('raw.metadata', detector=0, exposure=expId)
print(f"{expId} \t {mData['TIMESYS']} \t {mData['DATE']} \t {mData['DATE-BEG']} \t {mData['DATE-END']}")
```

```python
# Use these for finding the "allAxesInPosition" timestamp
# The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
before = 20.0
after = 20.0
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
# The result says that this timestamp is in UTC.
# This is NOT CORRECT, and stems from Pandas lack of knowledge of TAI
# These timestamps are actually in TAI
inPosition.index[0]
```

```python
# Use these for finding the shutter status timestamp
# The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
# They agree within a few milliseconds.

timestamp = f"time >= '{start}+00:00' AND time <= '{end}+00:00'"
query = f'SELECT "substate" FROM "efd"."autogen"."lsst.sal.ATCamera.logevent_shutterDetailedState"\
    WHERE {timestamp}'

shutter = await client.influx_client.query(query)

# These match within msec with the DATE-BEG and DATE-END timestamps in the header
print(shutter.index[0], mData['DATE-BEG'])
print(shutter.index[1], mData['DATE-END'])
```

```python
# Now get the mount tracking info for a time before and after the inPosition timestamp.
# We need to tell it that these timestamps are in TAI.
before = 20.0
after = 20.0
inPos = Time(inPosition.index[0], scale='tai')
tstart = inPos - TimeDelta(before, format='sec')
tend = inPos + TimeDelta(after, format='sec')
print(f"{inPos} \t {tstart} \t {tend}")
```

```python
# Get and plot the data
# I need to override the merge_packed_time_series internal_time_scale to 'utc' in order for it all to work.
# As I understand it, this is a bug in astropy that is being worked.
mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                          tstart, tend)
nasmyth_position = await client.select_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", ['*'],
                                          tstart, tend)

az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1, internal_time_scale="utc")

# Plot it
fig = plt.figure(figsize = (16,6))
plt.suptitle(f"Mount Tracking - ExpId {expId}", fontsize = 18)
# Azimuth axis
plt.subplot(1,1,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.set_title("Azimuth axis", fontsize=16)
ax1.axvline(inPos.isot, color="green", linestyle="--", label="All Axes In Position")
ax1.axvline(Time(mData['DATE-BEG'], scale='utc').isot, color='cyan', linestyle="--", label="Exp_Start")
ax1.axvline(Time(mData['DATE-END'], scale='utc').isot, color='magenta', linestyle="--", label="Exp_End")
ax1.set_ylabel("Degrees")
ax1.legend()
plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{expId}_13Aug21.pdf")

```

```python

```
