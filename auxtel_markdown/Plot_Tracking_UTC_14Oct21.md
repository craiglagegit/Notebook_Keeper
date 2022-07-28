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

# AuxTel Mount fails - 14-Oct-21

In this notebook, investigate again mount tracking on 14-Oct-21\
This is after the EFD was converted to UTC.

```python
import sys, time, os, asyncio

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.time import Time, TimeDelta
from lsst.daf.butler import Butler
```

```python
from lsst_efd_client import EfdClient, merge_packed_time_series, __version__
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

expId = 2021101400011
butler = Butler('/repo/main', collections="LATISS/raw/all")

mData = butler.get('raw.metadata', detector=0, exposure=expId)
print(f"{expId} \t {mData['TIMESYS']} \t {mData['DATE']} \t {mData['DATE-BEG']} \t {mData['DATE-END']}")
```

```python
# Need to convert DATE_BEG and DATE_END to UTC to sync up with the EFD
tai_delta = 37.0
date_beg_utc = Time(mData['DATE-BEG'],format='isot', scale='utc') - TimeDelta(tai_delta, format='sec')
date_end_utc = Time(mData['DATE-END'],format='isot', scale='utc') - TimeDelta(tai_delta, format='sec')
print(date_beg_utc, date_end_utc)
print(Time(mData['DATE-BEG'],format='isot', scale='tai'), Time(mData['DATE-END'],format='isot', scale='tai'))
```

```python
# Use these for finding the "allAxesInPosition" timestamp
# The inPosition timestamp makes sense with the DATE-BEG and DATE-END times
before = 5.0
after = 5.0
start = date_beg_utc - TimeDelta(before, format='sec')
end = date_end_utc + TimeDelta(after, format='sec')
print(start, end)
timestamp = f"time >= '{start}+00:00' AND time <= '{end}+00:00'"
query = f'SELECT "inPosition" FROM "efd"."autogen"."lsst.sal.ATMCS.logevent_allAxesInPosition"\
    WHERE {timestamp} and inPosition = true'

inPosition = await client.influx_client.query(query)
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

timestamp = f"time >= '{start}+00:00' AND time <= '{end}+00:00'"
query = f'SELECT "substate" FROM "efd"."autogen"."lsst.sal.ATCamera.logevent_shutterDetailedState"\
    WHERE {timestamp}'

shutter = await client.influx_client.query(query)

# These match within msec with the DATE-BEG and DATE-END timestamps in the header,
# after we have converted DATE_END and DATE_BEG to UTC
print(shutter.index[0], date_beg_utc)
print(shutter.index[1], date_end_utc)
```

```python
# Now get the mount tracking info for a time before and after the inPosition timestamp.
# We need to tell it that these timestamps are in TAI.
before = 5.0
after = 5.0
inPos = Time(inPosition.index[0])
tstart = inPos - TimeDelta(before, format='sec')
tend = inPos + TimeDelta(after, format='sec')
print(f"{inPos} \t {tstart} \t {tend}")
```

```python
# Get and plot the data
# It appears tha the packed time series are still in TAI
mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                          tstart, tend)

az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1)

# Plot it
fig = plt.figure(figsize = (16,6))
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
plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{expId}_14Oct21.pdf")

```

```python
# Get and plot the data
# I need to override the merge_packed_time_series internal_time_scale to 'utc' in order for it all to work.
# As I understand it, this is a bug in astropy that is being worked.
mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                          tstart - TimeDelta(tai_delta, format='sec'), tend - TimeDelta(tai_delta, format='sec'))

az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1)

# Plot it
fig = plt.figure(figsize = (16,6))
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
plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{expId}_TAI_Shift_14Oct21.pdf")

```

```python

```
