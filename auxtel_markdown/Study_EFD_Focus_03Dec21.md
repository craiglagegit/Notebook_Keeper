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

# AuxTel Focus Study - 03-Dec-21

In this notebook, investigate focus settings and temp on 03-Dec-21

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
# Get EFD client
from lsst_efd_client import EfdClient
from lsst_efd_client import  __version__ as efdVersion
print(efdVersion)
client = EfdClient('ldf_stable_efd')
```

```python
# Get one header data using Gen3 butler
# This confirms that the DATE_BEG and DATE_END timestamps remain in TAI, as specified.
before = 2.0
after = 2.0
tai_offset = 37.0

expId = 2021100500297
butler = Butler('/repo/main', collections="LATISS/raw/all")

mData = butler.get('raw.metadata', detector=0, exposure=expId)
print(f"{expId} \t {mData['TIMESYS']} \t {mData['DATE']} \t {mData['DATE-BEG']} \t {mData['DATE-END']}")
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
start = date_beg.utc - TimeDelta(before, format='sec') + TimeDelta(tai_offset, format='sec')
```

```python
# Get one header data using Gen3 butler
# This confirms that the DATE_BEG and DATE_END timestamps remain in TAI, as specified.

expId = 2021090800161
butler = Butler('/repo/main', collections="LATISS/raw/all")

mData = butler.get('raw.metadata', detector=0, exposure=expId)
print(f"{expId} \t {mData['TIMESYS']} \t {mData['DATE']} \t {mData['DATE-BEG']} \t {mData['DATE-END']}")
date_end = Time(mData['DATE-END'], format='isot', scale='tai')
end = date_end.utc + TimeDelta(after, format='sec') + TimeDelta(tai_offset, format='sec')
```

```python
# Use these for finding the various values
shutter = await client.select_time_series("lsst.sal.ATCamera.logevent_shutterDetailedState", "substate", start, end)
#shut_open = shutter[shutter['substate']==2]
#shut_closed = shutter[shutter['substate']==1]

#print(shut_open)
#print(shut_closed)
print(shutter)
```

```python
# Use these for finding the various values
command_z = await client.select_time_series("lsst.sal.ATAOS.command_offset", "z", start, end)
print(command_z)
```

```python
# Use these for finding the various values
corr_off = await client.select_time_series("lsst.sal.ATAOS.logevent_correctionOffsets", "z", start, end)
print(corr_off)
```

```python
# Use these for finding the various values
total_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "total", start, end)
print(total_off)
```

```python
# Use these for finding the various values
disp_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "disperser", start, end)
print(disp_off)
```

```python
# Use these for finding the various values
filter_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "filter", start, end)
print(filter_off)
```

```python
# Use these for finding the various values
user_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "userApplied", start, end)
print(user_off)
```

```python
# Use these for finding the various values
wave_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "wavelength", start, end)
print(wave_off)
```

```python
# Use these for finding the various values
pr_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "priority", start, end)
print(pr_off)
```

```python
# Use these for finding the various values
z_position = await client.select_time_series("lsst.sal.ATHexapod.command_moveToPosition", "z", start, end)
print(z_position)
```

```python
# Use these for finding the various values
shutter = await client.select_time_series("lsst.sal.ATCamera.logevent_shutterDetailedState", "substate", start, end)
shut_open = shutter[shutter['substate']==2]
shut_closed = shutter[shutter['substate']==1]

print(shut_open)
print(shut_closed)
#print(shutter)
```

```python
# Plot it
fig = plt.figure(figsize = (8,6))
#plt.suptitle(f"Mount Tracking - ExpId {expId}", fontsize = 18)
# Azimuth axis
plt.subplot(1,1,1)
ax1 = total_off['total'].plot(legend=True, color='red', label = 'Offset')
ax2 = z_position['z'].plot(legend=True, color='blue', label = 'Z-position')
#ax1.set_title("Azimuth axis", fontsize=16)
for i in range(len(shut_open)):
    ax1.axvline(shut_open.index[i], color='cyan', linestyle="--", label="Exp_Start")
for i in range(len(shut_closed)):
    ax1.axvline(shut_closed.index[i], color='magenta', linestyle="--", label="Exp_End")

#ax1.set_ylabel("Degrees")
#ax1.legend()
#plt.savefig(f"/project/cslage/AuxTel/offsets/Tracking_Timebase_{expId}_29Oct21.pdf")

```

```python
print(len(z_position), len(total_off))
```

```python
for i in range(len(z_position)):
    print(z_position['z'][i], total_off['total'][i+2], z_position['z'][i] - total_off['total'][i+2])
```

```python
# Offset is consistent with z_position.  Difference is constant
print(z_position['z'][-1], total_off['total'][-1], z_position['z'][-1] - total_off['total'][-1])
print(z_position['z'][-2], total_off['total'][-2], z_position['z'][-2] - total_off['total'][-2])
print(z_position['z'][-3], total_off['total'][-3], z_position['z'][-3] - total_off['total'][-3])
print(z_position['z'][-4], total_off['total'][-6], z_position['z'][-4] - total_off['total'][-6])
print(z_position['z'][-5], total_off['total'][-9], z_position['z'][-5] - total_off['total'][-9])

```

```python

```

```python
# Use these for finding the various values
temp_air = await client.select_time_series("lsst.sal.ESS.temperature4Ch", "temperatureC02", start, end)
temp_truss = await client.select_time_series("lsst.sal.ESS.temperature4Ch", "temperatureC03", start, end)
temp_m2 = await client.select_time_series("lsst.sal.ESS.temperature4Ch", "temperatureC04", start, end)
print(temp_air.tail(1))
print(temp_truss.tail(1))
print(temp_m2.tail(1))
```

```python
# Use these for finding the various values
temp_ext = await client.select_time_series("lsst.sal.WeatherStation.airTemperature", "avg1M", start, end)
print(temp_ext)
```

```python

```

```python
# Get one header data using Gen3 butler
# This confirms that the DATE_BEG and DATE_END timestamps remain in TAI, as specified.
before = 15.0
after = 10.0
tai_offset = 37.0

expId = 2021100500297
mData = butler.get('raw.metadata', detector=0, exposure=expId)
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
start = date_beg.utc - TimeDelta(before, format='sec') + TimeDelta(tai_offset, format='sec')
end = date_beg.utc + TimeDelta(after, format='sec') + TimeDelta(tai_offset, format='sec')
print(date_beg, start, end)
```

```python
# Use these for finding the various values
shutter = await client.select_time_series("lsst.sal.ATCamera.logevent_shutterDetailedState", "substate", start, end)
shut_open = shutter[shutter['substate']==2]
shut_closed = shutter[shutter['substate']==1]

print(shut_open)
print(shut_closed)
#print(shutter)
```

```python
# Use these for finding the various values
total_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "total", start, end)
print(total_off)
```

```python
# Use these for finding the various values
disp_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "disperser", start, end)
print(disp_off)
```

```python
# Use these for finding the various values
filter_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "filter", start, end)
print(filter_off)
```

```python
# Use these for finding the various values
user_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "userApplied", start, end)
print(user_off)
```

```python
# Use these for finding the various values
wave_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "wavelength", start, end)
print(wave_off)
```

```python
# Use these for finding the various values
pr_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "priority", start, end)
print(pr_off)
```

```python

```
