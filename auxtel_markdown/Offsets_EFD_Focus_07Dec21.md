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

# AuxTel Focus Study - 07-Dec-21

In this notebook, investigate focus settings and temp

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

```python tags=[]
# Gen3 butler
from lsst.daf.butler import Butler
dayObs = '2021-10-05'
dayObs = int(dayObs.replace('-', ''))
butler = Butler('/repo/main', collections="LATISS/raw/all")

exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
    exposureList.append([record.id, record])
exposureList.sort(key=lambda x: x[0])
for [id,record] in exposureList:
    if id > 2021100500292 and id < 2021100500298:
        print(record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name)

```

```python
# Get the header data
# The DATE_BEG and DATE_END timestamps remain in TAI, as specified.
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
