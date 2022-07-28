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

## ComCam status

Query ComCam camera status from NCSA\
Craig Lage - 28-May-21

```python
import sys, time, os, asyncio, glob
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
from lsst_efd_client import EfdClient
```

```python
# Get EFD client
client = EfdClient('summit_efd')
```

```python
t_end = Time(time.time(),format='unix', scale='tai')
t_start = t_end - TimeDelta(864000, format='sec') # Get transitions for last 10 days
archiver = await client.select_time_series("lsst.sal.CCArchiver.logevent_summaryState", ['*'],
                                      t_start, t_end)
camera = await client.select_time_series("lsst.sal.CCCamera.logevent_summaryState", ['*'],
                                      t_start, t_end)
headerService = await client.select_time_series("lsst.sal.CCHeaderService.logevent_summaryState", ['*'],
                                      t_start, t_end)

filter = await client.select_time_series("lsst.sal.CCCamera.logevent_endSetFilter", ['*'],
                                      t_start, t_end)

for [sal, name] in [[archiver, "CCArchiver"], [camera, "CCCamera"], [headerService, "CCHeaderService"]]:
    summaryState = sal['summaryState'][-1]
    if summaryState == 1:
        print(name + " is in state DISABLED")
    elif summaryState == 2:
        print(name + " is in state ENABLED")
    elif summaryState == 3:
        print(name + " is in state FAULT")
    elif summaryState == 4:
        print(name + " is in state OFFLINE")
    elif summaryState == 5:
        print(name + " is in state STANDBY")

print("Current filter is " + filter['filterName'][-1])
```

```python

```
