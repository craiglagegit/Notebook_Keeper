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

# AuxTel mount issue - 25-Feb-2022

In this notebook, investigate mount issue from 20220215\
Why was the mount issued a "full stop"?

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

```python jupyter={"outputs_hidden": true} tags=[]
for i, expId in enumerate(range(2022021600088, 2022021600150)):
    mData = butler.get('raw.metadata', detector=0, exposure=expId)
    date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
    if i > 0:
        delta = TimeDelta(date_beg.utc - last_date_beg.utc)
        print(i, expId, delta.sec)
    last_date_beg = date_beg
```

```python tags=[]
expStart = 2022021600088
expEnd = 2022021600150
indices = []
begs = []
ends = []
bads = []
for i, expId in enumerate(range(expStart, expEnd)):
    index = expId - 2022021600000
    indices.append(index)
    mData = butler.get('raw.metadata', detector=0, exposure=expId)
    date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
    date_end = Time(mData['DATE-END'], format='isot', scale='tai')
    if i == 0:
        tstart = date_beg.utc
    begs.append(date_beg.utc)
    ends.append(date_end.utc)
    if index in [200,206,218,224,230,236,242,262]:
        bads.append(date_beg.utc)
tstop = date_end.utc    
```

```python
az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", 
                                            "azimuthCalculatedAngle",  tstart, tstop)
```

```python
stopTracking = await client.select_time_series("lsst.sal.ATMCS.command_stopTracking", 
                                            "value",  tstart, tstop)
```

```python
len(stopTracking)
```

```python
fig = plt.figure(figsize = (8,4))
plt.suptitle(f"Mount Tracking - ExpIds {expStart} - {expEnd}", fontsize = 18)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
for i in range(len(begs)):
    ax1.axvline(begs[i].to_datetime(), ymin=0.1, ymax=0.25, color="blue")
ax1.axvline(begs[i].to_datetime(), ymin=0.1, ymax=0.25, color="blue", label="Shutter open")    
for j in range(len(bads)):
    ax1.axvline(bads[j].to_datetime(), ymin=0.27, ymax=0.42, color="green")  
if len(bads) > 0:
    ax1.axvline(bads[j].to_datetime(), ymin=0.27, ymax=0.42, color="green", label="Bad tracking")      
for k in range(len(stopTracking)):
    ax1.axvline(stopTracking.index[k], ymin=0.44, ymax=0.59, color="magenta")
if len(stopTracking) > 0:    
    ax1.axvline(stopTracking.index[k], ymin=0.44, ymax=0.59, color="magenta", label="StopTracking")    
#ax1.set_ylim(179.5, 180.5)
ax1.legend(loc='upper right')
plt.savefig("/project/cslage/AuxTel/mount_graphs/Mount_Fails_16Feb22.pdf")
```

```python

```
