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

```python tags=[]
expStart = 2022021500067
expId = expStart
mData = butler.get('raw.metadata', detector=0, exposure=expId)
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
date_end = Time(mData['DATE-END'], format='isot', scale='tai')
tstart = date_beg.utc
expEnd = 2022021500580
expId = expEnd
mData = butler.get('raw.metadata', detector=0, exposure=expId)
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
date_end = Time(mData['DATE-END'], format='isot', scale='tai')
tstop = date_end.utc    
```

```python
torque = await client.select_packed_time_series("lsst.sal.ATMCS.measuredTorque", 'nasmyth2MotorTorque',
                                              tstart, tstop)
```

```python
angle = await client.select_packed_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", 'nasmyth2CalculatedAngle',
                                              tstart, tstop)
```

```python
torqueList = torque.values.tolist()[::100]
angleList = angle.values.tolist()[::100]
plt.figure()
plt.scatter(np.array(angleList)[:,0],np.array(torqueList)[:,0])
plt.plot([-160,160],[3.0,3.0], color='red', ls='--')
plt.plot([-160,160],[-3.0,-3.0], color='red', ls='--')
plt.arrow(-140, 2.5, 50,0, width=0.1,head_length = 5.0, color='green')
plt.arrow(140, -2.5, -50,0, width=0.1,head_length = 5.0, color='green')
plt.xlabel("Rotator angle(degrees)")
plt.ylabel("Torque (amps)")
plt.savefig("/project/cslage/AuxTel/mount_graphs/Torque_vs_Angle_Observing_15Feb22.pdf")
```

```python
len(torqueList)
```

```python

```
