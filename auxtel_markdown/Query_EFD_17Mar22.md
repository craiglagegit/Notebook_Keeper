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

## AuxTel Image elongation due to azimuth oscillation

Craig Lage 17-Mar-22

```python
import sys, time, os, asyncio, glob

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pickle as pkl
import pandas as pd
import astropy.io.fits as pf
from astropy.time import Time, TimeDelta
from astropy.coordinates import AltAz, ICRS, EarthLocation, Angle, FK5
import astropy.units as u

from lsst.daf.butler import Butler
from lsst_efd_client import EfdClient

from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python
client = EfdClient('ldf_stable_efd')
```

```python tags=[]
dir(client)
```

```python tags=[]
await client.get_topics()
```

```python

```

```python
await client.get_fields('lsst.sal.ATMCS.logevent_azimuthDrive1Status')
```

```python
# Times to start looking at error codes
start = Time("2022-03-17 00:00:00Z", scale='utc')
end = Time("2022-03-17 09:00:00Z", scale='utc')
```

```python
errors = await client.select_time_series('lsst.sal.ATMCS.logevent_errorCode', \
                                                ['*'],  start, end)
```

```python tags=[]
errors
```

```python
# Look before and after the last error time
lastErrorTime = Time(errors.index[3])
before = 300.0
after = 300.0
start = lastErrorTime - TimeDelta(before, format='sec')
end = lastErrorTime + TimeDelta(after, format='sec')
print(start, end)
```

```python
status = await client.select_time_series('lsst.sal.ATMCS.logevent_azimuthDrive1Status', \
                                                ['*'],  start, end)
```

```python
len(status)
```

```python
status
```

```python
# Times to start looking at error codes
start = Time("2022-03-17 00:00:00Z", scale='utc')
end = Time("2022-03-17 00:05:00Z", scale='utc')
test = await client.select_packed_time_series('lsst.sal.ATPtg.mountPositions', \
                                                ['*'],  start, end)
```

```python

```
