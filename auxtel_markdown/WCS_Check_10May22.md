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

## WCS Check
Craig Lage 10-May-22

```python
import sys, time, os, asyncio, glob
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
from astropy.coordinates import AltAz, ICRS, EarthLocation, Angle, FK5, SkyCoord
import astropy.units as u
from lsst.daf.butler import Butler
```

```python
butler = Butler('/repo/main', collections="LATISS/raw/all")
```

```python tags=[]
%%capture --no-stdout
expIdList = [2022050500695, \
             2022040601000, 2022021700321, 2022040500400, 2022021700311]
print("ExpId          RA       Dec     ROTPA     WCS.x     WCS.y")
for expId in expIdList:
    exp = butler.get('raw', detector=0, exposure=expId)
    mData = exp.getMetadata()
    el = Angle(mData['ELSTART'] * u.deg)
    az = Angle(mData['AZSTART'] * u.deg)
    dec = Angle(mData['DECSTART'] * u.deg)
    ra = Angle(mData['RASTART'] * u.deg)
    rotpa = Angle(mData['ROTPA']*u.deg)
    origin = exp.getWcs().getPixelOrigin()
    print(f"{expId}  {ra.deg:.2f} {dec.deg:.2f}, {rotpa.deg:.2f}   {origin.getX():.2f}   {origin.getY():.2f}")
    
```

```python
exp.getBBox().getCenter()
```

```python

```
