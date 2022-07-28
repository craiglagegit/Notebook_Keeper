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

## Queries - Gen3

In this notebook, we show several ways to query the Gen3 data\
Craig Lage - 21-May-21

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
```

```python
from lsst.daf.butler import Butler
```

```python
butler = Butler('/repo/main', collections="LATISS/raw/all")
```

```python jupyter={"outputs_hidden": true} tags=[]
# Gen3 butler
dayObs = '2022-02-16'
dayObs = int(dayObs.replace('-', ''))

exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
    exposureList.append([record.id, record])
exposureList.sort(key=lambda x: x[0])
for [id,record] in exposureList:
    print(record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name)

```

```python
exposureList[0]
```

```python jupyter={"outputs_hidden": true} tags=[]
for key in mData.keys():
    print(key, mData[key])
```

```python tags=[]
for exposure in exposureList:
    expId = exposure[0]
    mData = butler.get('raw.metadata', detector=0, exposure=expId)
    if mData['IMGTYPE'] not in ['SCIENCE', 'ENGTEST']:
        continue
    speed = abs((mData['AZSTART'] - mData['AZEND']) / mData['EXPTIME']) * 3600.0
    print(expId, mData['ELSTART'],mData['AZSTART'], mData['HASTART'], mData['DECSTART'], speed)
```

```python

```
