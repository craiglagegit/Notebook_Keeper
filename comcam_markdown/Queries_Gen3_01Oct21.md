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
import lsst.daf.butler as dafButler
```

```python
butler = dafButler.Butler('/repo/main', collections="LSSTComCam/raw/all")
```

```python
dayObs = 20211014
detector = 4
exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
    #print(record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name)
    exposureList.append([record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name])
exposureList.sort()
for item in exposureList:
    print(item)
```

```python
dayObs = 20211014
detector = 4
exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
    #print(record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name)
    exposureList.append([record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name])
exposureList.sort()
for item in exposureList:
    print(item)
```

```python
temp = []
for i in range(2021093000079, 2021093000119):
    temp.append(i)
print(temp)
```

```python

```
