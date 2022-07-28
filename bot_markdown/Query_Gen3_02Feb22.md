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

## Query header - Gen3

In this notebook, we show several ways to query the fits headers\
Craig Lage - 16-Mar-21

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import lsst.daf.butler as dafButler
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
```

```python jupyter={"outputs_hidden": true} tags=[]
# Gen3 butler
REPO_DIR = '/repo/main'
butler = dafButler.Butler(REPO_DIR, collections=["LSSTCam/raw/all","LSSTCam/calib"])
 
exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.science_program='13144'"):
    exposureList.append([record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name])
exposureList = sorted(exposureList, key=lambda x: x[0])   
for item in exposureList:
    print(f"{item[0]} \t {item[1]} \t {item[2]} \t\t {item[3]} \t {item[4]}")
```

```python jupyter={"outputs_hidden": true} tags=[]
# Gen3 butler
REPO_DIR = '/repo/main'
butler = dafButler.Butler(REPO_DIR, collections=["LSSTCam/raw/all","LSSTCam/calib"])
 
exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.science_program='13117'"):
    exposureList.append([record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name])
exposureList = sorted(exposureList, key=lambda x: x[0])   
for item in exposureList:
    print(f"{item[0]} \t {item[1]} \t {item[2]} \t\t {item[3]} \t {item[4]}")
```

```python jupyter={"outputs_hidden": true} tags=[]
# To look at the header keywords
#expId = 3021020500605
#expId = 3021120600565
expId=3021120600612
exp = butler.get('raw',  exposure=expId, detector=55, instrument="LSSTCam")
mData = exp.getMetadata()
for key in mData.keys():
    print(key, mData[key])
```

```python jupyter={"outputs_hidden": true} tags=[]
hdulist = pf.open("/lsstdata/offline/instrument/LSSTCam-bot/storage/20211206/MC_C_20211206_000612/MC_C_20211206_000612_R13_S01.fits")
for i in range(len(hdulist)):
    hdr = hdulist[i].header
    for key in hdr.keys():
        print(i, key, hdr[key])
```

```python
# This is useful for finding what datasets the butler knows about
registry = butler.registry
for x in registry.queryDatasetTypes():
    if 'isr' in x.name:
        print(x)

```

```python

```
