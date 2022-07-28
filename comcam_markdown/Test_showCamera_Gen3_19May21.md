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

# Notebook for testing showCamera.
Initially written by Michael Reuter.\
Testing with comCam data from 2021-04-02.

```python
import logging
import sys

import matplotlib.pyplot as plt
import numpy as np

import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
import lsst.daf.butler as dafButler

%matplotlib inline 
%config InlineBackend.figure_format = 'retina'

afwDisplay.setDefaultBackend("matplotlib")
```

```python
def myCallback(im, ccd, imageSource):
    """Assemble the CCD image.  Just bias subtraction and gain correction"""
    oim = camGeomUtils.rawCallback(im, ccd, imageSource,
                                   subtractBias=False, correctGain=True)
    return oim
```

```python
dataPath = "/repo/main"
instrument = "LSSTComCam"
butler = dafButler.Butler(dataPath, 
                          collections=["LSSTComCam/raw/all", "LSSTComCam/calib/unbounded"],
                          instrument=instrument)

```

```python
day_obs = 20210401
seq_num = 27
raftName = "R22"
dataId = {"instrument": instrument, "detector.raft": raftName,
          "exposure.day_obs": day_obs, "exposure.seq_num": seq_num}
```

```python
camera = butler.get('camera', instrument=instrument)
```

```python
metadata = butler.get('raw.visitInfo', {**dataId, "detector.id": 0})
```

```python
print(metadata)
```

```python
for det in range(9):
    exp = butler.get('raw', {**dataId, "detector.id": det})
    arr = arr = exp.image.array
    print(det, arr.min(), arr.max(), arr.mean(), arr.std())
```

```python
fig = plt.figure(figsize=(12,12))
disp = afwDisplay.Display(1, "matplotlib")
#disp.scale('asinh', 'zscale')
#disp.scale('linear', 0, max=8000)

dataType='raw'
mos = camGeomUtils.showCamera(camera,
                              camGeomUtils.ButlerImage(butler, dataType, 
                                                       instrument=instrument, raft=raftName,
                                                       day_obs=day_obs, seq_num=seq_num,
                                                       verbose=True, callback=myCallback,
                                                       background=np.nan),
                              binSize=16, display=disp, overlay=False,
                              title="%d %d %s" % (day_obs, seq_num, dataType))
```

```python

```
