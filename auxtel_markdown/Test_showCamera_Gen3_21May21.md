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
import sys
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
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
                                   subtractBias=True, correctGain=False)

```

```python
dataPath = "/repo/main"
instrument = "LATISS"
butler = dafButler.Butler(dataPath, 
                          collections=["LATISS/raw/all",  "LATISS/calib/unbounded"],
                          instrument=instrument)
```

```python
day_obs = 20210521
seq_num = 22
dataId = {"instrument": instrument,
          "exposure.day_obs": day_obs, "exposure.seq_num": seq_num}
```

```python
camera = butler.get('camera', instrument=instrument)
```

```python
butler.get('raw.visitInfo', {**dataId, "detector.id": 0})
```

```python
fig = plt.figure(figsize=(12,12))
disp = afwDisplay.Display(1, "matplotlib")
#disp.scale('asinh', 'zscale')
#disp.scale('linear', 0, max=8000)

dataType='raw'
mos = camGeomUtils.showCamera(camera,
                              camGeomUtils.ButlerImage(butler, dataType, dataId=dataId,
                                                       verbose=True, callback=myCallback,
                                                       background=np.nan),
                              binSize=16, display=disp, overlay=False,
                              title="%d %d %s" % (day_obs, seq_num, dataType))
```

```python

```
