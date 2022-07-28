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
Initially written 09 Mar 2020 by Craig Lage.\
Testing w_2020_13 code.

```python
! eups list -s | grep lsst_distrib
! eups list -s | grep obs_lsst
```

```python
import numpy as np
import matplotlib.pyplot as plt
from lsst.daf.persistence import Butler
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
```

```python
def myCallback(im, ccd, imageSource):
    """Assemble the CCD image.  Just bias subtraction and gain correction"""
    oim = camGeomUtils.rawCallback(im, ccd, imageSource,
                                       subtractBias=True, correctGain=True)
    return oim
```

```python
REPO_DIR = '/project/cslage/ComCam/20200303/'
#REPO_DIR = '/lsstdata/offline/teststand/comcam/CCS/gen2repo'
butler = Butler(REPO_DIR)
#visit = 3020030300064
visit = 3020030300034
dataId = dict(visit=visit, detector=94)
camera = butler.get('camera')

fig = plt.figure(figsize=(16,16))
disp = afwDisplay.Display(1, "matplotlib")
disp.scale('linear', 0, max=8000)

dataType='raw'
mos = camGeomUtils.showCamera(camera, \
                              camGeomUtils.ButlerImage(butler, dataType, visit=dataId["visit"], \
                                                        verbose=False, callback = myCallback,\
#                                                        verbose=False, callback = camGeomUtils.rawCallback,\
                                                      background = np.nan),\
                              binSize=16, detectorNameList=[90,91,92,93,94,95,96,97,98], display=disp, overlay=False, \
                              title="%d %s" % (visit, dataType))
#fig.savefig(REPO_DIR+'plots/Pinhole_w_2020_13_04Apr20.png')
```

```python
#REPO_DIR = '/project/cslage/ComCam/20200303/'
REPO_DIR = '/lsstdata/offline/teststand/comcam/CCS/gen2repo'
butler = Butler(REPO_DIR)
visit = 3020030300034
camera = butler.get('camera')

fig = plt.figure(figsize=(8,8))
disp = afwDisplay.Display(1, "matplotlib")
disp.scale('linear', 0, max=8000)

dataType='raw'
mos = camGeomUtils.showCamera(camera, \
                              camGeomUtils.ButlerImage(butler, dataType, visit=visit, \
                                                        verbose=True, callback = myCallback,\
                                                      background = np.nan),\
                              binSize=16, display=disp, overlay=False, \
                              title="%d %s" % (visit, dataType))
#fig.savefig(REPO_DIR+'plots/Pinhole_w_2020_13_04Apr20.png')
```

```python

```
