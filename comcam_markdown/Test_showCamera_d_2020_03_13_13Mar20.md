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

Initially written 09 Mar 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
! eups list -s obs_lsst
```

```python
from lsst.daf.persistence import Butler
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay

REPO_DIR = '/project/cslage/ComCam/20200303/'
butler = Butler(REPO_DIR)
visit = 3020030300034

def myCallback(im, ccd, imageSource):
    """Assemble the CCD image.  Just bias subtraction and gain correction"""
    oim = camGeomUtils.rawCallback(im, ccd, imageSource,
                                       subtractBias=True, correctGain=True)
    return oim

camera = butler.get('camera')
disp = afwDisplay.Display(0, "firefly")
disp.scale('asinh', 'zscale')

dataType='raw'
mos = camGeomUtils.showCamera(camera, \
                              camGeomUtils.ButlerImage(butler, dataType, visit=visit, \
                                                        verbose=True, callback = myCallback,\
                                                      background = 10000),\
                              binSize=16, display=disp, overlay=False, \
                              title="%d %s" % (visit, dataType))
```

```python

```
