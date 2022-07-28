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
Testing with comCam data from 2020-08-12.

```python
! eups list -s | grep lsst_distrib
! eups list -s cp_pipe
! eups list -s obs_lsst
```

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
```

```python
def myCallback(im, ccd, imageSource):
    """Assemble the CCD image.  Just bias subtraction and gain correction"""
    oim = camGeomUtils.rawCallback(im, ccd, imageSource,
                                       subtractBias=False, correctGain=False)
    return oim
```

```python
REPO_DIR = '/project/shared/comCam-CCS'
butler = Butler(REPO_DIR)
visit=2020081200025
```

```python
test=butler.get('raw', detector=4, visit=visit)
test=butler.get('raw', detector=4, dayObs='2020-08-12', seqNum=25)
```

```python
camera = butler.get('camera')
fig = plt.figure(figsize=(16,16))
disp = afwDisplay.Display(1, "matplotlib")
disp.scale('linear', min=20000, max=60000)
dataType='raw'
mos = camGeomUtils.showCamera(camera, \
                              camGeomUtils.ButlerImage(butler, dataType, visit=visit, \
                                                        verbose=True, callback = myCallback,\
                                                      background = np.nan),\
                              binSize=16, display=disp, overlay=False, \
                              title="%d %s" % (visit, dataType))
fig.savefig('/project/cslage/ComCam/20200812/plots/Flat_25_12Aug20.png')
```

```python

```

```python

```
