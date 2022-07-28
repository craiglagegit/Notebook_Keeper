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
! eups list -s cp_pipe
! eups list -s obs_lsst
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
                                       subtractBias=False, correctGain=False)
    return oim
```

```python
REPO_DIR = '/project/shared/comCam/'
butler = Butler(REPO_DIR)
expId = 2021040200054
exp = butler.get("raw", detector=4, expId=expId)
```

```python
info = exp.getInfo()
```

```python jupyter={"outputs_hidden": true}
dir(exp)
```

```python jupyter={"outputs_hidden": true}
mData = exp.getMetadata()
for key in mData.keys():
    print(key, mData[key])
```

```python
camera = butler.get('camera')
fig = plt.figure()
disp = afwDisplay.Display(1, "matplotlib")
disp.scale('linear', 0, max=1000)

dataType='raw'
mos = camGeomUtils.showCamera(camera, \
                              camGeomUtils.ButlerImage(butler, dataType, visit=expId, \
                                                        verbose=True, callback = myCallback,\
                                                      background = np.nan),\
                              binSize=1, display=disp, overlay=False, \
                              title="%d %s" % (expId, dataType))
#fig.savefig('/project/cslage/ComCam/20200303/plots/Pinhole_w_2020_13_08Apr20.png')
```

```python

```
