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
REPO_DIR = '/project/shared/comCam-CCS'
butler = Butler(REPO_DIR)
visit=3020112300030
dayObs = '2020-11-20'
```

```python
def myCallback(im, ccd, imageSource):
    """Assemble the CCD image.  Just bias subtraction and gain correction"""
    oim = camGeomUtils.rawCallback(im, ccd, imageSource,
                                       subtractBias=False, correctGain=False)
    return oim
```

```python
camera = butler.get('camera')
fig = plt.figure(figsize=(16,16))
disp = afwDisplay.Display(1, "matplotlib")
disp.scale('linear', min='zscale', max=None)
dataType='raw'
mos = camGeomUtils.showCamera(camera, \
                              camGeomUtils.ButlerImage(butler, dataType, visit=visit, \
                                                        verbose=True,  callback=myCallback, \
                                                      background = np.nan),\
                              title='%s, %s'%(visit,dataType),\
                              binSize=1, display=disp, overlay=False)
#fig.savefig(REPO_DIR+'/plots/Bias_%d_%s.png'%(visit,dayObs))
```

```python
means = []
stds = []
for det in range(9):
    exp = butler.get('raw', detector=det, dayObs=dayObs, visit=visit)
    mean = np.mean(exp.getMaskedImage().getArrays()[0])
    std = np.std(exp.getMaskedImage().getArrays()[0])
    print(det, mean, std)
    means.append(mean)
    stds.append(std)
mean = np.mean(means)
std = np.std(stds)
print(mean, std)
```

```python

```

```python
camGeomUtils.showCamera?
```

```python
afwDisplay.Display.scale?
```

```python
camGeomUtils.rawCallback?
```

```python

```
