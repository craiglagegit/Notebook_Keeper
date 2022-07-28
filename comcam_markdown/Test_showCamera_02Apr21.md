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

```python jupyter={"outputs_hidden": true}
REPO_DIR = '/project/shared/comCam'
butler = Butler(REPO_DIR)
expId=2020102100043

```

```python
expId=2021040100055
dayObs = '2021-04-02'
```

```python
exp = butler.get('raw', detector=4,  expId=expId)
```

```python
def myCallback(im, ccd, imageSource):
    """Assemble the CCD image.  Just bias subtraction and gain correction"""
    oim = camGeomUtils.rawCallback(im, ccd, imageSource,
                                       subtractBias=True, correctGain=False)
    return oim
```

```python
camGeomUtils.rawCallback?
```

```python
camera = butler.get('camera')
fig = plt.figure(figsize=(16,16))
disp = afwDisplay.Display(1, "matplotlib")
#disp.scale('linear', min='zscale', max=None)
disp.scale('asinh', min='zscale', max=None)
dataType='raw'
mos = camGeomUtils.showCamera(camera, \
                              camGeomUtils.ButlerImage(butler, dataType, expId=expId, \
                                                        verbose=True,  callback=myCallback, \
                                                      background = np.nan),\
                              title='%s, %s'%(visit,dataType),\
                              binSize=4, display=disp, overlay=False)
#fig.savefig(REPO_DIR+'/plots/Bias_%d_%s.png'%(visit,dayObs))
```

```python
means = []
stds = []
for det in range(9):
    exp = butler.get('raw', detector=det, expId=expId)
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
