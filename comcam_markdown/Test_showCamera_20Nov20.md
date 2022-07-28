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
#REPO_DIR = '/project/shared/comCam-CCS/rerun/cslage/ISR_ComCam_2020-12-28'
REPO_DIR = '/project/shared/comCam-CCS/rerun/cslage/PTC_2020-12-29'
butler = Butler(REPO_DIR)
visit=3020122900076
dayObs = '2020-12-29'
```

```python
means = []
stds = []
for det in range(9):
    exp = butler.get('postISRCCD', detector=det, dayObs=dayObs, visit=visit)
    mean = np.mean(exp.getMaskedImage().getArrays()[0])
    std = np.std(exp.getMaskedImage().getArrays()[0])
    print(det, mean, std)
    means.append(mean)
    stds.append(std)
mean = np.mean(means)
std = np.std(stds)
testType = exp.getMetadata()['TESTTYPE']
print(mean, std, testType)
```

```python
camera = butler.get('camera')
fig = plt.figure(figsize=(16,16))
disp = afwDisplay.Display(1, "matplotlib")
disp.scale('linear', min=0, max=50000)
dataType='postISRCCD'
mos = camGeomUtils.showCamera(camera, \
                              camGeomUtils.ButlerImage(butler, dataType, visit=visit, \
                                                        verbose=True,  \
                                                      background = np.nan),\
                              title='%s, %s, Mean = %.3f, Std = %.3f'%(visit,dataType,mean,std),\
                              binSize=32, display=disp, overlay=False)
fig.savefig(REPO_DIR+'/plots/%s_%d_%s.png'%(testType,visit,dayObs))
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

```
