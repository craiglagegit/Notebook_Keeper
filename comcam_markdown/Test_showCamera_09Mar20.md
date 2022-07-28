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
import sys, os, glob
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf

from lsst.daf.persistence import Butler
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
```

```python
REPO_DIR = '/project/cslage/ComCam/20200303/'
#REPO_DIR = '/lsstdata/offline/teststand/comcam/'
#REPO_DIR = '/project/shared/comCam'
butler = Butler(REPO_DIR)
#visit = 2019111300120
visit = 3020030300034
```

```python
camera = butler.get('camera')
imageSource = camGeomUtils.ButlerImage(butler, 'raw', visit=visit)
test = camGeomUtils.showCamera(camera, imageSource=imageSource)
```

```python
print(test.array.shape)
print(test.array.max(), test.array.min())
```

```python
 
```

```python
plt.imshow(test.array,vmin=10, vmax=10000, cmap='gray')
```

```python
print(dir(camera))
```

```python
print(camera.getNameMap())
```

```python
print(dir(camera.getNameMap()['R22_S11'][0]))
```

```python
print(camera.getNameMap()['R22_S11'][0].getGain())
```

```python
test = camGeomUtils.ButlerImage(butler, 'raw', visit=visit, background = 10.0, \
                                                        verbose=True, callback = myCallback)
```

```python
print(camera.getNameMap()['R22_S11'][0].getRawBBox())
```

```python
ccd = camera.getNameMap()['R22_S11']
```

```python
test2 = test.getCcdImage(ccd)
```

```python
print(test2)
```

```python
def myCallback(im, ccd, imageSource):
    """Assemble the CCD image"""
    oim = camGeomUtils.rawCallback(im, ccd, imageSource,
                                       subtractBias=True, correctGain=True)
    return oim
```

```python
afwDisplay.Display?
```

```python
disp = afwDisplay.Display(1, "matplotlib")
disp.scale('asinh', 'zscale')

dataType='raw'
mos = camGeomUtils.showCamera(camera, \
                              camGeomUtils.ButlerImage(butler, dataType, visit=visit, \
                                                        verbose=True, callback = myCallback,\
                                                      background = np.nan),\
                              binSize=16, display=disp, overlay=False, \
                              title="%d %s" % (visit, dataType))

print(mos.array)
```

```python
raw = butler.get('raw', detector=4, visit=visit)
```

```python
print(raw.image.array.max())
```

```python
# Finally have it working.  Need to keep working on getting things in the right place and 
# getting the gaps right.

# It looks like the file to edit is in policy/comCam.yaml, but it seems this file gets built.
# Where does the offset get stored?
```

```python
#Check if numbering is correct
for detector in range(9):
    raw = butler.get('raw', detector=detector, visit=visit)
    ccd = raw.getDetector()
    ccdName = ccd.getName()
    print(detector, ccdName, raw.getMetadata()['LSST_NUM'])
```

```python
# From R22.yaml
"""
    S00 : ITL-3800C-229
    S01 : ITL-3800C-251
    S02 : ITL-3800C-215
    S10 : ITL-3800C-326
    S11 : ITL-3800C-283
    S12 : ITL-3800C-243
    S20 : ITL-3800C-319
    S21 : ITL-3800C-209
    S22 : ITL-3800C-206
"""
# Looks correct


```
