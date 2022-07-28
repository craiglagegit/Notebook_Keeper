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
Modified by Craig Lage - 13-Jan-22\
Testing with comCam data from 2021-09-23.

```python
import sys
import matplotlib.pyplot as plt
import numpy as np
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
import lsst.daf.butler as dafButler #Gen3 butler
```

```python
def myCallback(im, ccd, imageSource):
    """Assemble the CCD image.  Just bias subtraction and gain correction"""
    oim = camGeomUtils.rawCallback(im, ccd, imageSource,
                                   subtractBias=True, correctGain=False)
    return oim
```

```python
# Instantiate the Gen3 butler
dataPath = "/repo/main"
instrument = "LATISS"
butler = dafButler.Butler(dataPath, 
                          collections=["LATISS/raw/all", "LATISS/calib/unbounded","u/cslage/calib/latiss/calib.20210217"],
                          instrument=instrument)
```

```python
day_obs = 20210218
seq_num = 699

dataId = {"instrument": instrument, "exposure.day_obs": day_obs, "exposure.seq_num": seq_num}
```

```python
# camera has the info necessary to assemble the 9 CCDs
camera = butler.get('camera', instrument=instrument)
# Print the metadata just as a check.
metadata = butler.get('raw.visitInfo', {**dataId, "detector.id": 0})
print(metadata)
```

```python
# Print out mean and sigma for each detector
# This is before bias subtraction and gain adjustment.
for det in range(1):
    exp = butler.get('raw', {**dataId, "detector.id": det})
    arr = arr = exp.image.array
    print(det, arr.min(), arr.max(), arr.mean(), arr.std())
```

```python jupyter={"outputs_hidden": true} tags=[]
# Instantiate the firefly display.
# This should open a new tab with Firefly
disp = afwDisplay.Display(0, "firefly")
```

```python
disp = afwDisplay.Display(0, "matplotlib")
```

```python
disp.frame
```

```python
# This is bias subtracted and gain adjusted (see myCallback above), so gives a decent image.
# Note that it is only ~ 200-300 counts after bias subtraction
fig = plt.figure(figsize=(16,16))
disp = afwDisplay.Display(6, "matplotlib")
disp.scale('linear', 'zscale')

dataType='raw' # 'raw' will look at the raw image, 'bias' looks at the master bias
mos = camGeomUtils.showCamera(camera,
                              camGeomUtils.ButlerImage(butler, dataType, 
                                                       instrument=instrument,
                                                       day_obs=day_obs, seq_num=seq_num,
                                                       verbose=True, callback=myCallback,
                                                       background=np.nan),
                              binSize=1, display=disp, overlay=False,
                              title="%d %d %s" % (day_obs, seq_num, dataType))
plt.savefig("/project/cslage/AuxTel/crosstalk/Image_2021021800699.png")
```

```python

```
