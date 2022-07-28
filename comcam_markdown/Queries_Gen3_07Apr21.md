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

## Queries - Gen3

In this notebook, we show several ways to query the Gen3 data\
Craig Lage - 07-Apr-21

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import lsst.daf.butler as dafButler
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
```

```python jupyter={"outputs_hidden": true}
# Gen3 butler
REPO_DIR = '/repo/main'
butler = dafButler.Butler(REPO_DIR, collections="LSSTComCam/raw/all")

exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=20210402"):
    print(record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name)
    exposureList.append(record.id)
```

```python
dafButler.Butler?
```

```python
!ls /repo/main

```

```python
# Gen3 butler

butler = dafButler.Butler(REPO_DIR)
```

```python
list(butler.registry.queryCollections())
```

```python
dir(butler)
```

```python
# To look at the header keywords
expId = 2021040200005
exp = butler.get('raw', detector=0, exposure=expId, collections="u/cslage/ptc_20210402A", instrument='LSSTComCam')
mData = exp.getMetadata()
#for key in mData.keys():
#    print(key, mData[key])
```

```python

```

```python
# This also works, but the above is faster
for exposure in exposureList[7:8]:
    mData = butler.get('raw.metadata', detector=8, exposure=exposure)
    expTime = mData['EXPTIME']
    imgType = mData['IMGTYPE']
    obj = mData['OBJECT']
    print(exposure, expTime, imgType, obj)
    exp = butler.get('raw', detector=8, exposure=exposure)
    arr = exp.image.array
    print(arr.mean(), arr.std())
```

```python
expId = 2021031100134
exp = butler.get('raw', detector=0, exposure=expId)
```

```python
# Look at the data with matplotlib
# The raw data doesn't look very good, because of the large pedestal of about 15,000 ADU
from matplotlib.colors import LogNorm
# Now let's look at ithem
def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

plt.figure(figsize=(8,8))
plt.suptitle(f"Image",fontsize=18)
arr = exp.image.array
img = plt.imshow(arr, norm=LogNorm(vmin=14000, vmax=100000), interpolation='Nearest', cmap='gray')
colorbar(img)
plt.tight_layout(h_pad=1)
#plt.savefig(REPO_DIR+"/plots/NGC4755_17Feb21.pdf")
```
