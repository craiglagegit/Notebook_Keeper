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
Craig Lage - 21-May-21

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
```

```python
from lsst.daf.butler import Butler
```

```python
butler = Butler('/repo/main', collections="LATISS/raw/all")
```

```python tags=[]
# Gen3 butler
dayObs = '2022-02-23'
dayObs = int(dayObs.replace('-', ''))

exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
    exposureList.append([record.id, record])
exposureList.sort(key=lambda x: x[0])
for [id,record] in exposureList:
    print(record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name)

```

```python tags=[]
FWHM_list
```

```python
CD1_1 0.0552614153694827
CD1_2 0.089281442480291
CD2_1 -0.089281442480291
CD2_2 0.0552614153694827
```

```python tags=[]
expId = 2022050300927
mData = butler.get('raw.metadata', detector=0, exposure=expId)
for key in mData.keys():
    print(key, mData[key])

```

```python
expId = 2021021800691
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
img = plt.imshow(arr, norm=LogNorm(vmin=12000, vmax=20000), interpolation='Nearest', cmap='gray')
colorbar(img)
plt.tight_layout(h_pad=1)
#plt.savefig(REPO_DIR+"/plots/NGC4755_17Feb21.pdf")
```

```python
np.median(arr)
```

```python
dayObs = 20210402
expId = 2021040200034
butler = Butler('/repo/main', collections="LSSTComCam/raw/all")
mData = butler.get('raw.metadata', detector=4, exposure=expId)
```

```python
for key in mData.keys():
    print(key, mData[key])

```

```python
expId = 2022030400004
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
img = plt.imshow(arr, norm=LogNorm(vmin=10, vmax=100000), interpolation='Nearest', cmap='gray')
colorbar(img)
plt.tight_layout(h_pad=1)
#plt.savefig(REPO_DIR+"/plots/NGC4755_17Feb21.pdf")
```

```python
arr.std()
```

```python
butler = Butler('/repo/main', collections=['LATISS/raw/all','LATISS/calib','u/cslage/calib/latiss/calib.2021021'])
```

```python

```

```python
expId = 2021021700077
exp = butler.get('raw', detector=0, exposure=expId)
```

```python
expId = 2021021700090
bias = butler.get('bias', detector=0, exposure=expId)
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib","u/cslage/calib/13144/calib.20220103",\
                                           "u/cslage/tests/linearizer_dm33297_21jan22"])
```

```python
expId = 3021120600576
```

```python
bias = butler.get('bias', detector=55, exposure=expId)
```

```python
defect = butler.get('defects', detector=55, exposure=expId)
```

```python
ptc = butler.get('ptc', detector=55, exposure=expId)
```

```python
lin = butler.get('linearizer', detector=55, exposure=expId)
```

```python

```
