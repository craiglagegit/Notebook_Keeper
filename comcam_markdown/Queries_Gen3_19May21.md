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
# Gen3 butler
from lsst.daf.butler import Butler
dayObs = 20210525
expId = 2021052500170
butler = Butler('/repo/main', collections="LATISS/raw/all")
mData = butler.get('raw.metadata', detector=0, exposure=expId)
```

```python jupyter={"outputs_hidden": true}
for key in mData.keys():
    print(key, mData[key])

```

```python
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
img = plt.imshow(arr, norm=LogNorm(vmin=14000, vmax=20000), interpolation='Nearest', cmap='gray')
colorbar(img)
plt.tight_layout(h_pad=1)
#plt.savefig(REPO_DIR+"/plots/NGC4755_17Feb21.pdf")
```

```python
dayObs = 20210602
expId = 2021060200050
butler = Butler('/repo/main', collections="LSSTComCam/raw/all")
mData = butler.get('raw.metadata', detector=4, exposure=expId)
```

```python
for key in mData.keys():
    print(key, mData[key])

```

```python
exp = butler.get('raw', detector=4, exposure=expId)
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

```python

```
