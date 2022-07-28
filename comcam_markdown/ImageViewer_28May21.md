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

## Image viewer

In this notebook, we show several ways to look at images\
Craig Lage - 28-May-21

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
```

```python
# Instantiate the butler
REPO_DIR = '/project/shared/comCam'
butler = Butler(REPO_DIR)
```

```python
# Get the raw data
dayObs = '2021-05-28'
expId=2021052800003
for det in range(9):
    exp = butler.get('raw', detector=det, expId=expId)
    mean = np.mean(exp.getMaskedImage().getArrays()[0])
    std = np.std(exp.getMaskedImage().getArrays()[0])
    print(f"Detector {det}, Mean = {mean:.1f}, Std = {std:.1f}")
```

```python
# Look at the data with matplotlib
# The raw data doesn't look very good, because of the large pedestal of about 15,000 ADU
expId=2021052800003
det = 4
exp = butler.get('raw', detector=det, expId=expId)

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
plt.suptitle(f"Image {expId}",fontsize=18)
arr = exp.image.array
img = plt.imshow(arr, norm=LogNorm(vmin=14000, vmax=100000), interpolation='Nearest', cmap='gray')
colorbar(img)
plt.tight_layout(h_pad=1)
```

```python

```
