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
#REPO_DIR = '/project/shared/comCam-CCS/rerun/cslage/ISR_ComCam_2020-12-28'
#REPO_DIR = '/project/shared/auxTel/rerun/cslage/PTC_Defect_2021-02-17'
#butler = Butler(REPO_DIR)
dayObs = '2021-02-18'
```

```python jupyter={"outputs_hidden": true}
#dayObs = '2021-02-17'
dayObs = '2021-02-18'
#expId=2021021700349
expId=2021021800693
det = 0
REPO_DIR = '/project/shared/auxTel/rerun/cslage/PTC_OvBias+Only_%s'%dayObs
butler = Butler(REPO_DIR)
exp2 = butler.get('postISRCCD', detector=det, dayObs=dayObs, expId=expId)
mean2 = np.mean(exp2.getMaskedImage().getArrays()[0])
std2 = np.std(exp2.getMaskedImage().getArrays()[0])

#print(mean1, mean2)
```

```python
REPO_DIR = '/project/shared/auxTel/rerun/cslage/PTC_Defect_%s'%dayObs
butler = Butler(REPO_DIR)

exp1 = butler.get('postISRCCD', detector=det, dayObs=dayObs, expId=expId)
mean1 = np.mean(exp1.getMaskedImage().getArrays()[0])
std1 = np.std(exp1.getMaskedImage().getArrays()[0])
```

```python
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

plt.figure(figsize=(16,8))
#plt.suptitle("Jewel Box Cluster (NGC4755) - 17-Feb-2021 - Rubin AuxTel",fontsize=18)
plt.suptitle("Crosstalk Check - 18-Feb-2021 - Rubin AuxTel",fontsize=18)
arr1 = exp1.image.array
plt.subplot(1,2,2)
plt.title("Overscan, Bias, Defects, Gains, Saturation, and Crosstalk")
img1 = plt.imshow(arr1, norm=LogNorm(vmin=5, vmax=150000), cmap='gray')
colorbar(img1)
plt.subplot(1,2,1)
plt.title("Overscan and Bias only")
arr2 = exp2.image.array
img2 = plt.imshow(arr2, norm=LogNorm(vmin=5, vmax=150000), cmap='gray')
colorbar(img2)
plt.tight_layout(h_pad=1)
plt.savefig(REPO_DIR+"/plots/Crosstalk_%d_17Feb21.pdf"%expId)

```

```python

```
