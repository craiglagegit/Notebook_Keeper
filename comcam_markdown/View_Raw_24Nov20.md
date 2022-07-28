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

# Notebook for viewing postISRCCD images.

Initially written 29 May 2020 by Craig Lage\
Re-tested 08 Sep 2020 with latest code.

```python
! eups list -s | grep lsst_distrib
! eups list -s cp_pipe
```

```python
import sys, os, glob, subprocess
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
```

```python
DATA_DIR = '/project/shared/comCam-CCS'
RAFT = 'R22'
SENSOR ='S22'
DETECTOR = 4
expId1 = 3020112300030
expId2 = 3020112300044
butler = Butler(DATA_DIR)
```

```python
# Now let's get the raw files.
exp1 = butler.get('raw', raftName=RAFT, detectorName=SENSOR, expId=expId1)
exp2 = butler.get('raw', raftName=RAFT, detectorName=SENSOR, expId=expId2)
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

fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(16,8))
arr = exp1.image.array
mean = np.mean(arr)
std = np.std(arr)
img1 = ax1.imshow(arr, norm=LogNorm(vmin=16000, vmax=22000))
colorbar(img1)
ax1.set_title("Raw %s %s%s; Mean = %.3f, Std = %.3f"%(RAFT,SENSOR,expId1,mean, std))
arr = exp2.image.array
mean = np.mean(arr)
std = np.std(arr)
img2 = ax2.imshow(arr, norm=LogNorm(vmin=16000, vmax=22000))
colorbar(img1)
ax2.set_title("Raw %s %s%s; Mean = %.3f, Std = %.3f"%(RAFT,SENSOR,expId2,mean, std))

plt.tight_layout(h_pad=1)
#plt.savefig(E2V_REPO_DIR+"plots/Master_Biases_08Sep20.pdf")
#plt.savefig(ITL_REPO_DIR+"plots/Master_Biases_08Sep20.pdf")
```

```python
ccd1 = exp1.getDetector()
ampName = 'C04'
for ampObject in ccd1:

    amp = ampObject.getName()
    mask = exp1.getMask()
    maskimg1 = mask[ampObject.getBBox()].getArray()
    data = exp1[ampObject.getBBox()].image.array
    if amp == ampName:
        break
print(data.max(), data.min())
print(maskimg1[100,0:20])
plt.subplot(1,2,1)
plt.imshow(maskimg1[0:100,0:100])
plt.subplot(1,2,2)
plt.imshow(data[0:100,0:100])
plt.show()

```

```python
ccd1 = exp1.getDetector()
for ampObject in ccd1:
    amp = ampObject.getName()
    mask = exp1.getMask()
    maskimg1 = mask[ampObject.getBBox()].getArray()
    print(amp, maskimg1[100,0:15])

```

```python

```

```python
DATA_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_New_12606'
RAFT = 'R43'
SENSOR ='S12'
DETECTOR = 185
expId1 = 3020100800155
expId2 = 3020100800156
butler = Butler(DATA_DIR)
exp3 = butler.get('postISRCCD', raftName=RAFT, detectorName=SENSOR, expId=expId1)
```

```python
ccd3 = exp3.getDetector()
for ampObject in ccd1:
    amp = ampObject.getName()
    mask = exp3.getMask()
    maskimg3 = mask[ampObject.getBBox()].getArray()
    print(amp, maskimg3[100,0:15])

```

```python
print(type(exp1.getImage().array))
```

```python
def1 = butler.get('defects', raftName=RAFT, detectorName=SENSOR, expId=expId1)
```

```python
print(dir(def1))
```

```python
print(def1.toSimpleTable())
print(def1.getMetadata())
print(def1.toFitsRegionTable()[17])
print(dir(def1.toFitsRegionTable()))
```

```python
arr = def1.toFitsRegionTable().getAllBits().array
```

```python
print(arr.shape)
```

```python

```

```python
DATA_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_New_12606'
RAFT = 'R43'
SENSOR ='S01'
DETECTOR = 181
expId1 = 3020100800155
expId2 = 3020100800156
butler = Butler(DATA_DIR)
```

```python
def2 = butler.get('defects', raftName=RAFT, detectorName=SENSOR, expId=expId1)
print(def2.toSimpleTable())
```

```python
for i in range(1878):
    print(def1.toFitsRegionTable()[i])
```

```python

```
