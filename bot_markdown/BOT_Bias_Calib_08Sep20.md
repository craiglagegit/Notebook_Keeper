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

# Notebook for viewing master bias images.

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
RAFT = 'R12'
SENSOR = 'S02'
data = {}
flux = "Hi"
#flux = "Low"
```

```python
# First, data from the 9 raft run
REPO_DIR = '/project/shared/BOT/rerun/cslage/PTC_6790D_NewAll_Amp/'
RUN = '6790D'
if flux == "Low":
    expId_1 = 3019101200332
    expId_2 = 3019101200333
elif flux == "Hi":
    expId_1 = 3019101200414
    expId_2 = 3019101200415
butler = Butler(REPO_DIR)
postISRCCD_1 = butler.get('postISRCCD',  raftName=RAFT,run=RUN, detectorName=SENSOR, expId=expId_1)
postISRCCD_2 = butler.get('postISRCCD',  raftName=RAFT,run=RUN, detectorName=SENSOR, expId=expId_2)
raw_1 = butler.get('raw',  raftName=RAFT,run=RUN, detectorName=SENSOR, expId=expId_1)
raw_2 = butler.get('raw',  raftName=RAFT,run=RUN, detectorName=SENSOR, expId=expId_2)
data[RUN] = [expId_1, expId_2, raw_1, raw_2, postISRCCD_1, postISRCCD_2]
```

```python
# Next, data from the 13 raft run
REPO_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_12543/'
RUN = '12543'
if flux == "Low":
    expId_1 = 3020090200352
    expId_2 = 3020090200353
elif flux == "Hi":
    expId_1 = 3020090200370
    expId_2 = 3020090200371
butler = Butler(REPO_DIR)
postISRCCD_1 = butler.get('postISRCCD',  raftName=RAFT,run=RUN, detectorName=SENSOR, expId=expId_1)
postISRCCD_2 = butler.get('postISRCCD',  raftName=RAFT,run=RUN, detectorName=SENSOR, expId=expId_2)
raw_1 = butler.get('raw',  raftName=RAFT,run=RUN, detectorName=SENSOR, expId=expId_1)
raw_2 = butler.get('raw',  raftName=RAFT,run=RUN, detectorName=SENSOR, expId=expId_2)
data[RUN] = [expId_1, expId_2, raw_1, raw_2, postISRCCD_1, postISRCCD_2]
```

```python
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
arr = e2v_bias.image.array
bias_mean = np.mean(arr)
bias_std = np.std(arr)
img1 = ax1.imshow(arr, vmin=0, vmax=10)
colorbar(img1)
ax1.set_title("E2V master bias %s%s; Mean = %.3f, Std = %.3f"%(E2V_RAFT,E2V_SENSOR,bias_mean, bias_std))
arr = itl_bias.image.array
bias_mean = np.mean(arr)
bias_std = np.std(arr)
img2 = ax2.imshow(arr, vmin=0, vmax=40)
ax2.set_title("ITL master bias %s%s; Mean = %.3f, Std = %.3f"%(ITL_RAFT,ITL_SENSOR,bias_mean, bias_std))
colorbar(img2)
plt.tight_layout(h_pad=1)
plt.savefig(E2V_REPO_DIR+"plots/Master_Biases_08Sep20.pdf")
plt.savefig(ITL_REPO_DIR+"plots/Master_Biases_08Sep20.pdf")
```

```python

```
