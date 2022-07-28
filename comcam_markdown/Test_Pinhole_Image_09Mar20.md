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

# Notebook for ingesting ComCam pinhole images.

Initially written 04 Mar 2020 by Craig Lage\
This ingests the images into my own repo, \
and does assembly.  ISR is just CCD assembly, and bias subtraction.\
I applied the gains from last year's PTC measurements.\
The assembly of the raft is still manual at this point.

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
from lsst.ip.isr.isrTask import IsrTask
```

```python
# Set up the ISR task
# For now, this is just applying the bias and the gains
# For some reason, the darks are not working
isrConfig = IsrTask.ConfigClass()
isrConfig.doLinearize = False
isrConfig.doBias = True
isrConfig.doApplyGains = True
isrConfig.doFlat = False
isrConfig.doDark = False
isrConfig.doFringe = False
isrConfig.doDefect = False
isrConfig.doAddDistortionModel = False
isrConfig.doWrite = False
isrTask = IsrTask(config=isrConfig)

REPO_DIR = '/project/cslage/ComCam/20200303/'
butler = Butler(REPO_DIR)
```

```python
# Gains from PTC data from last year are now in the yaml file
# Assemble all 9 CCDs
# Run all three pinhole images in two colors with and without gaps.
# The gaps are realistic estimates of the gaps between CCDs

aSize = 0.32 # Size of the pixel array
xGap = 0.0072 # Gap between imaging arrays
yGap = 0.0032 # Gap between imaging arrays
visit = 3020030300034
color = 'gray'
plt.figure(figsize=(16,16))
xs = [0.0,aSize+xGap,2*(aSize+xGap),0.0,aSize+xGap,2*(aSize+xGap),0.0,aSize+xGap,2*(aSize+xGap)]
ys = [0.0,0.0,0.0,aSize+yGap,aSize+yGap,aSize+yGap,2*(aSize+yGap),2*(aSize+yGap),2*(aSize+yGap)]
for detector in range(9):
    dataRef = butler.dataRef('raw', detector=detector, visit=visit)
    postIsr = isrTask.runDataRef(dataRef).exposure
    ax=plt.axes([xs[detector],ys[detector],aSize*(509.0/500.0),aSize],aspect=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.imshow(np.log10(postIsr.image.array[0:4072,0:4000]),vmin=2.5, vmax=4.0, cmap=color)
plt.savefig(REPO_DIR+"images/Image_Log_09Mar20_yaml_Gain.png")
```

```python

```
