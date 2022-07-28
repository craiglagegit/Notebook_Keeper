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

# Notebook for Calculating ComCam CCD Values

Initially written 29 Sep 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
! eups list -s | grep cp_pipe
```

```python
import sys, os, glob, time
import numpy as np
import matplotlib.pyplot as plt
import collections
from lsst.daf.persistence import Butler
```

```python
DIR = '/project/shared/comCam/rerun/cslage/PTC_2020-08-27/'
RAFT = 'R22'
```

```python
# Read noise, from PTC
butler = Butler(DIR)
read_noise = []
for detector in range(9):
    ptcDataset = butler.get('photonTransferCurveDataset', raftName=RAFT, detector=detector)
    noise = ptcDataset.noise
    for amp in noise.keys():
        read_noise.append(noise[amp])
print(np.median(read_noise), min(read_noise), max(read_noise))
```

```python
plt.hist(read_noise)
```

```python
# Dark current
dark_current = []
for detector in range(9):
    dark = butler.get('dark', raftName=RAFT, detector=detector, expId=2020082700015)
    ccd = dark.getDetector()
    for amp in ccd:
        img1 = dark.image
        arr1 = img1.Factory(img1, amp.getBBox()).array
        dark_current.append(np.median(arr1))
        #print(detector, amp.getName(), np.median(arr1))
print(np.median(dark_current), min(dark_current), max(dark_current))
```

```python
plt.hist(dark_current)
```

```python
# Bad pixels
bad_pixels = []
for detector in range(9):
    postISR = butler.get('postISRCCD', raftName=RAFT, detector=detector, expId=2020082700067)
    ccd =postISR.getDetector()
    bad = 0
    for amp in ccd:
        img1 = postISR.getMaskedImage()
        mask = img1.Factory(img1, amp.getBBox()).getMask()
        bad += (collections.Counter(mask.array.flatten())[5]) # "Bad pixels"
    #print(detector, bad)
    bad /= (2000*509*16)
    bad *= 100 # To put in %
    bad_pixels.append(bad)

print(np.median(bad_pixels), min(bad_pixels), max(bad_pixels))
```

```python
# Now Serial CTI
# These are for the ITL sensor, identifying the end of the imaging region and the overscan region.
xstart = 505
xstop = 542
ov_start = 512
# Run most of the rows, but stay away from the edges.
ystart = 200
ystop = 1800
xaxis = np.linspace(xstart,xstop-1,xstop-xstart)
ctis = []
for detector in range(9):
    raw = butler.get('raw', raftName=RAFT, detector=detector, expId=2020082700082)
    ccd = raw.getDetector()
    for amp in ccd:
        img1 = raw.image
        data = img1.Factory(img1, amp.getRawBBox()).array
        data = np.flip(data, axis=1)
        flat_overscan = np.mean(data[:,xstop-8:xstop],axis = 1)
        cte_data = ((np.transpose(np.transpose(data) - flat_overscan))[ystart:ystop,:].mean(axis=0))[xstart:xstop]
        cte_std = ((np.transpose(np.transpose(data) - flat_overscan))[ystart:ystop,:].std(axis=0) / np.sqrt(float(ystop-ystart)))[xstart:xstop]
        cti = np.median((np.transpose(np.transpose(data) - flat_overscan))[ystart:ystop,ov_start]\
        / (np.transpose(np.transpose(data) - flat_overscan))[ystart:ystop,ov_start-1]) / ov_start
        ctis.append(cti)
print(np.median(ctis), max(ctis), min(ctis))
```

```python
plt.hist(ctis)
```

```python

```
