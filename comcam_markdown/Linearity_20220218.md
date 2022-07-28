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

# Notebook for investigating linearity corrections

Initially written 20 Dec 2021 by Craig Lage\
copying from Chris Waters.

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.butler import Butler
import lsst.afw.math as afwMath
```

```python
butler = Butler("/repo/main", collections=["LSSTComCam/raw/all","LSSTComCam/calib",\
                                                    "u/cslage/comcam/linearizerA_20220218"])
camera = butler.get('camera', instrument='LSSTComCam')
```

```python
expId = 2022021800078
fig = plt.figure(figsize=(12,6))
plt.suptitle("Linearizer spline knots - 20220218", fontsize=24)
plt.subplot(1,1,1)
plt.xlabel("Flux(ADU)", fontsize=18)
plt.ylabel("Departure from linearity (ADU)", fontsize=18)
plt.xticks([0,50000,100000])
plt.xlim(0,100000)
plt.ylim(-1000,1000)
plt.subplots_adjust(wspace=0.5)

for RAFT in ['R22']:
    for DETECTOR, SENSOR in enumerate(['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']):
        lin = butler.get('linearizer', detector=DETECTOR, exposure=expId, instrument='LSSTComCam')
        for amp in camera[0].getAmplifiers():
            ampName = amp.getName()
            centers, values = np.split(lin.linearityCoeffs[ampName], 2)
            plt.scatter(centers, values, marker='.')
plt.ylim(-500,500)
plt.savefig("/repo/main/u/cslage/comcam/ptc_20220218/plots/Spline_Knots_20220218.png")
```

```python

```
