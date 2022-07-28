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
butler = Butler("/repo/main", collections=[ "u/cslage/bps_13144G"])
expId=3021120600576
camera = butler.get('camera', instrument='LSSTCam')
```

```python
# Not succeeding in using the butler, so I'm just going directly to the fits files.

names = ["E2V", "ITL"]
for i, det in enumerate([55, 74]):
    plt.subplot(1,2,i+1)
    plt.title(names[i]+f" - {det}")
    if names[i] == "E2V":
        filename = "/repo/main/u/cslage/bps_13144K/20211230T012140Z/linearizer/linearizer_LSSTCam_R13_S01_u_cslage_bps_13144K_20211230T012140Z.fits"
    elif names[i] == "ITL":
        filename = "/repo/main/u/cslage/bps_13144K/20211230T012140Z/linearizer/linearizer_LSSTCam_R20_S02_u_cslage_bps_13144K_20211230T012140Z.fits"
    hdu = pf.open(filename)
    data = hdu[1].data

    for amp in range(16):

        centers, values = np.split(data["COEFFS"][amp], 2)
        plt.scatter(centers, values, marker='+')

    plt.subplots_adjust(wspace=0.5)
```

```python
# Not succeeding in using the butler, so I'm just going directly to the fits files.

names = ["E2V", "ITL"]
for i, det in enumerate([55, 74]):
    plt.subplot(1,2,i+1)
    plt.title(names[i]+f" - {det}")
    if names[i] == "E2V":
        filename = "/repo/main/u/cslage/bps_13144K/20211230T012140Z/linearizer/linearizer_LSSTCam_R13_S01_u_cslage_bps_13144K_20211230T012140Z.fits"
    elif names[i] == "ITL":
        filename = "/repo/main/u/cslage/bps_13144K/20211230T012140Z/linearizer/linearizer_LSSTCam_R20_S02_u_cslage_bps_13144K_20211230T012140Z.fits"
    hdu = pf.open(filename)
    data = hdu[1].data

    for amp in range(1
                    ):

        centers, values = np.split(data["COEFFS"][amp], 2)
        plt.scatter(centers, values, marker='+')

    plt.subplots_adjust(wspace=0.5)
```

```python

```
