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

# Notebook for plotting the spline knots

Initially written 25 Feb 2022 by Craig Lage

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
# Force low flux knots to zero
# butler = Butler("/repo/main", collections=["LSSTCam/raw/all", "u/cslage/tests/linearizer_spline_E_25feb22"])
# Alex check
butler = Butler("/repo/main", collections=["LSSTCam/raw/all", "u/cslage/linearizer_28jan22"])

```

```python
expId=3021120600576
camera = butler.get('camera', instrument='LSSTCam')
det = 9# 94
lin = butler.get('linearizer', detector=det, exposure=expId, instrument='LSSTCam')
plt.subplot(1,1,1)
plt.title("Spline knots - 13144M - Detector %d"%det)
offset = 0.0
for it, amp in enumerate(camera[0].getAmplifiers()):
    centers, values = np.split(lin.linearityCoeffs[amp.getName()], 2)
    plt.scatter(centers, values + it * offset, marker='+')
    print(amp.getName(), centers, values)
    #break
plt.savefig("/repo/main/u/cslage/bps_13144M/plots/Spline_Knots_ITL_Det9_11Mar22.pdf")



```

```python

```
