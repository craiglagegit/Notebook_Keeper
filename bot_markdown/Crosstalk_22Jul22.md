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

# Notebook for extracting crosstalk parameters from fits file.

Initially written 22 Jul 2022 by Craig Lage.

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.butler import Butler
```

```python
filename = '/project/cslage/BOT_LSSTCam/crosstalk/crosstalk_matrix_R32_S11_run13186.fits'
filename = '/project/cslage/BOT_LSSTCam/crosstalk/crosstalk_matrix_R10_S12_13175.fits'
```

```python
hdulist = pf.open(filename)
```

```python
hdulist[0].data.reshape((1,256))[0]
```

```python

```
