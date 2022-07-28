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
# First check the OCS dir
DATA_DIR = '/project/shared/auxTel/'
filedir = DATA_DIR+'_parent/raw/'
files = glob.glob(filedir+'*/2020031600212-det000.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)
```

```python
# All headers
file = files[0]
hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)

for i in range(18):
    hdr=hdulist[i].header
    for key in hdr.keys():
        print(i,key, hdr[key])

```

```python

```
