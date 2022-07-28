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

# Notebook for querying defect files.

Initially written 18 Jan 2021 by Craig Lage\

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
RAFT = 'R20'
SENSOR = 'S00'
DATA_DIR = '/project/shared/BOT/'
```

```python
# First, let's look at some old data
filedir = DATA_DIR+'rerun/cslage/PTC_LSSTCAM_12673/calibrations/LSSTCam/defects/r20_s00/'
files = glob.glob(filedir+'*.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)
```

```python
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    for i, hdu in enumerate(hdulist):
        hdr=hdu.header
        for key in hdr.keys():
            print(i, key, hdr[key])

```

```python
tableData = hdulist[1].data 
print(len(tableData))
```

```python
print(tableData[:10])  # show the first ten rows

```

```python
# Now let's look at the same data after ingest
filedir = DATA_DIR+'rerun/cslage/PTC_LSSTCAM_12673/CALIB/defects/*/*-R20-S00*'
files = glob.glob(filedir+'*.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)
```

```python
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    for i, hdu in enumerate(hdulist):
        hdr=hdu.header
        for key in hdr.keys():
            print(i, key, hdr[key])

```

```python
tableData = hdulist[1].data 
print(len(tableData))
```

```python
print(tableData[:10])  # show the first ten rows

```

```python

```

```python
# Now let's look at the new data
DATA_DIR = '/project/shared/comCam-CCS/'
filedir = DATA_DIR+'rerun/cslage/PTCTemp_2020-12-29/defects/*/*-R22-S11*'
files = glob.glob(filedir+'*.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)
```

```python
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    for i, hdu in enumerate(hdulist):
        hdr=hdu.header
        for key in hdr.keys():
            print(i, key, hdr[key])

```

```python
tableData = hdulist[1].data 
print(len(tableData))
```

```python
print(tableData[:10])  # show the first two rows

```

```python
# Now let's look at the new data after ingest
DATA_DIR = '/project/shared/comCam-CCS/'
filedir = DATA_DIR+'rerun/cslage/PTCTemp_2020-12-29/CALIB/defects/*/*-R22_S11*'
files = glob.glob(filedir+'*.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)
```

```python
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    for i, hdu in enumerate(hdulist):
        hdr=hdu.header
        for key in hdr.keys():
            print(i, key, hdr[key])

```

```python
tableData = hdulist[1].data 
print(len(tableData))
```

```python jupyter={"source_hidden": true}
print(tableData[:10])  # show the first two rows

```

```python

```
