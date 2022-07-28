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

# Notebook for querying ComCam data.

Initially written 13 Aug 2020 by Craig Lage

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
RAFT = 'R22'
SENSOR = 'S11'
```

```python
# First check the OCS dir
DATA_DIR = '/project/shared/comCam/'
filedir = DATA_DIR+'_parent/raw/'
files = glob.glob(filedir+'*/*/2020081200025-%s-%s-det004.fits'%(RAFT,SENSOR))
files = np.sort(files)
numFiles = len(files)
print(numFiles)
```

```python
# First check the exposure times
DATA_DIR = '/project/shared/comCam/'
filedir = DATA_DIR+'_parent/raw/'
files = glob.glob(filedir+'*/*/2020081300???-%s-%s-det004.fits'%(RAFT,SENSOR))
files = np.sort(files)
numFiles = len(files)
print(numFiles)

for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
        
    phdr=hdulist[0].header
    imhdr=hdulist[1].header
    filename = file.split('/')[8][8:13]#phdr['FILENAME']
    exptime = phdr['EXPTIME']
    imgtype = phdr['IMGTYPE'] 
    avg = 0.0#imhdr['AVERAGE']
    print("%s\t%s\t%f\t%f"%(filename, imgtype, exptime, avg))

```

```python
# OCS Primary headers
file = files[0]
hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
phdr_ocs = hdulist[0].header
dhdr_ocs = hdulist[1].header
dat_ocs = hdulist[1].data
for key in phdr_ocs.keys():
    print(key, phdr_ocs[key])

```

```python
# Now check the CCS dir
DATA_DIR = '/project/shared/comCam-CCS/'
filedir = DATA_DIR+'_parent/raw/'
files = glob.glob(filedir+'*/*/2020081200025-%s-%s-det004.fits'%(RAFT,SENSOR))
files = np.sort(files)
numFiles = len(files)
print(numFiles)
```

```python
# CCS Primary headers
file = files[0]
hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
phdr_ccs=hdulist[0].header
dhdr_ccs = hdulist[1].header
dat_ccs = hdulist[1].data
for key in phdr_ccs.keys():
    print(key, phdr_ccs[key])

```

```python
# CCS All headers
file = files[0]
hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)

for i in range(18):
    hdr=hdulist[i].header
    for key in hdr.keys():
        print(i,key, hdr[key])

```

```python
# Compare the primary headers
merged_keys = list(phdr_ocs.keys()) + list(set(phdr_ccs.keys()) - set(phdr_ocs.keys()))
for key in merged_keys:
    if (key in phdr_ocs.keys()) and  (key in phdr_ccs.keys()):
        print(key, 'OCS:', phdr_ocs[key],'CCS:', phdr_ccs[key])
    elif (key in phdr_ocs.keys()):
        print(key, 'OCS:', phdr_ocs[key])
    elif (key in phdr_ccs.keys()):
        print(key, 'CCS:', phdr_ccs[key])
```

```python
# Compare a pixel data value.  They should be the same and they are
print('OCS:',dat_ocs[500,500],'CCS:',dat_ccs[500,500])
```

```python
# Compare a data header
merged_keys = list(dhdr_ocs.keys()) + list(set(dhdr_ccs.keys()) - set(dhdr_ocs.keys()))
for key in merged_keys:
    if (key in dhdr_ocs.keys()) and  (key in dhdr_ccs.keys()):
        print(key, 'OCS:', dhdr_ocs[key],'CCS:', dhdr_ccs[key])
    elif (key in dhdr_ocs.keys()):
        print(key, 'OCS:', dhdr_ocs[key])
    elif (key in dhdr_ccs.keys()):
        print(key, 'CCS:', dhdr_ccs[key])
```

```python

```
