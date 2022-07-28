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

# Notebook for re-ingesting ComCam data.

Initially written 15 Nov 2019 by Craig Lage\
This ingests the images into my own repo, \
creates the master bias images, and ingests them.\
02-Dec-19 - I've added the master flats and master darks.
19-Dec-19  Adding code to add missing header keywords\
I'm setting up for BF correction, so this requires the following:\
obs_base: tickets/DM-18683\
obs_lsst: tickets/DM-18683\
cp_pipe: tickets/DM-18683\
ip_isr: tickets/DM-22659

```python
! eups list -s | grep lsst_distrib
```

```python
import eups
import sys, os, glob
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf

from lsst.daf.persistence import Butler
from lsst.ip.isr.isrTask import IsrTask
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python
DATA_DIR = '/project/shared/comCam/'
REPO_DIR = '/project/cslage/ComCam/20191218/'
OUTPUT_DIR = '/project/cslage/ComCam/20191218/'
DETECTOR = 4
raftName = 'R22'
```

```python
# First check the exposure times
filedir = DATA_DIR+'raw/20191218/'
files = glob.glob(filedir+'CC_C_20191218_00????/CC_C_20191218_00????_R22_S11.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)

for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
        
    phdr=hdulist[0].header
    filename = file.split('/')[6][14:]#phdr['FILENAME']
    exptime = phdr['EXPTIME']
    imgtype = phdr['IMGTYPE'] 
    print("%s\t%s\t%f"%(filename, imgtype, exptime))


```

```python

```
