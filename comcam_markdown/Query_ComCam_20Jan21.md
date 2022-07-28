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

# Notebook for querying BOT data.

Initially written 27 May 2020 by Craig Lage\
Allows inspecting the image type and exposure time of the \
BOT images used for characterizing BF.

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
DATA_DIR = '/project/shared/comCam/'
```

```python
DATA_DIR = '/project/shared/comCam/'
dayObs='2021-05-25'
detector=4
butler = Butler(DATA_DIR)
visits = butler.queryMetadata('raw', ['expId', 'EXPTIME', 'TESTTYPE', 'FILTER'], detector=detector, dayObs=dayObs)
visits.sort(key = lambda x: x[0]) 

for (expId, exptime, testtype, filter) in visits:
    print(expId, exptime, testtype, filter)
```

```python
expId = 2021040200028
RERUN_DIR = '/project/shared/comCam/rerun/cslage/PTC_2021-04-02/'
butler = Butler(RERUN_DIR)
test = butler.get('raw', detector=4,  expId=expId)
```

```python
test = butler.get('bias', detector=4,  expId=expId)
```

```python
test = butler.get('dark', detector=4,  expId=expId)
```

```python
test = butler.get('defects', detector=0,  expId=expId)
```

```python
test = butler.get('defects', detector=2,  expId=expId)
```

```python
test = butler.get('defects', dataId={'detector': 2, 'expId': 2021040200025, 'dayObs': '2021-04-02', 'raftName': 'R22', 'detectorName': 'S02'})
```

```python
test = butler.get('defects', dataId={'detector': 2, 'expId': 2021040200025, 'dayObs':'2021-04-02', 'raftName':'R22'})
```

```python
# RUN = '12597'
butler = Butler(DATA_DIR)
visits = butler.queryMetadata('raw', ['run', 'expId', 'EXPTIME', 'TESTTYPE'], raftName=RAFT,\
                              detectorName=SENSOR, dayObs='2020-10-06')
for (run, expId, exptime, testtype) in visits:
    print(run, expId, exptime, testtype)
```

```python
butler = Butler(DATA_DIR)
visits = butler.queryMetadata('raw', ['expId', 'EXPTIME', 'TESTTYPE'], raftName=RAFT,\
                              detectorName=SENSOR, dayObs='2020-08-28')
for (expId, exptime, testtype) in visits:
    print(expId, exptime, testtype)
```

```python
RUN = '12539'
butler = Butler(DATA_DIR)
dets = butler.queryMetadata('raw', ['raftName','detector'],run=RUN)
print("Total CCDs = %d"%len(dets))
imDets = []
for [raftName, det] in dets:
    if raftName in ['R00', 'R04', 'R40', 'R44']:
        continue
    else:
        imDets.append(det)
print("Imaging CCDs = %d"%len(imDets))
print(np.sort(imDets))
```

```python
# First check the flats
DATA_DIR = '/project/shared/comCam-CCS/'
filedir = DATA_DIR+'_parent/raw/'
files = glob.glob(filedir+'*/*/3020122800???-%s-%s-det004.fits'%(RAFT,SENSOR))
files = np.sort(files)
numFiles = len(files)
print(numFiles)
print(files[0])

for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    phdr=hdulist[0].header
    filenumber = file.split('/')[-1][0:13]
    seq = int(file.split('/')[-1][8:13])
    try:
        exptime = phdr['EXPTIME']
    except:
        exptime = 'None'
    try:
        filter = phdr['FILTER']
    except:
        filter = 'None'
    imgtype = phdr['IMGTYPE'] 
    print(filenumber, seq, imgtype, exptime, filter)
```

```python
for key in phdr.keys():
    print(key, phdr[key])
```

```python

```
