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
from astropy.time import Time
from lsst.daf.persistence import Butler
```

```python
RAFT = 'R22'
SENSOR = 'S11'
DATA_DIR = '/project/shared/comCam-CCS/'
```

```python jupyter={"outputs_hidden": true}
REPO_DIR = '/project/shared/comCam-CCS/rerun/cslage/PTCTemp_2020-12-29'
butler = Butler(REPO_DIR)
```

```python
expId = 3020122900047
butler.get('raw', raftName=RAFT,detector=3, expId=expId)
```

```python
expId = 3020122900047
butler.get('bias', raftName=RAFT,detectorName=SENSOR, expId=expId)
```

```python
expId = 3020122900047
butler.get('dark', raftName=RAFT,detectorName=SENSOR, expId=expId)
```

```python
expId = 3020122900047
butler.get('flat', raftName=RAFT,detectorName=SENSOR, expId=expId)
```

```python
expId = 3020122900047
butler.get('defects', raftName=RAFT,detectorName=SENSOR, expId=expId)
```

```python
dayObs='2021-01-11'
butler = Butler(DATA_DIR)
visits = butler.queryMetadata('raw', ['expId', 'EXPTIME', 'TESTTYPE', 'DATE'], raftName=RAFT,\
                              detectorName=SENSOR, dayObs=dayObs)
visits.sort(key = lambda x: x[3]) 

for (expId, exptime, testtype, date) in visits:
    print(date, expId, exptime, testtype)
```

```python
# CCS All headers
filedir = DATA_DIR+'_parent/raw/'
files = glob.glob(filedir+'*/*/202101110000?-%s-%s-det004.fits'%(RAFT,SENSOR))
files = np.sort(files)
numFiles = len(files)
print(numFiles)

expIds = []
CCS_unix = []
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    for i in range(1):
        hdr=hdulist[i].header
        for key in hdr.keys():
            if key == 'DATE-OBS':
                print(i,file[50:63],key, hdr[key])
                expIds.append(file[50:63])
                #CCS_unix.append(Time(hdr[key]).unix)
                CCS_unix.append(hdr[key])

```

```python
RAFT = 'R22'
SENSOR = 'S11'
DATA_DIR = '/project/shared/comCam/'
```

```python
dayObs='2021-01-11'
butler = Butler(DATA_DIR)
visits = butler.queryMetadata('raw', ['expId', 'EXPTIME', 'TESTTYPE', 'DATE'], raftName=RAFT,\
                              detectorName=SENSOR, dayObs=dayObs)
visits.sort(key = lambda x: x[3]) 

for (expId, exptime, testtype, date) in visits:
    print(date, expId, exptime, testtype)
```

```python
# OCS All headers
filedir = DATA_DIR+'_parent/raw/'
files = glob.glob(filedir+'*/*/202101110000?-%s-%s-det004.fits'%(RAFT,SENSOR))
files = np.sort(files)
numFiles = len(files)
print(numFiles)

OCS_unix = []
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    for i in range(1):
        hdr=hdulist[i].header
        for key in hdr.keys():
            if key == 'DATE-OBS':
                print(i,file[30:63],key, hdr[key])
                #OCS_unix.append(Time(hdr[key]).unix)
                OCS_unix.append(hdr[key])

```

```python
print('expId             OCS_DATE_OBS - CCS_DATE_OBS')
for i, ocs_t in enumerate(OCS_unix):
    print('%s \t \t %f'%(expIds[i+2], ocs_t - CCS_unix[i+2]))
```

```python
print('expId              OCS_DATE_OBS                    CCS_DATE_OBS')
for i, ocs_t in enumerate(OCS_unix):
    print('%s \t %s \t %s'%(expIds[i+2], ocs_t , CCS_unix[i+2]))
```

```python

```
