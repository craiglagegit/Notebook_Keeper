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

# Notebook for viewing postISRCCD images.

Initially written 28 Sep 2020 by Craig Lage

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
def detector(raft, sensor):
    # Subroutine to find vendor and detector number given raft and sensor                                                                                                                                                           
    startingCol = [1,0,0,0,1] # First raft column in each row                                                                                                                                                                       
    rows = [0,3,8,13,18] # Starting raft sequence number of each row                                                                                                                                                                
    if raft in ['R11','R12','R13','R14','R21','R22','R23','R24','R30',\
                'R31','R32','R33','R34']:
        vendor = 'E2V'
    else:
        vendor = 'ITL'
    raftRow = int(list(raft)[1])
    raftCol = int(list(raft)[2]) - startingCol[raftRow]
    sensorRow = int(list(sensor)[1])
    sensorCol = int(list(sensor)[2])
    detectorNum = (rows[raftRow] + raftCol) * 9
    detectorNum += 3 * sensorRow + sensorCol
    plotNum = 21 - 5 * raftRow + int(list(raft)[2])
    return vendor, detectorNum, plotNum

```

```python
REPO_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_12606'
butler = Butler(REPO_DIR)
```

```python
rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

expIds = [3020110200100,3020110200101,3020110200103,3020110200104,3020110200130,3020110200131,3020110200106,3020110200107,3020110200133,3020110200134,3020110200109,3020110200110,3020110200136,3020110200137,3020110200112,3020110200113,3020110200139,3020110200140,3020110200115,3020110200116,3020110200142,3020110200143,3020110200145,3020110200146,3020110200148,3020110200149,3020110200118,3020110200119,3020110200151,3020110200152,3020110200154,3020110200155,3020110200157,3020110200158,3020110200121,3020110200122,3020110200160,3020110200161,3020110200163,3020110200164,3020110200166,3020110200167,3020110200169,3020110200170,3020110200124,3020110200125,3020110200172,3020110200173,3020110200175,3020110200176,3020110200178,3020110200179,3020110200127,3020110200128,3020110200181,3020110200182,3020110200184,3020110200185,3020110200187,3020110200188,3020110200190,3020110200191,3020110200193,3020110200194,3020110200196,3020110200197,3020110200199,3020110200200,3020110200202,3020110200203,3020110200205,3020110200206]

expIds = [3020100800155,3020100800156,3020100800158,3020100800159,3020100800185,3020100800186,\
          3020100800161,3020100800162,3020100800188,3020100800189,3020100800164,3020100800165,\
          3020100800191,3020100800192,3020100800167,3020100800168,3020100800194,3020100800195,\
          3020100800170,3020100800171,3020100800197,3020100800198,3020100800173,3020100800174,\
          3020100800200,3020100800201,3020100800176,3020100800177,3020100800203,3020100800204,\
          3020100800179,3020100800180,3020100800206,3020100800207,3020100800182,3020100800183,\
          3020100800209,3020100800210,3020100800212,3020100800213,3020100800215,3020100800216,\
          3020100800218,3020100800219,3020100800221,3020100800222]
```

```python
for RAFT in rafts:
    for SENSOR in sensors:
        VENDOR, DETECTOR, plotNum = detector(RAFT,SENSOR)
        counter = 0
        for expId in expIds:
            try:
                postISRCCD = butler.get('postISRCCD',  raftName=RAFT,detectorName=SENSOR, expId=expId)
                counter += 1
            except:
                continue
        print("Detector %d had %d good postISRCCDs"%(DETECTOR, counter))

```

```python

```
