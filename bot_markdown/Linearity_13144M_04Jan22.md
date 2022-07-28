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
    return vendor, detectorNum
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib",\
                                                    "u/cslage/calib/13144/calib.20220103"])
camera = butler.get('camera', instrument='LSSTCam')
```

```python
expId = 3021120700200
fig = plt.figure(figsize=(12,6))
plt.suptitle("Linearizer spline knots - Run 13144M", fontsize=24)
plt.subplot(1,2,1)
plt.title("E2V", fontsize=18)
plt.xlabel("Flux(ADU)", fontsize=18)
plt.ylabel("Departure from linearity (ADU)", fontsize=18)
plt.xticks([0,50000,100000])
plt.xlim(0,100000)
plt.ylim(-1000,1000)
plt.subplot(1,2,2)
plt.title("ITL", fontsize=18)
plt.xlabel("Flux(ADU)", fontsize=18)
plt.ylabel("Departure from linearity (ADU)", fontsize=18)
plt.xticks([0,50000,100000])
plt.xlim(0,100000)
plt.ylim(-1000,1000)
plt.subplots_adjust(wspace=0.5)

for RAFT in ['R01',  'R02',  'R03', 'R10',  'R11',  'R12',  'R13', 'R14', 'R20',  'R21',  'R22',  'R23', 'R24', \
             'R30', 'R31', 'R32', 'R33', 'R34', 'R41', 'R42', 'R43']:
    for SENSOR in ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']:
        VENDOR, DETECTOR = detector(RAFT,SENSOR)
        try:
            lin = butler.get('linearizer', detector=DETECTOR, exposure=expId, instrument='LSSTCam')
        except:
            continue
        for amp in camera[0].getAmplifiers():
            ampName = amp.getName()
            centers, values = np.split(lin.linearityCoeffs[ampName], 2)
            if VENDOR == "E2V":
                plt.subplot(1,2,1)
                plt.scatter(centers, values, marker='.')
            elif VENDOR == "ITL":
                plt.subplot(1,2,2)
                plt.scatter(centers, values, marker='.')

plt.savefig("/repo/main/u/cslage/bps_13144M/plots/Spline_Knots_13144M_20Jan22.png")
```

```python

```
