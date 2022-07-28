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
REPO_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_12673S'
butler = Butler(REPO_DIR)
```

```python
rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

expIds = [3020110200100,3020110200101,3020110200103,3020110200104,3020110200130,3020110200131,\
          3020110200106,3020110200107,3020110200133,3020110200134,3020110200109,3020110200110,\
          3020110200136,3020110200137,3020110200112,3020110200113,3020110200139,3020110200140,\
          3020110200115,3020110200116,3020110200142,3020110200143,3020110200145,3020110200146,\
          3020110200148,3020110200149,3020110200118,3020110200119,3020110200151,3020110200152,\
          3020110200154,3020110200155,3020110200157,3020110200158,3020110200121,3020110200122,\
          3020110200160,3020110200161,3020110200163,3020110200164,3020110200166,3020110200167,\
          3020110200169,3020110200170,3020110200124,3020110200125,3020110200172,3020110200173,\
          3020110200175,3020110200176,3020110200178,3020110200179,3020110200127,3020110200128,\
          3020110200181,3020110200182,3020110200184,3020110200185,3020110200187,3020110200188,\
          3020110200190,3020110200191,3020110200193,3020110200194,3020110200196,3020110200197,\
          3020110200199,3020110200200,3020110200202,3020110200203,3020110200205,3020110200206]


```

```python
maskArrays = {}
imageArrays = {}
for RAFT in ['R11']:#rafts:
    for SENSOR in ['S12']:#sensors:
        VENDOR, DETECTOR, plotNum = detector(RAFT,SENSOR)
        for expId in [3020110200190]:#expIds:
            postISRCCD = butler.get('postISRCCD',  raftName=RAFT,detectorName=SENSOR, expId=expId)
            amps = postISRCCD.getDetector().getAmplifiers()
            im1 = postISRCCD.getMaskedImage()
            for ampObject in amps:
                ampName = ampObject.getName()
                im1Area = im1[ampObject.getBBox()]
                maskArr = im1Area.getMask().getArray()
                maskArrays[ampName] = maskArr
                imageArr = im1Area.getImage().getArray()
                imageArrays[ampName] = imageArr
                w1 = np.where(maskArr == 0, 1, 0)
                print(ampName, np.sum(w1))

```

```python
plt.subplot(2,3,1)
plt.imshow(maskArrays['C06'])
plt.subplot(2,3,2)
plt.imshow(maskArrays['C06'][0:100,0:100])
plt.subplot(2,3,3)
plt.imshow(maskArrays['C06'][750:1000,100:350])

plt.subplot(2,3,4)
plt.imshow(maskArrays['C01'])
plt.subplot(2,3,5)
plt.imshow(maskArrays['C01'][0:100,0:100])
plt.subplot(2,3,6)
plt.imshow(maskArrays['C01'][750:1000,100:350])

```

```python
print(im1.getMask().getMaskPlaneDict().items())
```

```python
print('C06 Mask Top', maskArrays['C06'][2,6:11])
print('C06 Image Top', imageArrays['C06'][2,6:11])

print('C06 Mask Mid', maskArrays['C06'][50,6:11])
print('C06 Image Mid', imageArrays['C06'][50,6:11])

print('C06 Mask Defect', maskArrays['C06'][895,114:119])
print('C06 Image Defect', imageArrays['C06'][895,114:119])

print('C01 Mask Top', maskArrays['C01'][2,6:11])
print('C01 Image Top', imageArrays['C01'][2,6:11])

print('C01 Mask Mid', maskArrays['C01'][50,6:11])
print('C01 Image Mid', imageArrays['C01'][50,6:11])

print('C01 Mask Defect', maskArrays['C01'][895,114:119])
print('C01 Image Defect', imageArrays['C01'][895,114:119])
```

```python
print('C06 Mask Good Region', maskArrays['C06'][1000,114:119])
print('C06 Image Good Region', imageArrays['C06'][1000,114:119])
print('C01 Mask Good Region', maskArrays['C01'][1000,114:119])
print('C01 Image Good Region', imageArrays['C01'][1000,114:119])

```

```python
calibButler = Butler('/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_12673')
def2 = calibButler.get('defects', raftName=RAFT, detectorName=SENSOR, expId=expId)
print(def2.toSimpleTable())
```

```python
# This would query them all, but will take a long time.
nPixels = []
nPixelsDict = {}

for RAFT in rafts:
    for SENSOR in sensors:
        VENDOR, DETECTOR, plotNum = detector(RAFT,SENSOR)
        for expId in expIds:
            postISRCCD = butler.get('postISRCCD',  raftName=RAFT,detectorName=SENSOR, expId=expId)
            amps = postISRCCD.getDetector().getAmplifiers()
            im1 = postISRCCD.getMaskedImage()
            for ampObject in amps:
                ampName = ampObject.getName()
                im1Area = im1[ampObject.getBBox()]
                key = "%s_%s_%s"%(RAFT,SENSOR,ampName)
                w1 = np.where(im1Area.getMask().getArray() == 0, 1, 0)
                nPixels.append(w1)
                nPixelsDict[key] = w1
        print("Detector %d done."%DETECTOR)

```

```python
amps = postISRCCD.getDetector().getAmplifiers()
im1 = postISRCCD.getMaskedImage()
#print(dir(det))
for ampObject in amps:
    #print(dir(ampObject))
    ampName = ampObject.getName()
    im1Area = im1[ampObject.getBBox()]
    w1 = np.where(im1Area.getMask().getArray() == 0, 1, 0)
    print(ampName, np.sum(w1))
```

```python

```
