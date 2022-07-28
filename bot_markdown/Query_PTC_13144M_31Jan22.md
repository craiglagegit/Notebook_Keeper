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

Initially written 27 Jan 2022 by Craig Lage

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.io.fits as pf
from lsst.daf.butler import Butler
import lsst.afw.math as afwMath
from focal_plane_plotting import plot_focal_plane
```

```python
!eups list -s | grep ip_isr
!eups list -s | grep cp_pipe
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib"])
camera = butler.get('camera', instrument='LSSTCam')
```

```python
ptcButler = Butler("/repo/main", collections=["u/cslage/bps_13144M"])
exposure=3021120600576
```

```python
# Get the eotest results
filename = "/project/cslage/BOT_LSSTCam/eotest/eotest_gain_13144_15dec21.pkl"
file = open(filename, 'rb')
#fe55_results = pkl.load(file)
ptc_results = pkl.load(file)
file.close()

rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

```

```python
def getDetector(raft, sensor):
    # Subroutine to find vendor and detector number given raft and sensor                                                                                                                                                           
    startingCol = [1,0,0,0,1] # First raft column in each row                                                                                                                                                                       
    rows = [0,3,8,13,18] # Starting raft sequence number of each row                                                                                                                                                                
    if raft in ['R11','R12','R13','R14','R21','R22','R23','R24','R30',\
                'R31','R32','R33','R34']:
        vendor = 'E2V'
    else:
        vendor = 'ITL'
    raftRow = int(list(raft)[1])
    raftCol = int(list(raft)[2])
    sensorRow = int(list(sensor)[1])
    sensorCol = int(list(sensor)[2])
    detectorNum = (rows[raftRow] + (raftCol - startingCol[raftRow])) * 9
    detectorNum += 3 * sensorRow + sensorCol
    return vendor, detectorNum

# This dictionary captures the amp naming correspondence
slacAmps = {'C10':'AMP01','C11':'AMP02','C12':'AMP03','C13':'AMP04',\
           'C14':'AMP05','C15':'AMP06','C16':'AMP07','C17':'AMP08',\
           'C07':'AMP09','C06':'AMP10','C05':'AMP11','C04':'AMP12',\
           'C03':'AMP13','C02':'AMP14','C01':'AMP15','C00':'AMP16'}
```

```python tags=[]

for ii, [RAFT,SENSOR,amp] in enumerate([['R02', 'S00', 'C10']]):
    slacAmp = slacAmps[amp]
    slacNum = int(slacAmp.strip('AMP')) - 1
    VENDOR, DETECTOR  = getDetector(RAFT, SENSOR)
    ptcDataset = ptcButler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')
    rawMeans = ptcDataset.rawMeans[amp]
    rawVars = ptcDataset.rawVars[amp]
    filename = "/project/cslage/BOT_LSSTCam/eotest/%s_%s_13144_ptc.fits"%(RAFT,SENSOR)
    hdu = pf.open(filename)
    slacData = hdu[1].data
    slacMeans = slacData['%s_MEAN'%slacAmp]
    slacVars = slacData['%s_VAR'%slacAmp]
    sortedList = sorted(zip(slacVars, slacMeans))
    slacMeans = [x[1] for x in sortedList]
    slacVars = [x[0] for x in sortedList]
    print(f"Detector: {DETECTOR}, amp: {amp}")
    for i in range(len(slacMeans)):
        print(f"Index={i}, SLAC mean = {slacMeans[i]:.2f}, SLAC Var = {slacVars[i]:.2f}")

```

```python tags=[]
slacMeans
```

```python
test = [x[0] for x in sorted(zip(slacVars, slacMeans))]
```

```python
sortedList = sorted(zip(slacVars, slacMeans))
slacMeans = [x[1] for x in sortedList]
slacVars = [x[0] for x in sortedList]
```

```python
type(test[7])
```

```python
exposure=3021120600576
linButler = Butler("/repo/main", collections=["u/cslage/linearizer_28jan22"])
corr = linButler.get('pdCorrection', exposure=exposure, instrument='LSSTCam')
```

```python
keyList = []
corrList = []
for i, key in enumerate(corr.abscissaCorrections.keys()):
    keyList.append(key)
    corrList.append(corr.abscissaCorrections[key])
    if i > 4:
        break
```

```python
keyList
```

```python
corrList
```

```python

```

```python jupyter={"outputs_hidden": true} tags=[]
dir(butler)
```

```python jupyter={"outputs_hidden": true} tags=[]
dir(butler.registry)
```

```python
types = butler.registry.queryDatasetTypes()
```

```python jupyter={"outputs_hidden": true} tags=[]
for type in types:
    print(type)
```

```python
butler = Butler("/repo/main", collections=["LSSTComCam/raw/all","LSSTComCam/calib"])
types = butler.registry.queryDatasetTypes()
for type in types:
    print(type)
```

```python

```
