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

# Notebook for extracting PTC and linearity data

Initially written 08 Jun 22 by Craig Lage

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.butler import Butler
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib",\
                                                    "u/cslage/calib/13144/calib.20220103"])
```

```python
nonlinPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144M"])
```

```python
def ExpApprox(mu, g, a00, n):
    if (g < 1.0E-6) or (abs(a00) < 1.0E-9):
        return np.zeros([len(mu)])
    else:
        expFactor = 2.0 * a00 * mu * g
        if max(expFactor) > 100.0:
            return np.zeros([len(mu)])
        else:
            preFactor = 1.0 / (2.0 * g * g * a00)
            noiseTerm = n / (g * g)
            return preFactor * (np.exp(expFactor) - 1.0) + noiseTerm
        
def calcMondiode(expId):
    factor = 5.0
    DATA_DIR = '/lsstdata/offline/teststand/BOT/storage/'
    date = int(expId/100000)
    seq = expId - date * 100000
    date = date - 10000000
    file = DATA_DIR + '%d/MC_C_%d_%06d/Photodiode_Readings_%d_%06d.txt'%(date,date,seq,date,seq)

    x, y = np.recfromtxt(file).transpose()
    # Threshold for finding baseline current values:                                                                                                                                                         
    ythresh = (min(y) + max(y))/factor + min(y)
    # Subtract the median of the baseline values to get a calibrated                                                                                                                                         
    # current.                                                                                                                                                                                               
    y -= np.median(y[np.where(y < ythresh)])
    integral = sum((y[1:] + y[:-1])/2*(x[1:] - x[:-1]))
    return integral

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
    raftCol = int(list(raft)[2])
    sensorRow = int(list(sensor)[1])
    sensorCol = int(list(sensor)[2])
    detectorNum = (rows[raftRow] + (raftCol - startingCol[raftRow])) * 9
    detectorNum += 3 * sensorRow + sensorCol
    return vendor, detectorNum, 4 - raftRow, raftCol

rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

```

```python tags=[]
expId=3021120700200

ampList = [['R03', 'S11', 'C06'], ['R22', 'S11', 'C06']]

for [R, S, C] in ampList:
    file = open(f"/project/cslage/BOT_LSSTCam/gen3/PTC_{R}_{S}_{C}.txt", "w")
    file.write(f"PTC data {R}_{S}_{C}\n")
    file.write("Exposure pair \t\t\t Mean \t Variance \t Exposure Time (sec) \t Monitor Diode (Coulombs)\n")
    for RAFT in [R]:#rafts:
        for SENSOR in [S]:#sensors:
            VENDOR, det, raftRow, raftCol = detector(RAFT, SENSOR)
            ptc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
            break
            #lin = linPtcButler.get('linearizer', detector=det, exposure=expId, instrument='LSSTCam')

            for ampName in ptc.gain.keys():
                if ampName != C:
                    continue
                gain = ptc.gain[ampName]
                a00 = ptc.ptcFitPars[ampName][0]
                noise = ptc.noise[ampName]
                mask = np.array(ptc.expIdMask[ampName], dtype=bool)
                maxDM = np.max(np.array(ptc.rawMeans[ampName])[mask])

                modExpTimes = []
                for ii, pair in enumerate(ptc.inputExpIdPairs[ampName]):
                    pair = pair[0]
                    modExpTime = 0.0
                    nExps = 0
                    for j in range(2):
                        expId = pair[j]
                        try:
                            monDiode = calcMondiode(expId)
                            modExpTime += monDiode
                            nExps += 1
                        except:
                            continue
                    if nExps > 0:
                        myMonDiode = modExpTime / nExps
                    else:
                        mask[ii] = False
                    file.write(f"{pair} \t {ptc.rawMeans[ampName][ii]} \t {ptc.rawVars[ampName][ii]} \t {ptc.rawExpTimes[ampName][ii]} \t  {myMonDiode} \n")
    break
    file.close()


    
```

```python
dir(ptc)
```

```python

```
