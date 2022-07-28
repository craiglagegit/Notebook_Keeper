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

# Notebook for comparing eotest gain with DM gain.

Initially written 20 Nov 2021 by Craig Lage.

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.butler import Butler
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib", "u/cslage/bps_13144M", "u/cslage/linearizer_20220104"])
exposure=3021120600576
```

```python
# Get the eotest results
filename = "/project/cslage/BOT_LSSTCam/eotest/eotest_gain_13144_15dec21.pkl"
file = open(filename, 'rb')
#fe55_results = pkl.load(file)
ptc_results = pkl.load(file)
file.close()
print(ptc_results.keys())

rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']
print(rafts)
print(ptc_results['ptc_gain']['R01']['S00'])
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
    raftCol = int(list(raft)[2])
    sensorRow = int(list(sensor)[1])
    sensorCol = int(list(sensor)[2])
    detectorNum = (rows[raftRow] + (raftCol - startingCol[raftRow])) * 9
    detectorNum += 3 * sensorRow + sensorCol
    return vendor, detectorNum, 4 - raftRow, raftCol

# This dictionary captures the amp naming correspondence
slacAmps = {'C10':'AMP01','C11':'AMP02','C12':'AMP03','C13':'AMP04',\
           'C14':'AMP05','C15':'AMP06','C16':'AMP07','C17':'AMP08',\
           'C07':'AMP09','C06':'AMP10','C05':'AMP11','C04':'AMP12',\
           'C03':'AMP13','C02':'AMP14','C01':'AMP15','C00':'AMP16'}

```

```python tags=[]
filename = "/repo/main/u/cslage/bps_13144M/plots/ptc_table_13144M_04jan22.txt"
file = open(filename, 'w')
header = "Amp\t\tGain_DM\t\tGain_EO\t\tA00_DM\t\tA00_EO\t     Noise_DM  Noise_EO  Turn_DM\tTurn_EO\t\tMaxNL_DM\n"
file.write(header)

for RAFT in rafts:
    for SENSOR in sensors:
        VENDOR, DETECTOR, raftRow, raftCol = detector(RAFT, SENSOR)
        try:
            ptcDataset = butler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')
        except:
            continue
        eoPTCGain = ptc_results['ptc_gain'][RAFT][SENSOR]
        eoPtcTurnoff = ptc_results['ptc_turnoff'][RAFT][SENSOR]
        eoA00 = ptc_results['ptc_a00'][RAFT][SENSOR]
        eoNoise = ptc_results['ptc_noise'][RAFT][SENSOR]
        for amp in ptcDataset.gain.keys():
            slacAmp = slacAmps[amp]
            slacNum = int(slacAmp.strip('AMP')) - 1
            if ptcDataset.ptcFitType == 'EXPAPPROXIMATION':
                dmA00 = ptcDataset.ptcFitPars[amp][0]
            if ptcDataset.ptcFitType == 'FULLCOVARIANCE':
                dmA00 = ptcDataset.aMatrix[amp][0][0]
            dmMeans = np.array(ptcDataset.finalMeans[amp])
            dmMeans = dmMeans[~np.isnan(dmMeans)]
            if len(dmMeans > 0):
                maxDM = dmMeans.max()
            else:
                maxDM = 0.0
            try:
                lin = butler.get('linearizer', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')
                centers, values = np.split(lin.linearityCoeffs[amp], 2)
                fluxMask = np.where(centers>20000.0)
                maxDeviation = np.max(abs((values/centers * 100.0)[fluxMask]))
            except:
                maxDeviation = np.nan
            
            data = f"{RAFT}_{SENSOR}_{amp}\t{ptcDataset.gain[amp]:.6f}\t{eoPTCGain[slacNum]:.6f}\t{dmA00:.6g}\t{-eoA00[slacNum]:.6g}\t{ptcDataset.noise[amp]:.2f}\t{eoNoise[slacNum]:.2f}\t{maxDM:.2f}\t{eoPtcTurnoff[slacNum]:.2f}\t{maxDeviation:.4f}\n"
            file.write(data)
file.close()
```

```python

```
