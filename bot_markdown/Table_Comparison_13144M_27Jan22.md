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
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib", "u/cslage/bps_13144M", "u/cslage/linearizerB_26jan22"])
camera = butler.get('camera', instrument='LSSTCam')
exposure=3021120600576
```

```python
# Get the eotest results
filename = "/project/cslage/BOT_LSSTCam/eotest/eotest_gain_13144_15dec21.pkl"
file = open(filename, 'rb')
#fe55_results = pkl.load(file)
ptc_results = pkl.load(file)
file.close()

# This dictionary captures the amp naming correspondence
slacAmps = {'C10':'AMP01','C11':'AMP02','C12':'AMP03','C13':'AMP04',\
           'C14':'AMP05','C15':'AMP06','C16':'AMP07','C17':'AMP08',\
           'C07':'AMP09','C06':'AMP10','C05':'AMP11','C04':'AMP12',\
           'C03':'AMP13','C02':'AMP14','C01':'AMP15','C00':'AMP16'}

```

```python tags=[]
filename = "/repo/main/u/cslage/bps_13144M/plots/ptc_table_13144M_27jan22.txt"
file = open(filename, 'w')
header = "Amp\t Gain_DM\t Gain_EO\t A00_DM\t A00_EO\t Noise_DM\t Noise_EO\t Turn_DM\t Turn_EO\t MaxNL_DM\t CorrStd\n"
file.write(header)
fluxMin = 10000.0

for detector in camera:
    if detector.getType().name != 'SCIENCE':
        continue
    detName = detector.getName()
    RAFT = detName.split('_')[0]
    SENSOR = detName.split('_')[1]
    DETECTOR = detector.getId()
    try:
        ptc = butler.get('ptc', detector=DETECTOR, instrument='LSSTCam')
        lin = butler.get('linearizer', detector=DETECTOR, instrument='LSSTCam')
    except:
        continue
    eoPTCGain = ptc_results['ptc_gain'][RAFT][SENSOR]
    eoPtcTurnoff = ptc_results['ptc_turnoff'][RAFT][SENSOR]
    eoA00 = ptc_results['ptc_a00'][RAFT][SENSOR]
    eoNoise = ptc_results['ptc_noise'][RAFT][SENSOR]
    for amp in detector.getAmplifiers():
        ampName = amp.getName()
        slacAmp = slacAmps[ampName]
        slacNum = int(slacAmp.strip('AMP')) - 1
        if ptc.ptcFitType == 'EXPAPPROXIMATION':
            dmA00 = ptc.ptcFitPars[ampName][0]
        if ptc.ptcFitType == 'FULLCOVARIANCE':
            dmA00 = ptc.aMatrix[ampName][0][0]
        dmMeans = np.array(ptc.finalMeans[ampName])
        dmMeans = dmMeans[~np.isnan(dmMeans)]
        if len(dmMeans > 0):
            maxDM = dmMeans.max()
        else:
            maxDM = 0.0            
        centers, values = np.split(lin.linearityCoeffs[ampName], 2)
        fluxMask = np.where(centers>fluxMin)
        try:
            maxDeviation = np.max(abs((values/centers * 100.0)[fluxMask]))
        except:
            maxDeviation = np.nan
        mask = np.array(ptc.expIdMask[ampName], dtype=bool)
        means = np.array(ptc.rawMeans[ampName])[mask]
        corrResiduals = np.array(lin.fitResiduals[ampName])[mask]
        fluxMask = means > fluxMin
        try:
            corrStd = np.nanstd((corrResiduals/means * 100.0)[fluxMask])
        except:
            corrStd = np.nan
        data = f"{RAFT}_{SENSOR}_{ampName}\t {ptc.gain[ampName]:.6f}\t {eoPTCGain[slacNum]:.6f}\
        \t {dmA00:.6g}\t {-eoA00[slacNum]:.6g}\t {ptc.noise[ampName]:.2f}\t {eoNoise[slacNum]:.2f}\
        \t {maxDM:.2f}\t {eoPtcTurnoff[slacNum]:.2f}\t {maxDeviation:.4f}\t {corrStd:.6f}\n"
        file.write(data)
        #break
    #break
file.close()
```

```python

```
