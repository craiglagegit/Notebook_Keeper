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

# Notebook for extracting PTC data.

Initially written 15 Jun 2022 by Craig Lage.

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.io.fits as pf
from lsst.daf.butler import Butler
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib", "u/cslage/bps_13144M"])
exposure=3021120600576
```

```python
ampList = [['R03', 'S11', 'C06']]

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

```

```python
def ExpApprox(mu, g, a00, n):
    # This is the equation for the ExpApprox for use in the notebook
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
```

```python tags=[]
for ii, [RAFT,SENSOR,amp] in enumerate(ampList):
    VENDOR, DETECTOR  = getDetector(RAFT, SENSOR)
    fig = plt.figure(figsize=(16,16))
    ax1 = plt.axes([0.1,0.1,0.8,0.8])
    ptcDataset = butler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')
    dmGain = ptcDataset.gain[amp]
    if ptcDataset.ptcFitType == 'EXPAPPROXIMATION':
        dmA00 = ptcDataset.ptcFitPars[amp][0]
    if ptcDataset.ptcFitType == 'FULLCOVARIANCE':
        dmA00 = ptcDataset.aMatrix[amp][0][0]
    dmNoise = ptcDataset.noise[amp]
    rawMeans = ptcDataset.rawMeans[amp]
    rawVars = ptcDataset.rawVars[amp]
    dmMeans = np.array(ptcDataset.finalMeans[amp])
    dmMeans = dmMeans[~np.isnan(dmMeans)]
    if len(dmMeans > 0):
        maxDM = dmMeans.max()
    else:
        maxDM = 0.0
    
    ax1.set_title("PTC, %s_%s_%s_Det_%d_%s"%(RAFT,SENSOR, amp, DETECTOR, VENDOR), fontsize = 18)
    ax1.scatter(rawMeans, rawVars, marker = 'x', s=200, color = 'red', label = 'Data')
    
    ax1.text(10000,40000,"Gain = %.4f"%dmGain, fontsize=18)
    ax1.text(10000,38000,"Noise = %.4f"%dmNoise, fontsize=18)
    ax1.text(10000,36000,"A00 = %.6g"%dmA00, fontsize=18)
    ax1.text(10000,34000,"Max ADU = %.1f"%maxDM, fontsize=18)
    xplot = np.linspace(0.0, 120000, 200)
    yplot = ExpApprox(xplot,dmGain,dmA00,dmNoise)
    linyplot = xplot * 1/dmGain
    ax1.plot(xplot, yplot, ls = '--', color = 'red', label = "ExpApprox")
    ax1.plot(xplot, linyplot, ls = '--', color = 'green', label = "Linear Fit")
    ax1.set_xlabel("Mean(ADU)", fontsize=18)
    ax1.set_ylabel("Variance(ADU)", fontsize=18)
    ax1.set_xlim(0,120000)
    ax1.set_ylim(0,50000)
    ax1.legend(fontsize = 18)
    file = open("/project/cslage/BOT_LSSTCam/gen3/PTC_R03_S11_C06_Det22.txt", "w")
    file.write("PTC data Detect 22, R03_S11_C06\n")
    file.write("Mean(ADU)                Variance(ADU)\n")
    for i in range(len(rawMeans)):
        file.write(f"{rawMeans[i]}         {rawVars[i]}\n")
    file.close()

plt.savefig("/project/cslage/BOT_LSSTCam/gen3/PTC_R03_S11_C06_Det22.png")

```

```python

```
