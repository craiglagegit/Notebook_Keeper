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

# Notebook for plotting long-range serial correlations.

Initially written 12 Nov 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
! eups list -s cp_pipe
! eups list -s ip_isr
```

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
from lsst.ip.isr import PhotonTransferCurveDataset
```

```python
run = '12673'
REPO_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_FullCov_%sA'%run
DATA_DIR = '/project/shared/BOT/'
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


rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

```

```python jupyter={"outputs_hidden": true}
covariances = {}
means = {}
detectors = []
for RAFT in rafts:
    for SENSOR in sensors:
        if (RAFT=='R32' and SENSOR in ['S00','S01','S02']) or (RAFT=='R33' and SENSOR in ['S20','S21','S22']):
            continue

        VENDOR, DETECTOR, raftRow, raftCol = detector(RAFT, SENSOR)
        if VENDOR == 'ITL':
            continue
        try:
            datasetFile = REPO_DIR+'/calibrations/ptc/ptcDataset-det%03d.fits'%DETECTOR
            ptcDataset = PhotonTransferCurveDataset.readFits(datasetFile)
            detectors.append(DETECTOR)
            covar = ptcDataset.covariances
            covariances[DETECTOR] = covar
            mus = ptcDataset.finalMeans
            means[DETECTOR] = mus
            print("Found detector %d"%DETECTOR, ptcDataset.ptcFitType)
        except:
            print("Didn't find detector %d"%DETECTOR)
            continue

        #for amp in cov.keys():
        #    covariances[DETECTOR].append(cov[amp])

print(covar.keys())
```

```python
# Plot the serial covariances vs distance 
PlotDelta = 8
jj = 1
#plt.figure(figsize=(16,16))
fig, [ax1, ax2] = plt.subplots(ncols=1, nrows=2, figsize=(16,16))
plt.suptitle("SLAC run 12673/12674 - Long Range Serial Covariances - E2V",fontsize = 24)
goodPoints = 0
badPoints = 0
badAmps = []
for RAFT in rafts:
    for SENSOR in sensors:
        if (RAFT=='R32' and SENSOR in ['S00','S01','S02']) or (RAFT=='R33' and SENSOR in ['S20','S21','S22']):
            continue

        VENDOR, DETECTOR, raftRow, raftCol = detector(RAFT, SENSOR)
        det = DETECTOR
        if VENDOR == 'ITL':
            continue
        for amp in covar.keys():
            highPoints = 0
            NumFluxes = int(len(means[det][amp]))
            xaxis = []
            yaxis = []
            for ii in range(1,PlotDelta):
                cov = []
                flux = []
                for n in range(NumFluxes):
                    xcorr = covariances[det][amp][n][ii][jj]
                    mean = means[det][amp][n]
                    if np.isnan(xcorr):
                        continue
                    #print(mean, xcorr)
                    cov.append(xcorr)
                    flux.append(mean)
                cov = np.array(cov)
                flux = np.array(flux)
                if len(flux) == 0:
                    badPoints += 1
                    continue
                coefs = np.polyfit(flux*flux, cov, 1)
                #print(det, amp, ii, coefs[0])
                if coefs[0] > 0:
                    goodPoints += 1
                    if coefs[0] > 10**(-7.0):
                        highPoints += 1
                    xaxis.append(ii)
                    yaxis.append(coefs[0])
                else:
                    badPoints += 1
                    #print(det,amp, flux, cov)
            
            ax1.plot(xaxis, np.log10(yaxis), label=amp)
            if highPoints > 2:
                badAmps.append("%s_%s_%s"%(RAFT,SENSOR,amp))
                ax2.plot(xaxis, np.log10(yaxis), label="%s_%s_%s"%(RAFT,SENSOR,amp))
ax1.set_title("All amps", fontsize = 18)
ax1.set_ylim(-9.0,-6.0)
ax1.set_ylabel("Log Covariance",fontsize=18)
ax1.set_xlabel("X pixel", fontsize=18)
ax2.set_title("Worst amps", fontsize = 18)
ax2.set_ylim(-9.0,-6.0)
ax2.set_ylabel("Log Covariance",fontsize=18)
ax2.set_xlabel("X pixel", fontsize=18)
ax2.legend()
plt.savefig(REPO_DIR+'/plots/Cov_vs_X_Log_13Nov20.pdf')

print("%d Good Points, %d Bad Points"%(goodPoints, badPoints))
```

```python
# Plot the serial covariances vs distance 
PlotDelta = 8
jj = 1
#plt.figure(figsize=(16,16))
fig, ax1 = plt.subplots(ncols=1, nrows=1, figsize=(16,8))
plt.suptitle("SLAC run 12673/12674 - Long Range Serial Covariances - R21_S11",fontsize = 24)
goodPoints = 0
badPoints = 0
badAmps = []
for RAFT in rafts:
    for SENSOR in sensors:
        if RAFT != 'R21' or SENSOR != 'S11':
            continue

        VENDOR, DETECTOR, raftRow, raftCol = detector(RAFT, SENSOR)
        det = DETECTOR
        if VENDOR == 'ITL':
            continue
        for amp in covar.keys():
            highPoints = 0
            NumFluxes = int(len(means[det][amp]))
            xaxis = []
            yaxis = []
            for ii in range(1,PlotDelta):
                cov = []
                flux = []
                for n in range(NumFluxes):
                    xcorr = covariances[det][amp][n][ii][jj]
                    mean = means[det][amp][n]
                    if np.isnan(xcorr):
                        continue
                    #print(mean, xcorr)
                    cov.append(xcorr)
                    flux.append(mean)
                cov = np.array(cov)
                flux = np.array(flux)
                if len(flux) == 0:
                    badPoints += 1
                    continue
                coefs = np.polyfit(flux*flux, cov, 1)
                #print(det, amp, ii, coefs[0])
                if coefs[0] > 0:
                    goodPoints += 1
                    if coefs[0] > 10**(-7.0):
                        highPoints += 1
                    xaxis.append(ii)
                    yaxis.append(coefs[0])
                else:
                    badPoints += 1
                    #print(det,amp, flux, cov)
            
            ax1.plot(xaxis, np.log10(yaxis), label=amp)
ax1.set_title("All amps", fontsize = 18)
ax1.set_ylim(-9.0,-6.0)
ax1.set_ylabel("Log Covariance",fontsize=18)
ax1.set_xlabel("X pixel", fontsize=18)
ax1.legend()
plt.savefig(REPO_DIR+'/plots/Cov_vs_X_Log_13Nov20_85.pdf')

print("%d Good Points, %d Bad Points"%(goodPoints, badPoints))
```

```python

```

```python

```
