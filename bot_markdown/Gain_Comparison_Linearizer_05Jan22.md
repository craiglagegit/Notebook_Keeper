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

# Notebook for comparing the impact of the linearizer on gain and other parameters.

Initially written 05 Jan 2022 by Craig Lage.

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.butler import Butler
```

```python
exposure=3021120600576
linPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144N"])
nonlinPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144M"])
```

```python
rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']



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
```

```python tags=[]
Lin_gains = {}
Non_Lin_gains = {}
Lin_noise = {}
Non_Lin_noise = {}
Lin_a00 = {}
Non_Lin_a00 = {}
Lin_turnoff = {}
Non_Lin_turnoff = {}

detectors = {}
for RAFT in rafts:
    for SENSOR in sensors:
        VENDOR, DETECTOR, raftRow, raftCol = detector(RAFT, SENSOR)
        try:
            linPTC = linPtcButler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')
            nonlinPTC = nonlinPtcButler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')            
        except:
            continue

        detectors[DETECTOR] = [VENDOR, RAFT, SENSOR]
        Lin_gains[DETECTOR] = []
        Non_Lin_gains[DETECTOR] = []
        Lin_noise[DETECTOR] = []
        Non_Lin_noise[DETECTOR] = []
        Lin_a00[DETECTOR] = []
        Non_Lin_a00[DETECTOR] = []
        Lin_turnoff[DETECTOR] = []
        Non_Lin_turnoff[DETECTOR] = []

        for amp in linPTC.gain.keys():
            if linPTC.ptcFitType == 'EXPAPPROXIMATION':
                linA00 = linPTC.ptcFitPars[amp][0]
            if linPTC.ptcFitType == 'FULLCOVARIANCE':
                linA00 = linPTC.aMatrix[amp][0][0]
            if nonlinPTC.ptcFitType == 'EXPAPPROXIMATION':
                nonlinA00 = nonlinPTC.ptcFitPars[amp][0]
            if nonlinPTC.ptcFitType == 'FULLCOVARIANCE':
                nonlinA00 = nonlinPTC.aMatrix[amp][0][0]
                
            linMeans = np.array(linPTC.finalMeans[amp])
            linMeans = linMeans[~np.isnan(linMeans)]
            if len(linMeans > 0):
                linmaxDM = linMeans.max()
            else:
                linmaxDM = 0.0
                
            nonlinMeans = np.array(nonlinPTC.finalMeans[amp])
            nonlinMeans = nonlinMeans[~np.isnan(nonlinMeans)]
            if len(nonlinMeans > 0):
                nonlinmaxDM = nonlinMeans.max()
            else:
                nonlinmaxDM = 0.0
            Lin_gains[DETECTOR].append(linPTC.gain[amp])
            Non_Lin_gains[DETECTOR].append(nonlinPTC.gain[amp])
            Lin_noise[DETECTOR].append(linPTC.noise[amp])
            Non_Lin_noise[DETECTOR].append(nonlinPTC.noise[amp])
            Lin_a00[DETECTOR].append(linA00)
            Non_Lin_a00[DETECTOR].append(nonlinA00)
            Lin_turnoff[DETECTOR].append(linmaxDM)
            Non_Lin_turnoff[DETECTOR].append(nonlinmaxDM)

```

```python
len(Lin_gains)
```

```python
minGain = 1.2
maxGain = 2.0
plotCounter = 0
numCCDs = len(detectors)
numAmps = 0
num = 343
plt.figure(figsize=(16,8))

ratios = []
plotCounter += 1
plt.subplot(1,2,plotCounter)
plt.title("Gain Comparison", fontsize = 18)
for [vendor, color] in [['E2V', 'red'], ['ITL', 'green']]:
    xplot = []
    yplot = []
    for DETECTOR in detectors.keys():
        [VENDOR, RAFT, SENSOR] = detectors[DETECTOR]
        for i in range(16):
            if Non_Lin_gains[DETECTOR][i] >.001:
                ratio = Lin_gains[DETECTOR][i] / Non_Lin_gains[DETECTOR][i]
            else:
                ratio = 1.0
            ratios.append(ratio)
            if vendor == VENDOR:
                numAmps += 1
                xplot.append(Non_Lin_gains[DETECTOR][i])
                yplot.append(Lin_gains[DETECTOR][i])
    plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
plt.xlim(minGain, maxGain)
plt.ylim(minGain, maxGain)
plt.xlabel("NonLinearized Gains(e-/ADU) (Run13144M) ", fontsize = 18)
plt.ylabel("Linearized Gains(e-/ADU) (Run 13144N)", fontsize = 18)
xplot = np.linspace(minGain, maxGain,100)
plt.plot(xplot, xplot, ls = '--', color='blue')
slope1 = 1.05
plt.plot(xplot,slope1*xplot, ls = '--', color='red', label = '%.3f'%slope1)
slope2 = 0.95
plt.plot(xplot,slope2*xplot, ls = '--', color='red', label = '%.3f'%slope2)

plt.text(1.3, 1.85, "%d CCDs"%numCCDs, fontsize = 18)
plt.text(1.3, 1.80, "%d amps"%numAmps, fontsize = 18)
plt.legend(loc='lower right', fontsize = 18)
plt.subplot(1,2,plotCounter+1)
plt.title("Gain Ratio (Lin / Nonlin Run 13144)", fontsize = 18)
n, bins, patches = plt.hist(ratios, bins = 50, range=(0.90,1.10))
ymax = n.max() * 1.10
plt.xlim(0.90, 1.10)
plt.ylim(0, ymax)
plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
plt.text(1.01, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
plt.text(1.01, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
plt.text(1.01, 0.85*n.max(), "%d amps"%numAmps, fontsize = 18)
plt.savefig("/repo/main/u/cslage/bps_13144N/plots/Lin_Nonlin_Gains_13144N_05Jan22.pdf")

```

```python
minPlot = 1.0
maxPlot = 5.0
plotCounter = 0
numCCDs = len(detectors)
numAmps = 0
num = 343
plt.figure(figsize=(16,8))

ratios = []
plotCounter += 1
plt.subplot(1,2,plotCounter)
plt.title("A00 Comparison", fontsize = 18)
badx = []
bady = []
for [vendor, color] in [['E2V', 'red'], ['ITL', 'green']]:
    xplot = []
    yplot = []
    for DETECTOR in detectors.keys():
        [VENDOR, RAFT, SENSOR] = detectors[DETECTOR]
        for i in range(16):
            if vendor == VENDOR:
                numAmps += 1
                non_a00 = -Non_Lin_a00[DETECTOR][i] * 1.0E6
                if non_a00 > minPlot and non_a00 < maxPlot:

                    ratio = Lin_a00[DETECTOR][i] / Non_Lin_a00[DETECTOR][i]
                    ratios.append(ratio)
                xplot.append(-Non_Lin_a00[DETECTOR][i] * 1.0E6)
                yplot.append(-Lin_a00[DETECTOR][i] * 1.0E6)
    plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
plt.xlim(minPlot, maxPlot)
plt.ylim(minPlot, maxPlot)
plt.xlabel("NonLinearized A00 (Run13144M)", fontsize = 18)
plt.ylabel("Linearized A00 (Run 13144N)", fontsize = 18)
xplot = np.linspace(minPlot, maxPlot,100)
plt.plot(xplot, xplot, ls = '--', color='blue')
slope1 = 1.10
plt.plot(xplot,slope1*xplot, ls = '--', color='red', label = '%.3f'%slope1)
slope2 = 0.90
plt.plot(xplot,slope2*xplot, ls = '--', color='red', label = '%.3f'%slope2)

#plt.text(1.05, 1.55, "%d CCDs"%numCCDs, fontsize = 18)
#plt.text(1.05, 1.525, "%d amps"%numAmps, fontsize = 18)
plt.legend(loc='lower right', fontsize = 18)
plt.subplot(1,2,plotCounter+1)
plt.title("A00 Ratio (Lin / NonLin Run 13144)", fontsize = 18)
n, bins, patches = plt.hist(ratios, bins = 50, range=(0.80,1.20))
ymax = n.max() * 1.10
plt.xlim(0.80, 1.20)
plt.ylim(0, ymax)
plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
plt.text(1.03, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
plt.text(1.03, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
plt.text(1.03, 0.85*n.max(), "%d amps"%len(ratios), fontsize = 18)
plt.savefig("/repo/main/u/cslage/bps_13144N/plots/Lin_NonLin_A00_13144N_05Jan22.pdf")
```

```python
minPlot = 1.0
maxPlot = 20.0
plotCounter = 0
numCCDs = len(detectors)
numAmps = 0
num =343
plt.figure(figsize=(16,8))
ratios = []
plotCounter += 1
plt.subplot(1,2,plotCounter)
plt.title("Noise Comparison", fontsize = 18)
badx = []
bady = []
for [vendor, color] in [['E2V', 'red'], ['ITL', 'green']]:
    xplot = []
    yplot = []
    for DETECTOR in detectors.keys():
        [VENDOR, RAFT, SENSOR] = detectors[DETECTOR]
        for i in range(16):
            if vendor == VENDOR:
                numAmps += 1
                if Non_Lin_noise[DETECTOR][i] > minPlot and Non_Lin_noise[DETECTOR][i] < maxPlot:
                    ratio = Lin_noise[DETECTOR][i] / Non_Lin_noise[DETECTOR][i]
                    ratios.append(ratio)

                xplot.append(Non_Lin_noise[DETECTOR][i])
                yplot.append(Lin_noise[DETECTOR][i])
    plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
plt.xlim(minPlot, maxPlot)
plt.ylim(minPlot, maxPlot)
plt.xlabel("NonLinearized Noise (e-) Run 13144M ", fontsize = 18)
plt.ylabel("Linearized Noise (e-) Run 13144N", fontsize = 18)
xplot = np.linspace(minPlot, maxPlot,100)
plt.plot(xplot, xplot, ls = '--', color='blue')
slope1 = 1.10
plt.plot(xplot,slope1*xplot, ls = '--', color='red', label = '%.3f'%slope1)
slope2 = 0.90
plt.plot(xplot,slope2*xplot, ls = '--', color='red', label = '%.3f'%slope2)

#plt.text(1.05, 1.55, "%d CCDs"%numCCDs, fontsize = 18)
#plt.text(1.05, 1.525, "%d amps"%numAmps, fontsize = 18)
plt.legend(loc='lower right', fontsize = 18)
plt.subplot(1,2,plotCounter+1)
plt.title("Noise Ratio (Lin / NonLin)", fontsize = 18)
n, bins, patches = plt.hist(ratios, bins = 50, range=(0.80,1.40))
ymax = n.max() * 1.10
plt.xlim(0.80, 1.40)
plt.ylim(0, ymax)
plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
plt.text(1.03, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
plt.text(1.03, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
plt.text(1.03, 0.85*n.max(), "%d amps"%len(ratios), fontsize = 18)
plt.savefig("/repo/main/u/cslage/bps_13144N/plots/Lin_NonLin_Noise_13144N_05Jan22.pdf")

```

```python
minPlot = 10000.0
maxPlot = 200000.0
plotCounter = 0
numCCDs = len(detectors)
numAmps = 0
num =343
plt.figure(figsize=(16,8))

ratios = []
plotCounter += 1
plt.subplot(1,2,plotCounter)
plt.title("Turnoff Comparison", fontsize = 18)
badx = []
bady = []
for [vendor, color] in [['E2V', 'red'], ['ITL', 'green']]:
    xplot = []
    yplot = []
    for DETECTOR in detectors.keys():
        [VENDOR, RAFT, SENSOR] = detectors[DETECTOR]
        for i in range(16):
            if vendor == VENDOR:
                numAmps += 1
                eo_turnoff = Non_Lin_turnoff[DETECTOR][i]
                if eo_turnoff > minPlot and eo_turnoff < maxPlot:
                    ratio = Lin_turnoff[DETECTOR][i] / Non_Lin_turnoff[DETECTOR][i]
                    ratios.append(ratio)

                xplot.append(Non_Lin_turnoff[DETECTOR][i])
                yplot.append(Lin_turnoff[DETECTOR][i])

    plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
plt.xlim(minPlot, maxPlot)
plt.ylim(minPlot, maxPlot)
plt.xlabel("NonLinearized Turnoff (ADU) Run 13144M", fontsize = 18)
plt.ylabel("Linearized Turnoff (ADU) Run 13144N", fontsize = 18)
xplot = np.linspace(minPlot, maxPlot,100)
plt.plot(xplot, xplot, ls = '--', color='blue')
slope1 = 1.10
plt.plot(xplot,slope1*xplot, ls = '--', color='red', label = '%.3f'%slope1)
slope2 = 0.90
plt.plot(xplot,slope2*xplot, ls = '--', color='red', label = '%.3f'%slope2)

#plt.text(1.05, 1.55, "%d CCDs"%numCCDs, fontsize = 18)
#plt.text(1.05, 1.525, "%d amps"%numAmps, fontsize = 18)
plt.legend(loc='lower right', fontsize = 18)
plt.subplot(1,2,plotCounter+1)
plt.title("Turnoff Ratio (Lin / NonLin)", fontsize = 18)
n, bins, patches = plt.hist(ratios, bins = 50, range=(0.80,1.20))
ymax = n.max() * 1.10
plt.xlim(0.80, 1.20)
plt.ylim(0, ymax)
plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
plt.text(1.03, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
plt.text(1.03, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
plt.text(1.03, 0.85*n.max(), "%d amps"%len(ratios), fontsize = 18)
plt.savefig("/repo/main/u/cslage/bps_13144N/plots/Lin_NonLin_Turnoff_13144N_05Jan22.pdf")

```

```python

```
