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

# Notebook for comparing the impact of FULLCOVARIANCE on gain and other parameters.

Initially written 22 Jul 2022 by Craig Lage.

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
runA = "u/cslage/bps_13144R"
runB = "u/cslage/bps_13144S"
A_Butler = Butler("/repo/main", collections=[runA])
B_Butler = Butler("/repo/main", collections=[runB])
nameA = 'EXPAPPROX'
nameB = 'FULLCOV'
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
A_gains = {}
B_gains = {}
A_noise = {}
B_noise = {}
A_a00 = {}
B_a00 = {}
A_turnoff = {}
B_turnoff = {}

detectors = {}
for RAFT in rafts:
    for SENSOR in sensors:
        VENDOR, DETECTOR, raftRow, raftCol = detector(RAFT, SENSOR)
        try:
            A_PTC = A_Butler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')
            B_PTC = B_Butler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')            
        except:
            continue

        detectors[DETECTOR] = [VENDOR, RAFT, SENSOR]
        A_gains[DETECTOR] = []
        B_gains[DETECTOR] = []
        A_noise[DETECTOR] = []
        B_noise[DETECTOR] = []
        A_a00[DETECTOR] = []
        B_a00[DETECTOR] = []
        A_turnoff[DETECTOR] = []
        B_turnoff[DETECTOR] = []

        for amp in A_PTC.gain.keys():
            if A_PTC.ptcFitType == 'EXPAPPROXIMATION':
                A_A00 = A_PTC.ptcFitPars[amp][0]
            if A_PTC.ptcFitType == 'FULLCOVARIANCE':
                A_A00 = A_PTC.aMatrix[amp][0][0]
            if B_PTC.ptcFitType == 'EXPAPPROXIMATION':
                B_A00 = B_PTC.ptcFitPars[amp][0]
            if B_PTC.ptcFitType == 'FULLCOVARIANCE':
                B_A00 = B_PTC.aMatrix[amp][0][0]
                
            A_Means = np.array(A_PTC.finalMeans[amp])
            A_Means = A_Means[~np.isnan(A_Means)]
            if len(A_Means > 0):
                A_maxDM = A_Means.max()
            else:
                A_maxDM = 0.0
                
            B_Means = np.array(B_PTC.finalMeans[amp])
            B_Means = B_Means[~np.isnan(B_Means)]
            if len(B_Means > 0):
                B_maxDM = B_Means.max()
            else:
                B_maxDM = 0.0
            A_gains[DETECTOR].append(A_PTC.gain[amp])
            B_gains[DETECTOR].append(B_PTC.gain[amp])
            A_noise[DETECTOR].append(A_PTC.noise[amp])
            B_noise[DETECTOR].append(B_PTC.noise[amp])
            A_a00[DETECTOR].append(A_A00)
            B_a00[DETECTOR].append(B_A00)
            A_turnoff[DETECTOR].append(A_maxDM)
            B_turnoff[DETECTOR].append(B_maxDM)

```

```python
print(len(A_noise), len(B_noise))
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
            if B_gains[DETECTOR][i] >.001:
                ratio = A_gains[DETECTOR][i] / B_gains[DETECTOR][i]
            else:
                ratio = 1.0
            ratios.append(ratio)
            if vendor == VENDOR:
                numAmps += 1
                xplot.append(B_gains[DETECTOR][i])
                yplot.append(A_gains[DETECTOR][i])
    plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
plt.xlim(minGain, maxGain)
plt.ylim(minGain, maxGain)
plt.xlabel("Gains(e-/ADU) (Run %s)"%nameB, fontsize = 18)
plt.ylabel("Gains(e-/ADU) (Run %s)"%nameA, fontsize = 18)
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
plt.title("Gain Ratio", fontsize = 18)
n, bins, patches = plt.hist(ratios, bins = 50, range=(0.90,1.10))
ymax = n.max() * 1.10
plt.xlim(0.90, 1.10)
plt.ylim(0, ymax)
plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
plt.text(1.01, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
plt.text(1.01, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
plt.text(1.01, 0.85*n.max(), "%d amps"%numAmps, fontsize = 18)
plt.savefig("/repo/main/%s/plots/Gains_Comparison_13144_FullCov_22Jul22.pdf"%runB)

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
                b_a00 = -B_a00[DETECTOR][i] * 1.0E6
                if b_a00 > minPlot and b_a00 < maxPlot:

                    ratio = A_a00[DETECTOR][i] / B_a00[DETECTOR][i]
                    ratios.append(ratio)
                xplot.append(-B_a00[DETECTOR][i] * 1.0E6)
                yplot.append(-A_a00[DETECTOR][i] * 1.0E6)
    plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
plt.xlim(minPlot, maxPlot)
plt.ylim(minPlot, maxPlot)
plt.xlabel("A00 (Run %s)"%nameB, fontsize = 18)
plt.ylabel("(Run %s)"%nameA, fontsize = 18)
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
plt.title("A00 Ratio", fontsize = 18)
n, bins, patches = plt.hist(ratios, bins = 50, range=(0.80,1.20))
ymax = n.max() * 1.10
plt.xlim(0.80, 1.20)
plt.ylim(0, ymax)
plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
plt.text(1.03, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
plt.text(1.03, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
plt.text(1.03, 0.85*n.max(), "%d amps"%len(ratios), fontsize = 18)
plt.savefig("/repo/main/%s/plots/A00_Comparison__13144_FullCov_22Jul22.pdf"%runB)
```

```python
minPlot = 1.0
maxPlot = 200.0
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
                #print(DETECTOR, i, A_noise[DETECTOR][i], B_noise[DETECTOR][i])
                if B_noise[DETECTOR][i] > minPlot and B_noise[DETECTOR][i] < maxPlot:
                    ratio = A_noise[DETECTOR][i] / B_noise[DETECTOR][i]
                    #print(DETECTOR, i, ratio)
                    ratios.append(ratio)

                xplot.append(B_noise[DETECTOR][i])
                yplot.append(A_noise[DETECTOR][i])
    plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
plt.xlim(minPlot, maxPlot)
plt.ylim(minPlot, maxPlot)
plt.xlabel("Noise (e-) Run %s"%nameB, fontsize = 18)
plt.ylabel("Noise (e-) Run %s"%nameA, fontsize = 18)
xplot = np.linspace(minPlot, maxPlot,100)
plt.plot(xplot, xplot, ls = '--', color='blue')
slope1 = 1.10
plt.plot(xplot,slope1*xplot, ls = '--', color='red', label = '%.3f'%slope1)
slope2 = 0.90
plt.plot(xplot,slope2*xplot, ls = '--', color='red', label = '%.3f'%slope2)

#plt.text(1.05, 1.55, "%d CCDs"%numCCDs, fontsize = 18)
#plt.text(1.05, 1.525, "%d amps"%numAmps, fontsize = 18)
plt.legend(loc='upper left', fontsize = 18)
plt.subplot(1,2,plotCounter+1)
plt.title("Noise Ratio ", fontsize = 18)
n, bins, patches = plt.hist(ratios, bins = 50, range=(0.0,0.5))
ymax = n.max() * 1.10
plt.xlim(0.0, 0.5)
plt.ylim(0, ymax)
plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
plt.text(0.3, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
plt.text(0.3, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
plt.text(0.3, 0.85*n.max(), "%d amps"%len(ratios), fontsize = 18)
plt.savefig("/repo/main/%s/plots/Noise_Comparison__13144_FullCov_22Jul22.pdf"%runB)
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
                turnoff = B_turnoff[DETECTOR][i]
                if turnoff > minPlot and turnoff < maxPlot:
                    ratio = A_turnoff[DETECTOR][i] / B_turnoff[DETECTOR][i]
                    ratios.append(ratio)

                xplot.append(B_turnoff[DETECTOR][i])
                yplot.append(A_turnoff[DETECTOR][i])

    plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
plt.xlim(minPlot, maxPlot)
plt.ylim(minPlot, maxPlot)
plt.xlabel("Turnoff (ADU) Run %s"%nameB, fontsize = 18)
plt.ylabel("Turnoff (ADU) Run %s"%nameA, fontsize = 18)
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
plt.savefig("/repo/main/%s/plots/Turnoff_Comparison__13144_FullCov_22Jul22.pdf"%runB)
```

```python
from scipy.stats import sigmaclip

minPlot = 1.0
maxPlot = 4.0
plotCounter = 0
numCCDs = len(detectors)
numAmps = 0
names = ['EXPAPPROX', 'FULLCOV']
plt.figure(figsize=(16,16))

plt.suptitle("A00 Histograms", fontsize = 18)
for vendor in ['E2V', 'ITL']:
    for i, data in enumerate([A_a00, B_a00]):
        plotCounter += 1
        plt.subplot(2,2,plotCounter)
        plt.title(f"{vendor} - {names[i]}")
        plotData = []
        for DETECTOR in detectors.keys():
            [VENDOR, RAFT, SENSOR] = detectors[DETECTOR]
            for ii in range(16):
                if vendor == VENDOR:
                    plotData.append(-data[DETECTOR][ii] * 1.0E6)
        plotData = np.array(plotData)
        totalNum = len(plotData)
        plotData = plotData[~np.isnan(plotData)]
        goodNum = len(plotData)
        plotData = sigmaclip(plotData, low=4.0, high=4.0)[0]
        unclippedNum = len(plotData)
        print(f"{vendor}, {names[i]}, {totalNum-goodNum} NaN values, {goodNum-unclippedNum} clipped values.")
        n, bins, patches = plt.hist(plotData, bins = 50, range=(minPlot,maxPlot))
        ymax = n.max() * 1.10
        plt.xlim(minPlot, maxPlot)
        plt.ylim(0, ymax)
        if vendor == 'E2V':
            plotx = 1.2
        else:
            plotx = 2.3
        plt.text(plotx, n.max(), "Clipped mean = %.4f"%np.nanmean(plotData), fontsize = 18)
        plt.text(plotx, 0.9*n.max(), "Clipped std = %.4f"%np.nanstd(plotData), fontsize = 18)
        plt.text(plotx, 0.8*n.max(), "%d Amps"%len(plotData), fontsize = 18)
plt.savefig("/repo/main/%s/plots/A00_Histograms_13144_25Jul22.pdf"%runB)
```

```python
goodAmps = 0
for RAFT in rafts:
    for SENSOR in sensors:
        VENDOR, DETECTOR, raftRow, raftCol = detector(RAFT, SENSOR)
        try:
            B_PTC = B_Butler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')            
        except:
            continue
        for amp in B_PTC.gain.keys():
            ampGood = True
            numGood = len(B_PTC.finalMeans[amp])
            for i in range(numGood):
                if not np.isnan(B_PTC.finalMeans[amp][i]):
                    if np.isnan(B_PTC.covariancesModel[amp][i]).any():
                        print(DETECTOR, amp, i, "has NaNs")
                        ampGood = False
            if ampGood:
                goodAmps += 1
print(goodAmps, "good amps.")
```

```python

```
