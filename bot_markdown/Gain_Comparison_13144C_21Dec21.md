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
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib",\
                                                    "LSSTCam/calib/u/cslage/13144", "u/cslage/bps_13144C"])
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
PTC_gains = {}
EO_PTC_gains = {}
PTC_noise = {}
EO_PTC_noise = {}
PTC_a00 = {}
EO_PTC_a00 = {}
PTC_turnoff = {}
EO_PTC_turnoff = {}

detectors = {}
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

        detectors[DETECTOR] = [VENDOR, RAFT, SENSOR]
        PTC_gains[DETECTOR] = []
        EO_PTC_gains[DETECTOR] = []
        PTC_noise[DETECTOR] = []
        EO_PTC_noise[DETECTOR] = []
        PTC_a00[DETECTOR] = []
        EO_PTC_a00[DETECTOR] = []
        PTC_turnoff[DETECTOR] = []
        EO_PTC_turnoff[DETECTOR] = []

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
            PTC_gains[DETECTOR].append(ptcDataset.gain[amp])
            EO_PTC_gains[DETECTOR].append(eoPTCGain[slacNum])
            PTC_noise[DETECTOR].append(ptcDataset.noise[amp])
            EO_PTC_noise[DETECTOR].append(eoNoise[slacNum])
            PTC_a00[DETECTOR].append(dmA00)
            EO_PTC_a00[DETECTOR].append(-eoA00[slacNum])
            PTC_turnoff[DETECTOR].append(maxDM)
            EO_PTC_turnoff[DETECTOR].append(eoPtcTurnoff[slacNum])

```

```python
minGain = 1.2
maxGain = 2.0
plotCounter = 0
numCCDs = len(detectors)
numAmps = 0
num = 343
plt.figure(figsize=(16,8))
badCounts = np.zeros([189])
goodCount = 0
badList = open('/project/cslage/BOT_LSSTCam/eotest/badList_13144C.txt', 'w')
zeroList = open('/project/cslage/BOT_LSSTCam/eotest/zeroList_13144C.txt', 'w')
goodList = open('/project/cslage/BOT_LSSTCam/eotest/goodList_13144C.txt', 'w')
for [EO, name] in [[EO_PTC_gains, 'PTC-Run 13144C']]:
    ratios = []
    plotCounter += 1
    plt.subplot(1,2,plotCounter)
    plt.title("Gain Comparison", fontsize = 18)
    badx = []
    bady = []
    for [vendor, color] in [['E2V', 'red'], ['ITL', 'green']]:
        xplot = []
        yplot = []
        for DETECTOR in detectors.keys():
            [VENDOR, RAFT, SENSOR] = detectors[DETECTOR]
            for i in range(16):
                if EO[DETECTOR][i] >.001:
                    ratio = PTC_gains[DETECTOR][i] / EO[DETECTOR][i]
                else:
                    ratio = 1.0
                ratios.append(ratio)
                if vendor == VENDOR:
                    numAmps += 1
                    xplot.append(EO[DETECTOR][i])
                    yplot.append(PTC_gains[DETECTOR][i])
                    if name == 'PTC-Run 13144C' and (ratio < 0.001):
                        zeroList.write("%s\t%s\t%s\n"%(RAFT,SENSOR,list(ptcDataset.gain.keys())[i]))
                        print(ratio, VENDOR, RAFT, SENSOR, DETECTOR, list(ptcDataset.gain.keys())[i],PTC_gains[DETECTOR][i],EO[DETECTOR][i])
                    elif name == 'PTC-Run 13144C' and (ratio < 0.95 or ratio > 1.05):
                        # 20% - 5 bad, 10% - 8 bad, 5% - 16 bad
                        badCounts[int(DETECTOR)] += 1
                        badList.write("%s\t%s\t%s\n"%(RAFT,SENSOR,list(ptcDataset.gain.keys())[i]))
                        print(ratio, VENDOR, RAFT, SENSOR, DETECTOR, list(ptcDataset.gain.keys())[i],PTC_gains[DETECTOR][i],EO[DETECTOR][i])
                    elif name == 'PTC-Run 13144C' and (ratio > 0.95 and ratio < 1.05):
                        if np.random.rand() < .01:
                            goodCount += 1
                            goodList.write("%s\t%s\t%s\n"%(RAFT,SENSOR,list(ptcDataset.gain.keys())[i]))
                        #print(ratio, VENDOR, RAFT, SENSOR, DETECTOR, list(ptcGain.keys())[i])
        plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
    plt.xlim(minGain, maxGain)
    plt.ylim(minGain, maxGain)
    plt.xlabel("EOtest Gains(e-/ADU) %s "%name, fontsize = 18)
    plt.ylabel("DM PTC Gains(e-/ADU) (Run 13144); %s Flat Pairs"%num, fontsize = 18)
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
    plt.title("Gain Ratio (DM PTC / EO %s)"%name, fontsize = 18)
    n, bins, patches = plt.hist(ratios, bins = 50, range=(0.90,1.10))
    ymax = n.max() * 1.10
    plt.xlim(0.90, 1.10)
    plt.ylim(0, ymax)
    plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
    plt.text(1.01, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
    plt.text(1.01, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
    plt.text(1.01, 0.85*n.max(), "%d amps"%numAmps, fontsize = 18)
plt.savefig("/repo/main/u/cslage/bps_13144C/plots/PTC_Eotest_Gains_13144C_21Dec21.pdf")
badList.close()
goodList.close()
print(goodCount)
```

```python
minPlot = 1.0
maxPlot = 5.0
plotCounter = 0
numCCDs = len(detectors)
numAmps = 0
num = 343
plt.figure(figsize=(16,8))
for [EO, name] in [[EO_PTC_a00, 'PTC-Run 13144']]:
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
                #print(i, DETECTOR, EO[DETECTOR][i], PTC_a00[DETECTOR][i])
                if vendor == VENDOR:
                    numAmps += 1
                    eo_a00 = -EO[DETECTOR][i] * 1.0E6
                    if eo_a00 > minPlot and eo_a00 < maxPlot:

                        ratio = PTC_a00[DETECTOR][i] / EO[DETECTOR][i]
                        #print(eo_a00, PTC_a00[DETECTOR][i], ratio)
                        ratios.append(ratio)
                    #else:
                    #    ratio = 1.0
                    #ratios.append(ratio)

                    xplot.append(-EO[DETECTOR][i] * 1.0E6)
                    yplot.append(-PTC_a00[DETECTOR][i] * 1.0E6)
                    """
                    if name == 'PTC-Run 13038' and (ratio < 0.95 or ratio > 1.05):
                        # 20% - 5 bad, 10% - 8 bad, 5% - 16 bad
                        badCounts[int(DETECTOR)] += 1
                        badList.write("%s\t%s\t%s\n"%(RAFT,SENSOR,list(ptcGain.keys())[i]))
                        print(ratio, VENDOR, RAFT, SENSOR, DETECTOR, list(ptcGain.keys())[i],PTC_gains[DETECTOR][i],EO[DETECTOR][i])
                    if name == 'PTC-Run 13038' and (ratio > 0.95 and ratio < 1.05):
                        if np.random.rand() < .01:
                            goodCount += 1
                            goodList.write("%s\t%s\t%s\n"%(RAFT,SENSOR,list(ptcGain.keys())[i]))
                        #print(ratio, VENDOR, RAFT, SENSOR, DETECTOR, list(ptcGain.keys())[i])
                    """
        plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
    plt.xlim(minPlot, maxPlot)
    plt.ylim(minPlot, maxPlot)
    plt.xlabel("EOtest A00 %s "%name, fontsize = 18)
    plt.ylabel("DM PTC A00 (Run 13144); %s Flat Pairs"%num, fontsize = 18)
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
    plt.title("A00 Ratio (DM PTC / EO %s)"%name, fontsize = 18)
    n, bins, patches = plt.hist(ratios, bins = 50, range=(0.80,1.20))
    ymax = n.max() * 1.10
    plt.xlim(0.80, 1.20)
    plt.ylim(0, ymax)
    plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
    plt.text(1.03, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
    plt.text(1.03, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
    plt.text(1.03, 0.85*n.max(), "%d amps"%len(ratios), fontsize = 18)
plt.savefig("/repo/main/u/cslage/bps_13144C/plots/PTC_Eotest_A00_13144C_21Dec21.pdf")
#badList.close()
#goodList.close()
#print(goodCount)
```

```python
minPlot = 1.0
maxPlot = 20.0
plotCounter = 0
numCCDs = len(detectors)
numAmps = 0
num =343
plt.figure(figsize=(16,8))
for [EO, name] in [[EO_PTC_noise, 'PTC-Run 13144']]:
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
                #print(i, DETECTOR, EO[DETECTOR][i], PTC_a00[DETECTOR][i])
                if vendor == VENDOR:
                    numAmps += 1
                    eo_noise = EO[DETECTOR][i]
                    if eo_noise > minPlot and eo_noise < maxPlot:

                        ratio = PTC_noise[DETECTOR][i] / EO[DETECTOR][i]
                        #print(eo_a00, PTC_a00[DETECTOR][i], ratio)
                        ratios.append(ratio)
                    #else:
                    #    ratio = 1.0
                    #ratios.append(ratio)

                    xplot.append(EO[DETECTOR][i])
                    yplot.append(PTC_noise[DETECTOR][i])
                    """
                    if name == 'PTC-Run 13038' and (ratio < 0.95 or ratio > 1.05):
                        # 20% - 5 bad, 10% - 8 bad, 5% - 16 bad
                        badCounts[int(DETECTOR)] += 1
                        badList.write("%s\t%s\t%s\n"%(RAFT,SENSOR,list(ptcGain.keys())[i]))
                        print(ratio, VENDOR, RAFT, SENSOR, DETECTOR, list(ptcGain.keys())[i],PTC_gains[DETECTOR][i],EO[DETECTOR][i])
                    if name == 'PTC-Run 13038' and (ratio > 0.95 and ratio < 1.05):
                        if np.random.rand() < .01:
                            goodCount += 1
                            goodList.write("%s\t%s\t%s\n"%(RAFT,SENSOR,list(ptcGain.keys())[i]))
                        #print(ratio, VENDOR, RAFT, SENSOR, DETECTOR, list(ptcGain.keys())[i])
                    """
        plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
    plt.xlim(minPlot, maxPlot)
    plt.ylim(minPlot, maxPlot)
    plt.xlabel("EOtest Noise (e-) %s "%name, fontsize = 18)
    plt.ylabel("DM PTC Noise (e-) (Run 13144); %s Flat Pairs"%num, fontsize = 18)
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
    plt.title("Noise Ratio (DM PTC / EO %s)"%name, fontsize = 18)
    n, bins, patches = plt.hist(ratios, bins = 50, range=(0.80,1.20))
    ymax = n.max() * 1.10
    plt.xlim(0.80, 1.20)
    plt.ylim(0, ymax)
    plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
    plt.text(1.03, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
    plt.text(1.03, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
    plt.text(1.03, 0.85*n.max(), "%d amps"%len(ratios), fontsize = 18)
plt.savefig("/repo/main/u/cslage/bps_13144C/plots/PTC_Eotest_Noise_13144C_21Dec21.pdf")
#badList.close()
#goodList.close()
#print(goodCount)
```

```python
minPlot = 10000.0
maxPlot = 200000.0
plotCounter = 0
numCCDs = len(detectors)
numAmps = 0
num =343
badCount = 0
badList = open('/project/cslage/BOT_LSSTCam/eotest/turnoff_clump_itl_13144C.txt', 'w')
plt.figure(figsize=(16,8))
for [EO, name] in [[EO_PTC_turnoff, 'PTC-Run 13144']]:
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
                #print(i, DETECTOR, EO[DETECTOR][i], PTC_a00[DETECTOR][i])
                if vendor == VENDOR:
                    numAmps += 1
                    eo_turnoff = EO[DETECTOR][i]
                    if eo_turnoff > minPlot and eo_turnoff < maxPlot:

                        ratio = PTC_turnoff[DETECTOR][i] / EO[DETECTOR][i]
                        #print(eo_a00, PTC_a00[DETECTOR][i], ratio)
                        ratios.append(ratio)
                    #else:
                    #    ratio = 1.0
                    #ratios.append(ratio)

                    xplot.append(EO[DETECTOR][i])
                    yplot.append(PTC_turnoff[DETECTOR][i])
                    
                    if name == 'PTC-Run 13144' and eo_turnoff > 75000.0 and \
                    ratio < 0.85 and VENDOR == 'ITL':
                        badList.write("%s\t%s\t%s\n"%(RAFT,SENSOR,list(ptcDataset.gain.keys())[i]))
                        badCount += 1
                        print(ratio, VENDOR, RAFT, SENSOR, DETECTOR ,PTC_turnoff[DETECTOR][i],eo_turnoff)
                    """
                    if name == 'PTC-Run 13038' and (ratio > 0.95 and ratio < 1.05):
                        if np.random.rand() < .01:
                            goodCount += 1
                            goodList.write("%s\t%s\t%s\n"%(RAFT,SENSOR,list(ptcGain.keys())[i]))
                        #print(ratio, VENDOR, RAFT, SENSOR, DETECTOR, list(ptcGain.keys())[i])
                    """
        plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
    plt.xlim(minPlot, maxPlot)
    plt.ylim(minPlot, maxPlot)
    plt.xlabel("EOtest Turnoff (ADU) %s "%name, fontsize = 18)
    plt.ylabel("DM PTC Turnoff (ADU) (Run 13144); %s Flat Pairs"%num, fontsize = 18)
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
    plt.title("Turnoff Ratio (DM PTC / EO %s)"%name, fontsize = 18)
    n, bins, patches = plt.hist(ratios, bins = 50, range=(0.80,1.20))
    ymax = n.max() * 1.10
    plt.xlim(0.80, 1.20)
    plt.ylim(0, ymax)
    plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
    plt.text(1.03, n.max(), "Median = %.4f"%np.nanmedian(ratios), fontsize = 18)
    plt.text(1.03, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
    plt.text(1.03, 0.85*n.max(), "%d amps"%len(ratios), fontsize = 18)
plt.savefig("/repo/main/u/cslage/bps_13144C/plots/PTC_Eotest_Turnoff_13144C_21Dec21.pdf")
badList.close()
#goodList.close()
print(badCount)
```

```python

```
