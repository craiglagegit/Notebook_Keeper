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

Initially written 21 May 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
! eups list -s cp_pipe
```

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
```

```python
num = 'All'
REPO_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_12597'
DATA_DIR = '/project/shared/BOT/'
run = '12597'
```

```python
# Get the eotest results
filename = "/project/cslage/BOT_LSSTCam/eotest/eotest_gain_07oct20.pkl"
file = open(filename, 'rb')
fe55_results = pkl.load(file)
ptc_results = pkl.load(file)
file.close()

rafts = ['R34', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']
print(rafts)
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
    raftCol = int(list(raft)[2]) - startingCol[raftRow]
    sensorRow = int(list(sensor)[1])
    sensorCol = int(list(sensor)[2])
    detectorNum = (rows[raftRow] + raftCol) * 9
    detectorNum += 3 * sensorRow + sensorCol
    plotNum = 21 - 5 * raftRow + int(list(raft)[2])
    return vendor, detectorNum, plotNum

```

```python
butler = Butler(REPO_DIR)

BF_gains = {}
PTC_gains = {}
EO_PTC_gains = {}
EO_Fe55_gains = {}
detectors = {}
for RAFT in rafts:
    for SENSOR in sensors:
        VENDOR, DETECTOR, plotNum = getDetector(RAFT, SENSOR)
        try:
            ptc_data = butler.get('photonTransferCurveDataset', dataId={'detector':DETECTOR, 'run':run})
            print("Found detector %d"%DETECTOR)
        except:
            print("Didn't find detector %d"%DETECTOR)
            continue
        #gain = bf_kernel.gain
        #gainErr = bf_kernel.gainErr
        ptcGain = ptc_data.gain
        ptcGainErr = ptc_data.gainErr
        eotestPTCGain = ptc_results['ptc_gain'][RAFT][SENSOR]
        eotestPTCGainErr = ptc_results['ptc_gain_error'][RAFT][SENSOR]
        eotestFe55Gain = fe55_results['gain'][RAFT][SENSOR]
        eotestFe55GainErr = fe55_results['gain_error'][RAFT][SENSOR]
        detectors[DETECTOR] = [VENDOR, RAFT, SENSOR]
        PTC_gains[DETECTOR] = []
        EO_PTC_gains[DETECTOR] = []
        EO_Fe55_gains[DETECTOR] = []

        for i, amp in enumerate(ptcGain.keys()):
            PTC_gains[DETECTOR].append(ptcGain[amp])
            EO_PTC_gains[DETECTOR].append(eotestPTCGain[i])
            EO_Fe55_gains[DETECTOR].append(eotestFe55Gain[i])

print(ptcGain.keys())
#print(eotestPTCGain.keys())
```

```python
minGain = 1.0
maxGain = 1.6
plotCounter = 0
numCCDs = len(detectors)
plt.figure(figsize=(16,16))
for [EO, name] in [[EO_PTC_gains, 'PTC-Run 12597'], [EO_Fe55_gains, 'Fe55-Run 12534']]:
    ratios = []
    plotCounter += 1
    plt.subplot(2,2,plotCounter)
    plt.title("Gain Comparison", fontsize = 18)
    badx = []
    bady = []
    for [vendor, color] in [['E2V', 'red'], ['ITL', 'green']]:
        xplot = []
        yplot = []
        for DETECTOR in detectors.keys():
            #if DETECTOR !=98:
            #    continue
            [VENDOR, RAFT, SENSOR] = detectors[DETECTOR]
            for i in range(16):
                #if list(ptcGain.keys())[i] in ['C07', 'C10', 'C00', 'C10']:
                #    badx.append(EO[DETECTOR][i])
                #    bady.append(PTC_gains[DETECTOR][i])
                #    continue
                ratio = PTC_gains[DETECTOR][i] / EO[DETECTOR][i]
                ratios.append(ratio)
                if vendor == VENDOR:
                    xplot.append(EO[DETECTOR][i])
                    yplot.append(PTC_gains[DETECTOR][i])
                    if name == 'PTC-Run 12597' and (ratio < 0.975 or ratio > 1.025):
                        print(EO[DETECTOR][i], PTC_gains[DETECTOR][i], VENDOR, RAFT, SENSOR, DETECTOR, list(ptcGain.keys())[i])

        plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
    #plt.scatter(badx, bady, color='magenta', marker='*',label='Edge')
    plt.xlim(minGain, maxGain)
    plt.ylim(minGain, maxGain)
    plt.xlabel("EOtest Gains(e-/ADU) %s "%name, fontsize = 18)
    plt.ylabel("DM PTC Gains(e-/ADU) (Run 12597); %s Flat Pairs"%num, fontsize = 18)
    xplot = np.linspace(minGain, maxGain,100)
    plt.plot(xplot, xplot, ls = '--', color='blue')
    slope1 = 1.02
    plt.plot(xplot,slope1*xplot, ls = '--', color='red', label = '%.3f'%slope1)
    slope2 = 0.98
    plt.plot(xplot,slope2*xplot, ls = '--', color='red', label = '%.3f'%slope2)

    plt.text(1.05, 1.55, "%d CCDs"%numCCDs, fontsize = 18)
    plt.legend(loc='lower right', fontsize = 18)
    plt.subplot(2,2,plotCounter+2)
    plt.title("Gain Ratio (DM PTC / EO %s)"%name, fontsize = 18)
    n, bins, patches = plt.hist(ratios, bins = 50, range=(0.90,1.10))
    ymax = n.max() * 1.10
    plt.xlim(0.90, 1.10)
    plt.ylim(0, ymax)
    plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
    plt.text(1.01, n.max(), "Mean = %.4f"%np.nanmean(ratios), fontsize = 18)
    plt.text(1.01, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
plt.savefig(REPO_DIR+"/plots/PTC_Eotest_Gains_12597_07Oct20.pdf")
```

```python
EO_PTC_gains_2 = {}
EO_Fe55_gains_2 = {}
detectors_2 = {}

for RAFT in rafts:
    for SENSOR in ptc_results['ptc_gain'][RAFT].keys():
        VENDOR, DETECTOR, plotNum = getDetector(RAFT, SENSOR)
        eotestPTCGain = ptc_results['ptc_gain'][RAFT][SENSOR]
        eotestPTCGainErr = ptc_results['ptc_gain_error'][RAFT][SENSOR]
        eotestFe55Gain = fe55_results['gain'][RAFT][SENSOR]
        eotestFe55GainErr = fe55_results['gain_error'][RAFT][SENSOR]
        detectors_2[DETECTOR] = [VENDOR, RAFT, SENSOR]
        EO_PTC_gains_2[DETECTOR] = []
        EO_Fe55_gains_2[DETECTOR] = []
        for i in range(16):
            EO_PTC_gains_2[DETECTOR].append(eotestPTCGain[i])
            EO_Fe55_gains_2[DETECTOR].append(eotestFe55Gain[i])

numCCDs = len(detectors_2)
minGain = 1.0
maxGain = 1.6

plt.figure(figsize=(16,8))

ratios = []
plt.subplot(1,2,1)
plt.title("Gain Comparison", fontsize = 18)
for [vendor, color] in [['E2V', 'red'], ['ITL', 'green']]:
    xplot = []
    yplot = []
    for DETECTOR in detectors.keys():
        [VENDOR, RAFT, SENSOR] = detectors[DETECTOR]
        for i in range(16):
            ratio = EO_PTC_gains_2[DETECTOR][i] / EO_Fe55_gains_2[DETECTOR][i]
            ratios.append(ratio)
            if vendor == VENDOR:
                xplot += EO_Fe55_gains_2[DETECTOR]
                yplot += EO_PTC_gains_2[DETECTOR]
    plt.scatter(xplot, yplot, color=color, marker='*',label=vendor)
plt.xlim(minGain, maxGain)
plt.ylim(minGain, maxGain)
plt.xlabel("EOtest Fe55-Run 12534 Gains(e-/ADU)", fontsize = 18)
plt.ylabel("EOtest PTC-Run-12543 Gains(e-/ADU); All Flat Pairs", fontsize = 18)
xplot = np.linspace(minGain, maxGain,100)
plt.plot(xplot, xplot, ls = '--', color='blue')
#plt.plot(xplot, 1.06*xplot, ls = '--', color='black')
plt.legend(loc='lower right', fontsize = 18)
plt.subplot(1,2,2)
plt.title("Gain Ratio (EO PTC / EO Fe55)", fontsize = 18)
n, bins, patches = plt.hist(ratios, bins = 50, range=(0.90,1.10))
ymax = n.max() * 1.10
plt.xlim(0.90, 1.10)
plt.ylim(0, ymax)
plt.plot([1.0,1.0], [0.0,ymax], color = 'black', ls = '--')
plt.text(1.01, n.max(), "Mean = %.4f"%np.nanmean(ratios), fontsize = 18)
plt.text(1.01, 0.9*n.max(), "%d CCDs"%numCCDs, fontsize = 18)
plt.savefig(REPO_DIR+"/plots/Eotest_PTC_Eotest_Fe55_Gains_12543_03Sep20.pdf")
```

```python

```

```python
min = 0.5
max = 2.5


plt.figure(figsize=(16,16))
plt.subplots_adjust(hspace=0.3, wspace = 0.3)
plt.suptitle("BOT Gains - Run 12543 - 2020-09-02", fontsize=24)

badDetectors = []                                                                                                                                                                        
rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

theseRafts = ['R02', 'R03','R12', 'R13', 'R14', 'R21', 'R22', 'R23', 'R24',\
             'R32', 'R33', 'R34', 'R42', 'R43']

#rafts = ['R31']
e2vs = []
itls = []
missingE2vs = []
missingItls = []
for RAFT in rafts:
    gains = []
    gain_err = []
    xaxis = []

    for xValue, SENSOR in enumerate(sensors):
        VENDOR, DETECTOR, plotNum = getDetector(RAFT,SENSOR)
        plt.subplot(5,5,plotNum)
        plt.title("%s"%RAFT, fontsize = 12)
        plt.ylabel("Gain", fontsize = 12)

        if VENDOR == 'E2V':
            color = 'green'
        elif VENDOR == 'ITL':
            color = 'blue'
        if DETECTOR in badDetectors:
            continue
        try:
            ptcDataset = butler.get('photonTransferCurveDataset', detector=DETECTOR)
        except:
            continue
        gain_data = ptcDataset.gain
        gain_err_data = ptcDataset.gainErr
        amps = gain_data.keys()
        for amp in amps:
            gain = gain_data[amp]
            if RAFT in theseRafts and (gain < min or gain > max or np.isnan(gain)):
                missing = "%s%s_%s"%(RAFT,SENSOR,amp)
                if VENDOR == 'E2V':
                    missingE2vs.append(missing)
                elif VENDOR == 'ITL':
                    missingItls.append(missing)
            if np.isnan(gain):
                continue
            gains.append(gain)
            gain_err.append(gain_err_data[amp])
            xaxis.append(xValue)
            if VENDOR == 'E2V':
                e2vs.append(gain)
            elif VENDOR == 'ITL':
                itls.append(gain)
    plt.scatter(xaxis, gains, marker = 'x', color = color)
    plt.xticks(list(range(9)),sensors, fontsize=8)
    plt.ylim(0.0,2.0)
plt.savefig(REPO_DIR+'/plots/Gain_Summary_12543_03Sep20.pdf')

```

```python
print(len(e2vs), len(itls))
print(missingE2vs)
print(missingItls)
```

```python
hists = [['E2V', e2vs, 9*9*16, missingE2vs], ['ITL', itls, 4*9*16, missingItls]]


plt.figure(figsize=(16,8))
plt.subplots_adjust(hspace=0.5)
plt.suptitle("BOT Gains - Run 12543 - 2020-09-02", fontsize=24)
plotcounter = 0
for [name, data, numAmps, missing] in hists:
    data = np.array(data)
    plotcounter += 1
    plt.subplot(1,2,plotcounter)
    plt.title("%s"%(name), fontsize = 24)
    n, bins, patches = plt.hist(data, bins = 50, range=(min, max))
    ymax = n.max() * 1.10
    index = np.where((data>min) & (data<max))
    plt.text(0.7*max, n.max(), "Mean = %.3f"%np.nanmean(data[index]), fontsize = 18)
    plt.text(0.7*max, 0.95*n.max(), "Std = %.3f"%np.nanstd(data[index]), fontsize = 18)
    plt.text(0.7*max, 0.90*n.max(), "n = %d / %d"%(len(data[index]), numAmps), fontsize = 18)
    for i, miss in enumerate(missing):
        plt.text(0.7*max, (0.85-i*0.05)*n.max(), "Missing = %s"%miss, fontsize = 12)
    plt.xlim(min, max)
    plt.ylim(0, ymax)
    plt.xlabel("Gain (e-/ADU)", fontsize = 18)
plt.savefig(REPO_DIR+'/plots/Gain_Histograms_12543_03Sep20.pdf')
```

```python
print(list(ptcGain.keys())[7])
```

```python
1.1077557802200317 1.0807223562147568 E2V R22 S22 98 C10
1.1147732734680176 1.1105905713832454 E2V R22 S22 98 C11
1.1244839429855347 0.7182354957659461 E2V R22 S22 98 C12
1.1306294202804565 1.1243714384411552 E2V R22 S22 98 C13
1.1285284757614136 0.7225705274503295 E2V R22 S22 98 C14
1.1288763284683228 1.1199200829709957 E2V R22 S22 98 C15
1.1300760507583618 1.1431124521479339 E2V R22 S22 98 C16
1.1291757822036743 1.1247472603409925 E2V R22 S22 98 C17
1.1418300867080688 1.119140498700793 E2V R22 S22 98 C07
1.1493253707885742 1.1387970154027613 E2V R22 S22 98 C06
1.1544439792633057 1.148634931310022 E2V R22 S22 98 C05
1.1525565385818481 1.147266900605653 E2V R22 S22 98 C04
1.1454555988311768 1.1474838349757484 E2V R22 S22 98 C03
1.1598831415176392 1.1513788283318356 E2V R22 S22 98 C02
1.1605957746505737 1.1572257881715746 E2V R22 S22 98 C01
1.157067894935608 1.1656160439519492 E2V R22 S22 98 C00

1.1077557802200317 1.0849886382128175 E2V R22 S22 98 C10
1.1147732734680176 1.1105905713832454 E2V R22 S22 98 C11
1.1244839429855347 1.1099566671326877 E2V R22 S22 98 C12
1.1306294202804565 1.1243714384411552 E2V R22 S22 98 C13
1.1285284757614136 1.1115686417382036 E2V R22 S22 98 C14
1.1288763284683228 1.1199200829709957 E2V R22 S22 98 C15
1.1300760507583618 1.1431124521479339 E2V R22 S22 98 C16
1.1291757822036743 1.1247472603409925 E2V R22 S22 98 C17
1.1418300867080688 1.119140498700793 E2V R22 S22 98 C07
1.1493253707885742 1.1508901724165932 E2V R22 S22 98 C06
1.1544439792633057 1.148634931310022 E2V R22 S22 98 C05
1.1525565385818481 1.147266900605653 E2V R22 S22 98 C04
1.1454555988311768 1.1474838349757484 E2V R22 S22 98 C03
1.1598831415176392 1.1513788283318356 E2V R22 S22 98 C02
1.1605957746505737 1.1572257881715746 E2V R22 S22 98 C01
1.157067894935608 1.1656160439519492 E2V R22 S22 98 C00

```
