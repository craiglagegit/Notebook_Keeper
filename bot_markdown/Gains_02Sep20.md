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

# Notebook for plotting BOT gains.

Initially written 28 Aug 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
! eups list -s | grep cp_pipe
```

```python
import sys, os, glob, time
import numpy as np
import matplotlib.pyplot as plt
from lsst.daf.persistence import Butler
```

```python
REPO_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_12606'
butler = Butler(REPO_DIR)
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
    raftCol = int(list(raft)[2]) - startingCol[raftRow]
    sensorRow = int(list(sensor)[1])
    sensorCol = int(list(sensor)[2])
    detectorNum = (rows[raftRow] + raftCol) * 9
    detectorNum += 3 * sensorRow + sensorCol
    plotNum = 21 - 5 * raftRow + int(list(raft)[2])
    return vendor, detectorNum, plotNum

```

```python
min = 0.5
max = 2.5


plt.figure(figsize=(16,16))
plt.subplots_adjust(hspace=0.3, wspace = 0.3)
plt.suptitle("BOT Gains - Run 12606 - 2020-10-09", fontsize=24)

badDetectors = [0, 166]                                                                                                                                                                        
rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']



e2vs = []
itls = []
missingE2vs = []
missingItls = []
for RAFT in rafts:
    gains = []
    gain_err = []
    xaxis = []

    for xValue, SENSOR in enumerate(sensors):
        VENDOR, DETECTOR, plotNum = detector(RAFT,SENSOR)
        plt.subplot(5,5,plotNum)
        plt.title("%s"%RAFT, fontsize = 12)
        plt.ylabel("Gain", fontsize = 12)

        if VENDOR == 'E2V':
            color = 'green'
        elif VENDOR == 'ITL':
            color = 'blue'
        if DETECTOR in badDetectors:
            plt.scatter([xValue,xValue], [1.5,1.5], marker = 'o', color = 'red')
            continue
        try:
            ptcDataset = butler.get('photonTransferCurveDataset', detector=DETECTOR)
        except:
            plt.scatter([xValue,xValue], [1.5,1.5], marker = 'o', color = 'red')
            continue
        gain_data = ptcDataset.gain
        gain_err_data = ptcDataset.gainErr
        amps = gain_data.keys()
        for amp in amps:
            gain = gain_data[amp]
            if RAFT in rafts and (gain < min or gain > max or np.isnan(gain)):
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
#plt.savefig(REPO_DIR+'/plots/Gain_Summary_03Sep20.pdf')

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
plt.suptitle("BOT Gains - Run 12539 - 2020-09-01", fontsize=24)
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
#plt.savefig(REPO_DIR+'/plots/Gain_Histograms_03Sep20.pdf')
```

```python

```
