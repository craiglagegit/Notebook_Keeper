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

# Notebook for investigating linearity corrections

Initially written 20 Dec 2021 by Craig Lage\
copying from Chris Waters.

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.io.fits as pf
from lsst.daf.butler import Butler
import lsst.afw.math as afwMath
from focal_plane_plotting import plot_focal_plane
```

```python
!eups list -s | grep ip_isr
!eups list -s | grep cp_pipe
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib"])
camera = butler.get('camera', instrument='LSSTCam')
```

```python
ptcButler = Butler("/repo/main", collections=["u/cslage/bps_13144M"])
```

```python
linButler = Butler("/repo/main", collections=["u/cslage/linearizerB_26jan22"])
```

```python
numDetectors = 0
for item in linButler.registry.queryDatasets('unCorrectedLinearizer'):
    numDetectors += 1
print(numDetectors)    
```

```python
# With the whole focal plane contributing to the correction, the residuals are much improved
det = 74
ptc = ptcButler.get('ptc', detector=det, instrument='LSSTCam')
uncorrLin = linButler.get('unCorrectedLinearizer', detector=det, instrument='LSSTCam')
corrLin = linButler.get('linearizer', detector=det, instrument='LSSTCam')

fig = plt.figure(figsize = (8,4))

for amp in camera[0].getAmplifiers():
    ampName = amp.getName()
    mask = np.array(ptc.expIdMask[ampName], dtype=bool)
    means = np.array(ptc.rawMeans[ampName])[mask]
    uncorrResiduals = np.array(uncorrLin.fitResiduals[ampName])[mask]
    corrResiduals = np.array(corrLin.fitResiduals[ampName])[mask]
    plt.title(f"Residuals - Det {det} - {ampName}", fontsize = 18)
    plt.scatter(means, uncorrResiduals/means * 100.0, label = "Uncorrected")
    plt.scatter(means, corrResiduals/means * 100.0, label = "Corrected")
    plt.plot([0.0,100000.0], [0.0,0.0], ls = '--', color='black')
    plt.ylim(-0.1,0.1)
    plt.legend()
    plt.xlabel("Flux(ADU)", fontsize = 12)
    plt.ylabel("Linearizer residual (%%)", fontsize = 12)
    break
plt.savefig("/repo/main/u/cslage/linearizerB_26jan22/plots/Residuals_26Jan22.pdf")
```

So this all seems to be working as intended.  On to the next step.

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
    return vendor, detectorNum
```

```python
names = []
corrStds = []
uncorrStds = []
fluxMin = 10000.0

for RAFT in ['R01',  'R02',  'R03', 'R10',  'R11',  'R12',  'R13', 'R14', 'R20',  'R21',  'R22',  'R23', 'R24', \
             'R30', 'R31', 'R32', 'R33', 'R34', 'R41', 'R42', 'R43']:
    for SENSOR in ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']:
        VENDOR, DETECTOR = detector(RAFT,SENSOR)
        try:
            ptc = ptcButler.get('ptc', detector=DETECTOR, instrument='LSSTCam')
            uncorrLin = linButler.get('unCorrectedLinearizer', detector=DETECTOR, instrument='LSSTCam')
            corrLin = linButler.get('linearizer', detector=DETECTOR, instrument='LSSTCam')
        except:
            continue

        for amp in camera[0].getAmplifiers():
            ampName = amp.getName()
            mask = np.array(ptc.expIdMask[ampName], dtype=bool)
            means = np.array(ptc.rawMeans[ampName])[mask]
            uncorrResiduals = np.array(uncorrLin.fitResiduals[ampName])[mask]
            corrResiduals = np.array(corrLin.fitResiduals[ampName])[mask]
            fluxMask = means > fluxMin
            corrStd = np.nanstd((corrResiduals/means * 100.0)[fluxMask])
            uncorrStd = np.nanstd((uncorrResiduals/means * 100.0)[fluxMask])
            names.append(f"{RAFT}_{SENSOR}_{ampName}")
            corrStds.append(corrStd)
            uncorrStds.append(uncorrStd)

```

```python
# Now plot it
xaxis = list(range(len(names)))
fig = plt.figure(figsize=(20,4))
plt.title(f"Standard Deviation of Residuals Run 13144, Flux > {fluxMin} ADU", fontsize = 24)
plt.scatter(xaxis, uncorrStds, marker = ".", s = 10, label = "Uncorrected")
plt.scatter(xaxis, corrStds, marker = ".", s = 10, label = "Corrected")
plt.ylim(0,0.05)
plt.legend(fontsize=12)
plt.xlabel("Amplifier Index", fontsize = 18)
plt.ylabel("Residual Standard Deviation (%)", fontsize = 18)

plt.savefig("/repo/main/u/cslage/linearizerB_26jan22/plots/Residual_Std_26Jan22.pdf")
```

```python
len(names)
```

```python

```
