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

Initially written 27 Jan 2022 by Craig Lage

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
ptcButler = Butler("/repo/main", collections=["u/cslage/bps_13144R"])
```

```python
# Re-testing after edit to remove "exposure" from dimensions
linButler = Butler("/repo/main", collections=["u/cslage/linearizer_28jan22"])
```

```python
numDetectors = 0
for item in linButler.registry.queryDatasets('unCorrectedLinearizer'):
    numDetectors += 1
print(numDetectors)    
```

```python tags=[]
# First, get the data and put it in a dict as the heatmap plotting needs
# Identify thr sequence number that gives the appropriate fluxes
# Pick an amp near the center
detId = 94
ampName = 'C14'
fluxHi = 20000.0
fluxLo = 1000.0
ptc = ptcButler.get('ptc', detector=detId, instrument='LSSTCam')
mask = np.array(ptc.expIdMask[ampName], dtype=bool)
means = np.array(ptc.rawMeans[ampName])[mask]
gain = ptc.gain[ampName]
foundFluxLo = False
for i, mean in enumerate(means):
    if mean * gain > fluxLo and not foundFluxLo:
        expIdLo = ptc.inputExpIdPairs[ampName][i][0][0]
        seqNoLo = i
        print(f"ExpIdLo = {expIdLo}, seqNoLo = {seqNoLo}")
        foundFluxLo = True
    if mean * gain > fluxHi:
        seqNoHi = i
        expIdHi = ptc.inputExpIdPairs[ampName][i][0][0]
        print(f"ExpIdHi = {expIdHi}, seqNoHi = {seqNoHi}")
        break

fluxMin = 10000.0
corrStds = dict()
uncorrStds = dict()
maxNL = dict()
fluxesLo = dict()
fluxesHi = dict()
noises = dict()
gains = dict()

badAmps = []
weakAmps = []

for detector in camera:
    if detector.getType().name != 'SCIENCE':
        continue
    detName = detector.getName()
    detId = detector.getId()
    corrStds[detName] = dict()
    uncorrStds[detName] = dict()
    maxNL[detName] = dict()
    fluxesLo[detName] = dict()
    fluxesHi[detName] = dict()
    noises[detName] = dict()
    gains[detName] = dict()
    try:
        ptc = ptcButler.get('ptc', detector=detId, instrument='LSSTCam')
        uncorrLin = linButler.get('unCorrectedLinearizer', detector=detId, instrument='LSSTCam')
        corrLin = linButler.get('linearizer', detector=detId, instrument='LSSTCam')
    except:
        continue

    for amp in detector.getAmplifiers():
        ampName = amp.getName()
        mask = np.array(ptc.expIdMask[ampName], dtype=bool)
        means = np.array(ptc.rawMeans[ampName])[mask]
        uncorrResiduals = np.array(uncorrLin.fitResiduals[ampName])[mask]
        corrResiduals = np.array(corrLin.fitResiduals[ampName])[mask]
        fluxMask = means > fluxMin
        corrStd = np.nanstd((corrResiduals/means * 100.0)[fluxMask])
        uncorrStd = np.nanstd((uncorrResiduals/means * 100.0)[fluxMask])
        corrStds[detName][ampName] = corrStd
        uncorrStds[detName][ampName] = uncorrStd
        centers, values = np.split(corrLin.linearityCoeffs[ampName], 2)
        fluxMask = np.where(centers>fluxMin)
        try:
            maxDeviation = np.max(abs((values/centers * 100.0)[fluxMask]))
        except:
            maxDeviation = np.nan
        maxNL[detName][ampName] = maxDeviation
        loFlux = ptc.rawMeans[ampName][seqNoLo] * ptc.gain[ampName]
        hiFlux = ptc.rawMeans[ampName][seqNoHi] * ptc.gain[ampName]
        if loFlux < fluxLo / 10.0 or np.isnan(loFlux) or hiFlux < fluxHi / 10.0 or np.isnan(hiFlux):
            badAmps.append(f"{detName}_{ampName}")
        elif loFlux < fluxLo * 0.6 or hiFlux < fluxHi * 0.6:
            weakAmps.append(f"{detName}_{ampName}")
        fluxesLo[detName][ampName] = loFlux
        fluxesHi[detName][ampName] = hiFlux
        noises[detName][ampName] = ptc.noise[ampName]
        gains[detName][ampName] = ptc.gain[ampName]

print("badAmps", badAmps)
print("weakAmps", weakAmps)
```

```python
# Check it
fluxesLo['R22_S11']['C12']
```

```python
# Now plot the corrStd heatmap. Great code from Jim, as usual.
fig = plt.figure(figsize=(11,8.5))
ax = plt.axes([0.1, 0.1, 0.8, 0.8])
title = "Standard deviation of corrected linearizer residuals.  Run 13144"
plot_focal_plane(ax, corrStds, camera=camera, z_range=(0.0, 0.025), title=title)
plt.savefig("/repo/main/u/cslage/linearizer_28jan22/plots/Std_HeatMap_28Jan22.pdf")
```

```python
# Now plot the MaxNL heatmap
fig = plt.figure(figsize=(11,8.5))
ax = plt.axes([0.1, 0.1, 0.8, 0.8])
title = "Max deviation from linearity (%) based on spline knots.  Run 13144"
plot_focal_plane(ax, maxNL, camera=camera, z_range=(0.0, 1.0), title=title)
plt.savefig("/repo/main/u/cslage/linearizer_28jan22/plots/MaxNL_HeatMap_28Jan22.pdf")
```

```python
# Now plot the flux heatmap low flux
fig = plt.figure(figsize=(11,8.5))
ax = plt.axes([0.1, 0.1, 0.8, 0.8])
title = f"Flux distribution in electrons.  Run 13144, Flux ~ {fluxLo}"
plot_focal_plane(ax, fluxesLo, camera=camera, z_range=(0.0, fluxLo * 1.4), title=title)
plt.text(-322, -250, "Bad amps")
plt.text(-322, -270, badAmps[0])
plt.text(-322, -290, badAmps[1])
plt.text(-322, -310, badAmps[2])
plt.text(235, -250, "Weak amps")
plt.text(235, -270, weakAmps[0])
plt.text(235, -290, weakAmps[1])
plt.text(235, -310, weakAmps[2])

plt.savefig("/repo/main/u/cslage/linearizer_28jan22/plots/Low_Flux_HeatMap_02Feb22.pdf")

# Need to add a list of weak/bad amps
```

```python
# Now plot the flux heatmap high flux
fig = plt.figure(figsize=(11,8.5))
ax = plt.axes([0.1, 0.1, 0.8, 0.8])
title = f"Flux distribution in electrons.  Run 13144, flux ~ {fluxHi}"
plot_focal_plane(ax, fluxesHi, camera=camera, z_range=(0.0, fluxHi * 1.4), title=title)
plt.text(-322, -250, "Bad amps")
plt.text(-322, -270, badAmps[0])
plt.text(-322, -290, badAmps[1])
plt.text(-322, -310, badAmps[2])
plt.text(235, -250, "Weak amps")
plt.text(235, -270, weakAmps[0])
plt.text(235, -290, weakAmps[1])
plt.text(235, -310, weakAmps[2])

plt.savefig("/repo/main/u/cslage/linearizer_28jan22/plots/High_Flux_HeatMap_31Jan22.pdf")

# Need to add a list of weak/bad amps
```

```python
# Now plot the noise heatmap
fig = plt.figure(figsize=(11,8.5))
ax = plt.axes([0.1, 0.1, 0.8, 0.8])
title = "Flux distribution in electrons.  Run 13144, SeqNo xxx"
plot_focal_plane(ax, noises, camera=camera, z_range=(0.0, 20.0), title=title)
plt.savefig("/repo/main/u/cslage/linearizer_28jan22/plots/Noise_HeatMap_31Jan22.pdf")

# Need to add a list of weak/bad amps
```

```python
# Now plot the gain heatmap
fig = plt.figure(figsize=(11,8.5))
ax = plt.axes([0.1, 0.1, 0.8, 0.8])
title = "Flux distribution in electrons.  Run 13144, SeqNo xxx"
plot_focal_plane(ax, gains, camera=camera, z_range=(1.0, 2.0), title=title)
plt.savefig("/repo/main/u/cslage/linearizer_28jan22/plots/Gain_HeatMap_31Jan22.pdf")

# Need to add a list of weak/bad amps
```

```python

```
