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
from lsst.cp.pipe.utils import (funcPolynomial, irlsFit)
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib",\
                                                    "u/cslage/calib/13144/calib.20220103"])
camera = butler.get('camera', instrument='LSSTCam')
```

```python
ptcButler_1 = Butler("/repo/main", collections=["u/cslage/bps_13144M"])
```

```python
ptcButler_2 = Butler("/repo/main", collections=["u/cslage/ptc_13177"])
```

```python
def ExpApprox(mu, g, a00, n):
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
        
def calcMondiode(expId):
    factor = 5.0
    DATA_DIR = '/lsstdata/offline/teststand/BOT/storage/'
    date = int(expId/100000)
    seq = expId - date * 100000
    date = date - 10000000
    file = DATA_DIR + '%d/MC_C_%d_%06d/Photodiode_Readings_%d_%06d.txt'%(date,date,seq,date,seq)

    x, y = np.recfromtxt(file).transpose()
    # Threshold for finding baseline current values:                                                                                                                                                         
    ythresh = (min(y) + max(y))/factor + min(y)
    # Subtract the median of the baseline values to get a calibrated                                                                                                                                         
    # current.                                                                                                                                                                                               
    y -= np.median(y[np.where(y < ythresh)])
    integral = sum((y[1:] + y[:-1])/2*(x[1:] - x[:-1]))
    return integral

```

```python tags=[]
pdf = PdfPages("/repo/main/u/cslage/bps_13144N/plots/Linearity_LoFlux_13144N_02Feb21.pdf")

names = ["E2V", "ITL"]
linNames = ["Not Linearized", "Linearized"]

for i, det in enumerate([55, 74]):
    expId=3021120700200
    ptc_1 = ptcButler_1.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    expId=3021112900404
    ptc_2 = ptcButler_2.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    for amp in camera[0].getAmplifiers():
        ampName = amp.getName()
        if [det, ampName] not in [[55, 'C04']]:
            continue
        fig = plt.figure(figsize=(16,8))
        plt.subplots_adjust(wspace = 0.5, hspace = 0.5)
        for n, ptc in enumerate([ptc_1, ptc_2]):
            gain = ptc.gain[ampName]
            a00 = ptc.ptcFitPars[ampName][0]
            noise = ptc.noise[ampName]
            mask = np.array(ptc.expIdMask[ampName], dtype=bool)
            maxDM = np.max(np.array(ptc.rawMeans[ampName])[mask])
            print(f"{names[i]}-{det}-{ampName}-{linNames[n]} Gain={gain:.4f}, A00={a00:.6g}, Noise={noise:.2f}, Turnoff={maxDM:.2f}")
            yplot = ExpApprox(np.array(ptc.rawMeans[ampName]), gain, a00, noise)
            plt.subplot(2,2,2*n+1)
            plt.title(f"{names[i]} - {det} - {ampName}\n{linNames[n]}")
            plt.scatter(ptc.rawMeans[ampName], ptc.rawVars[ampName], marker='+', label="Raw Data")
            plt.plot(ptc.rawMeans[ampName], yplot, ls = '--', color = 'red', label = 'ExpApprox')
            #plt.plot([maxDM, maxDM], [0, 50000], ls = '--', color='black', label = "PTC Turnoff")
            plt.legend()
            plt.xlim(0, 1000)
            plt.xticks([0,250,500,750,1000])
            plt.xlabel("Flux (ADU)")
            plt.ylabel("Variance (ADU^2)")
            plt.ylim(0,1000)
            plt.subplot(2,2,2*n+2)
            plt.title(f"{names[i]} - {det} - {ampName} PTC Residual\n{linNames[n]}")
            plt.scatter(np.array(ptc.rawMeans[ampName]) * gain, np.array(yplot - ptc.rawVars[ampName]) * gain**2, marker='+', label="Raw")
            #plt.plot([maxDM, maxDM], [-1000, 1000], ls = '--', color='black', label = "PTC Turnoff")
            plt.xlim(0,1000)
            plt.xticks([0,250,500,750,1000])
            plt.ylim(-25, 0)
            plt.xlabel("Flux (e-)")
            plt.ylabel("PTC Residual (e-^2)")
        #pdf.savefig(fig)
        #plt.close(fig)            
pdf.close()



```

```python tags=[]
pdf = PdfPages("/repo/main/u/cslage/ptc_13177/plots/Linearity_LoFlux_02Feb21.pdf")

names = ["E2V", "ITL"]
ptcNames = ["13144", "13117"]

for i, det in enumerate([55, 74]):
    expId=3021120700200
    ptc_1 = ptcButler_1.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    expId=3021112900404
    ptc_2 = ptcButler_2.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    for amp in camera[0].getAmplifiers():
        ampName = amp.getName()
        #if [det, ampName] not in [[55, 'C01']]:
        #    continue
        fig = plt.figure(figsize=(12,4))
        plt.subplots_adjust(wspace = 0.5, hspace = 0.5)
        for n, ptc in enumerate([ptc_1, ptc_2]):
            gain = ptc.gain[ampName]
            a00 = ptc.ptcFitPars[ampName][0]
            noise = ptc.noise[ampName]
            mask = np.array(ptc.expIdMask[ampName], dtype=bool)
            maxDM = np.max(np.array(ptc.rawMeans[ampName])[mask])
            print(f"{names[i]}-{det}-{ampName}-{ptcNames[n]} Gain={gain:.4f}, A00={a00:.6g}, Noise={noise:.2f}, Turnoff={maxDM:.2f}")
            yplot = ExpApprox(np.array(ptc.rawMeans[ampName]), gain, a00, noise)
            plt.subplot(1,3,n+1)
            plt.title(f"{names[i]} - {det} - {ampName}\n{ptcNames[n]}")
            plt.scatter(ptc.rawMeans[ampName], ptc.rawVars[ampName], marker='+', label="Raw Data")
            plt.plot(ptc.rawMeans[ampName], yplot, ls = '--', color = 'red', label = 'ExpApprox')
            #plt.plot([maxDM, maxDM], [0, 50000], ls = '--', color='black', label = "PTC Turnoff")
            plt.legend()
            plt.xlim(0, 1000)
            plt.xticks([0,250,500,750,1000])
            plt.xlabel("Flux (ADU)")
            plt.ylabel("Variance (ADU^2)")
            plt.ylim(0,1000)
            plt.subplot(1,3,3)
            plt.title(f"{names[i]} - {det} - {ampName} PTC Residual")
            plt.scatter(np.array(ptc.rawMeans[ampName]) * gain, np.array(yplot - ptc.rawVars[ampName]) * gain**2, marker='+', label=ptcNames[n])
            #plt.plot([maxDM, maxDM], [-1000, 1000], ls = '--', color='black', label = "PTC Turnoff")
            plt.xlim(0,1000)
            plt.xticks([0,250,500,750,1000])
            plt.ylim(-50,0)
            plt.xlabel("Flux (e-)")
            plt.ylabel("PTC Residual (e-^2)")
        plt.subplot(1,3,3)
        plt.legend()

        pdf.savefig(fig)
        plt.close(fig)            
pdf.close()



```

```python

```
