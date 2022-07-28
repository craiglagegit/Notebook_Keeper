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
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib",\
                                                    "LSSTCam/calib/u/cslage/13144"])
camera = butler.get('camera', instrument='LSSTCam')
```

```python
linPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144E"])
```

```python
nonlinPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144B"])
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
```

```python jupyter={"outputs_hidden": true} tags=[]
expId=3021120700200
pdf = PdfPages("/repo/main/u/cslage/bps_13144E/plots/Linearity_Results_13144E_22Dec21.pdf")

plt.subplots_adjust(wspace = 0.5)


names = ["E2V", "ITL"]
linNames = ["Not Linearized", "Linearized"]

for i, det in enumerate([55, 74]):
    #lin = butler.get('linearizer', detector=det, exposure=expId, instrument='LSSTCam')
    linPtc = linPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    nonlinPtc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    for amp in camera[0].getAmplifiers():
        ampName = amp.getName()
        fig = plt.figure(figsize=(16,4))
        for n, ptc in enumerate([nonlinPtc, linPtc]):
            gain = ptc.gain[ampName]
            a00 = ptc.ptcFitPars[ampName][0]
            noise = ptc.noise[ampName]
            yplot = ExpApprox(np.array(ptc.rawMeans[ampName]), gain, a00, noise)
            plt.subplot(1,4,2*n+1)
            plt.title(f"{names[i]} - {det} - {ampName}\n{linNames[n]}")
            plt.scatter(ptc.rawMeans[ampName], ptc.rawVars[ampName], marker='+', label="Raw Data")
            plt.plot(ptc.rawMeans[ampName], yplot, ls = '--', color = 'red', label = 'ExpApprox')
            plt.legend()
            plt.xlim(0, 100000)
            plt.xlabel("Mean (ADU)")
            plt.ylabel("Variance (ADU)")
            #plt.ylim(30000, 40000)
            plt.subplot(1,4,2*n+2)
            plt.title(f"{names[i]} - {det} - {ampName} Residual\n{linNames[n]}")
            plt.scatter(ptc.rawMeans[ampName], yplot - ptc.rawVars[ampName], marker='+', label="Raw")
            plt.xlim(0,100000)
            plt.ylim(-1000,1000)
            plt.xlabel("Mean (ADU)")
            plt.ylabel("Residual (ADU)")
        pdf.savefig(fig)
        plt.close(fig)
        #print(f"Finished {det} {ampName}")
pdf.close()



```

```python

```
