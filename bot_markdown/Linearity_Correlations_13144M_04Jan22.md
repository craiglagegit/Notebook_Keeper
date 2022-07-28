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
nonlinPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144M"])
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

```python
# Set parameters
minLinearAdu = 2000.0
maxLinearAdu = 20000.0
nSigmaClipLinear = 5.0
fitOrder = 10 # Number of spline knots
```

```python
filename = '/project/cslage/BOT_LSSTCam/linearizer/corrections_13144M_06jan22.pkl'
infile = open(filename,'rb')
abscissaCorrections = pkl.load(infile)
infile.close()
```

```python
expId=3021120700200

names = ["E2V", "ITL"]

offset = 0.15
fig = plt.figure(figsize=(8,8))

for i, det in enumerate([55]):#, 74]):
    ptc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    lin = butler.get('linearizer', detector=det, exposure=expId, instrument='LSSTCam')
    plt.subplot(2,1,i+1)
    plt.title(f"Spline fit residual - {names[i]} - {det}") 
    #plt.subplots_adjust(wspace = 0.5, hspace = 0.5)

    for jj, amp in enumerate(camera[0].getAmplifiers()):
        ampName = amp.getName()
        gain = ptc.gain[ampName]
        a00 = ptc.ptcFitPars[ampName][0]
        noise = ptc.noise[ampName]
        mask = np.array(ptc.expIdMask[ampName], dtype=bool)
        maxDM = np.max(np.array(ptc.rawMeans[ampName])[mask])
        # Now get and plot the linearizer fit
        # This code is copied from cp_pipe/linearity.py
        if ampName not in ['C10', 'C11', 'C12', 'C13']:
            continue
        modExpTimes = []
        corrections = []
        for ii, pair in enumerate(ptc.inputExpIdPairs[ampName]):
            pair = pair[0]
            modExpTime = 0.0
            nExps = 0
            for j in range(2):
                expId = pair[j]
                try:
                    monDiode = calcMondiode(expId)
                    modExpTime += monDiode
                    nExps += 1
                except:
                    continue
            if nExps > 0:
                # The 5E8 factor bring the modExpTimes back to about the same order as the expTimes                                                                                                          
                modExpTime = 5.0E8 * modExpTime / nExps
            else:
                mask[ii] = False
            try:
                correction = np.nanmedian(abscissaCorrections[str(pair)])
            except:
                correction = 0.0
                mask[ii] = False
            corrections.append(correction)
            modExpTimes.append(modExpTime)
        inputAbscissa = np.array(modExpTimes)[mask]        
        inputOrdinate = np.array(ptc.rawMeans[ampName])[mask]
        
        fluxMask = inputOrdinate < maxLinearAdu
        lowMask = inputOrdinate > minLinearAdu
        fluxMask = fluxMask & lowMask
        linearAbscissa = inputAbscissa[fluxMask]
        linearOrdinate = inputOrdinate[fluxMask]
        linearFit, linearFitErr, chiSq, weights = irlsFit([0.0, 100.0], linearAbscissa,
                                                          linearOrdinate, funcPolynomial)
        # Convert this proxy-to-flux fit into an expected linear flux
        linearOrdinate = linearFit[0] + linearFit[1] * inputAbscissa
        # Get the spline coordinates from the stored linearizer
        binCenters, values = np.split(lin.linearityCoeffs[amp.getName()], 2)

        interp = afwMath.makeInterpolate(binCenters.tolist(), values.tolist(),
                                         afwMath.stringToInterpStyle("AKIMA_SPLINE"))
        modelOrdinate = linearOrdinate + interp.interpolate(linearOrdinate)
        
        corrections = np.array(corrections)[mask]
        corrections =  - corrections * linearFit[1] / modelOrdinate * 100.0
        
        plotMask = (inputOrdinate - linearOrdinate) < 1000.0
        plt.scatter(inputAbscissa, (modelOrdinate - inputOrdinate) / modelOrdinate * 100.0 + offset * jj, label = ampName)
        plt.plot([0,100],[offset * jj, offset * jj], ls = '--', color='black')
        if ampName == 'C10':
            plt.scatter(inputAbscissa, corrections - 2 * offset, label = "Correction")
            plt.plot([0,100],[-offset * 2, -offset * 2], ls = '--', color='black')
        plt.xlabel("Effective exposure time(sec)")
        plt.ylabel("Residual (%) + Offset")
        plt.ylim(-.50, 0.6)
        plt.xlim(0, 100)
        #plt.xticks([0,25000,50000,75000,100000])
        plt.legend(loc = 'lower center', bbox_to_anchor = (0.5, -0.6))
        
         
plt.savefig("/repo/main/u/cslage/bps_13144M/plots/Linearity_Corrections_13144M_07Jan21.pdf")



```

```python tags=[]
for ii, pair in enumerate(ptc.inputExpIdPairs[ampName]):
    pair = pair[0]
    correction = np.nanmedian(abscissaCorrections[str(pair)])
    print(pair, correction)
```

```python

```
