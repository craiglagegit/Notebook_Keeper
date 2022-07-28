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
                                                    "LSSTCam/calib/u/cslage/13144"])
camera = butler.get('camera', instrument='LSSTCam')
```

```python
linPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144F"])
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

```python
# Set parameters
minLinearAdu = 2000.0
maxLinearAdu = 20000.0
nSigmaClipLinear = 5.0
fitOrder = 10 # Number of spline knots
```

```python
expId=3021120700200
pdf = PdfPages("/repo/main/u/cslage/bps_13144F/plots/Linearity_Results_13144F_23Dec21.pdf")

names = ["E2V", "ITL"]
linNames = ["Not Linearized", "Linearized"]

for i, det in enumerate([55, 74]):
    #lin = butler.get('linearizer', detector=det, exposure=expId, instrument='LSSTCam')
    linPtc = linPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    nonlinPtc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    for amp in camera[0].getAmplifiers():
        ampName = amp.getName()
        fig = plt.figure(figsize=(16,8))
        plt.subplots_adjust(wspace = 0.5, hspace = 0.5)
        for n, ptc in enumerate([nonlinPtc, linPtc]):
            gain = ptc.gain[ampName]
            a00 = ptc.ptcFitPars[ampName][0]
            noise = ptc.noise[ampName]
            yplot = ExpApprox(np.array(ptc.rawMeans[ampName]), gain, a00, noise)
            plt.subplot(2,4,2*n+1)
            plt.title(f"{names[i]} - {det} - {ampName}\n{linNames[n]}")
            plt.scatter(ptc.rawMeans[ampName], ptc.rawVars[ampName], marker='+', label="Raw Data")
            plt.plot(ptc.rawMeans[ampName], yplot, ls = '--', color = 'red', label = 'ExpApprox')
            plt.legend()
            plt.xlim(0, 100000)
            plt.xticks([0,25000,50000,75000,100000])
            plt.xlabel("Flux (ADU)")
            plt.ylabel("Variance (ADU)")
            #plt.ylim(30000, 40000)
            plt.subplot(2,4,2*n+2)
            plt.title(f"{names[i]} - {det} - {ampName} PTC Residual\n{linNames[n]}")
            plt.scatter(ptc.rawMeans[ampName], yplot - ptc.rawVars[ampName], marker='+', label="Raw")
            plt.xlim(0,100000)
            plt.xticks([0,25000,50000,75000,100000])
            plt.ylim(-1000,1000)
            plt.xlabel("Flux (ADU)")
            plt.ylabel("PTC Residual (ADU)")
            
        # Now get and plot the linearizer fit
        # This code is copied from cp_pipe/linearity.py
        mask = np.repeat(True, len(nonlinPtc.expIdMask[ampName])) # if ignorePtcMask=True
        inputAbscissa = np.array(nonlinPtc.rawExpTimes[ampName])[mask]
        inputOrdinate = np.array(nonlinPtc.rawMeans[ampName])[mask]
        fluxMask = inputOrdinate < maxLinearAdu
        lowMask = inputOrdinate > minLinearAdu
        fluxMask = fluxMask & lowMask
        linearAbscissa = inputAbscissa[fluxMask]
        linearOrdinate = inputOrdinate[fluxMask]
        linearFit, linearFitErr, chiSq, weights = irlsFit([0.0, 100.0], linearAbscissa,
                                                          linearOrdinate, funcPolynomial)
        # Convert this proxy-to-flux fit into an expected linear flux
        linearOrdinate = linearFit[0] + linearFit[1] * inputAbscissa
        # Exclude low end outliers
        threshold = nSigmaClipLinear * np.sqrt(linearOrdinate)
        fluxMask = np.abs(inputOrdinate - linearOrdinate) < threshold
        linearOrdinate = linearOrdinate[fluxMask]
        fitOrdinate = inputOrdinate[fluxMask]
        numPerBin, binEdges = np.histogram(linearOrdinate, bins=fitOrder)
        # Algorithm note: With the counts of points per
        # bin above, the next histogram calculates the
        # values to put in each bin by weighting each
        # point by the correction value.
        values = np.histogram(linearOrdinate, bins=fitOrder,
                              weights=(inputOrdinate[fluxMask] - linearOrdinate))[0]/numPerBin
        # After this is done, the binCenters are
        # calculated by weighting by the value we're
        # binning over.  This ensures that widely
        # spaced/poorly sampled data aren't assigned to
        # the midpoint of the bin (as could be done using
        # the binEdges above), but to the weighted mean of
        # the inputs.  Note that both histograms are
        # scaled by the count per bin to normalize what
        # the histogram returns (a sum of the points
        # inside) into an average.
        binCenters = np.histogram(linearOrdinate, bins=fitOrder,
                                  weights=linearOrdinate)[0]/numPerBin
        values = values[numPerBin > 0]
        binCenters = binCenters[numPerBin > 0]
        interp = afwMath.makeInterpolate(binCenters.tolist(), values.tolist(),
                                         afwMath.stringToInterpStyle("AKIMA_SPLINE"))
        modelOrdinate = linearOrdinate + interp.interpolate(linearOrdinate)
        plt.subplot(2,4,5)
        plt.title("Spline fit to exposure data")
        plt.plot(inputAbscissa[fluxMask], (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
        plt.scatter(binCenters / linearFit[1], values, marker = 'x', s = 200, color='red', label="Spline knots")
        plt.scatter(inputAbscissa[fluxMask], (inputOrdinate[fluxMask] - linearOrdinate), label="Input data")
        plt.xlabel("Exposure Time (sec)")
        plt.ylabel("Deviation from Linearity(ADU)")
        plt.legend()    
        plt.subplot(2,4,6)
        plt.title("Spline fit residual")
        plt.scatter(inputAbscissa[fluxMask], (modelOrdinate - inputOrdinate[fluxMask]))
        plt.xlabel("Exposure Time (sec)")
        plt.ylabel("Residual (ADU)")    
        plt.subplot(2,4,7)
        plt.title("Spline fit to exposure data")
        plt.plot(linearOrdinate, (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
        plt.scatter(binCenters, values, marker = 'x', s = 200, color='red', label="Spline knots")
        plt.scatter(linearOrdinate, (inputOrdinate[fluxMask] - linearOrdinate), label="Input data")
        plt.xlabel("Flux (ADU)")
        plt.ylabel("Deviation from Linearity(ADU)")
        plt.legend()        
        plt.subplot(2,4,8)
        plt.title("Spline fit residual")
        plt.scatter(linearOrdinate, (modelOrdinate - inputOrdinate[fluxMask]))
        plt.xlabel("Flux (ADU)")
        plt.ylabel("Residual (ADU)")      
         
        pdf.savefig(fig)
        plt.close(fig)
        #print(f"Finished {det} {ampName}")
pdf.close()



```

```python

```
