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
butler = Butler("/repo/main", collections=["LSSTComCam/raw/all","LSSTComCam/calib",\
                                                    "u/cslage/comcam/calib_202201218"])
camera = butler.get('camera', instrument='LSSTComCam')
```

```python
linPtcButler = Butler("/repo/main", collections=["u/cslage/comcam/ptc_linearized_20220218"])
```

```python
nonlinPtcButler = Butler("/repo/main", collections=["u/cslage/comcam/ptc_20220218"])
```

```python
linButler = Butler("/repo/main", collections=["u/cslage/comcam/linearizerA_20220218"])
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

```python tags=[]
expId=2022021800078
pdf = PdfPages("/repo/main/u/cslage/comcam/ptc_linearized_20220218/plots/Linearity_Results_20220218.pdf")

linNames = ["Not Linearized", "Linearized"]

for det in range(9):
    #if det > 0:
    #    continue
    linPtc = linPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTComCam')
    nonlinPtc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTComCam')
    lin = linButler.get('linearizer', detector=det, exposure=expId, instrument='LSSTComCam')
    for amp in camera[0].getAmplifiers():
        ampName = amp.getName()
        #if [det, ampName] not in [[0, 'C17']]:
        #    continue
        mask = np.array(linPtc.expIdMask[ampName], dtype=bool)
        maxDM = np.max(np.array(linPtc.rawMeans[ampName])[mask])            
        print(det, ampName, maxDM)
        fig = plt.figure(figsize=(16,8))
        plt.subplots_adjust(wspace = 0.5, hspace = 0.5)

        for n, ptc in enumerate([nonlinPtc, linPtc]):
            gain = ptc.gain[ampName]
            a00 = ptc.ptcFitPars[ampName][0]
            noise = ptc.noise[ampName]
            mask = np.array(ptc.expIdMask[ampName], dtype=bool)
            maxDM = np.max(np.array(ptc.rawMeans[ampName])[mask])
            print(f"Detector-{det}-{ampName}-{linNames[n]} Gain={gain:.4f}, A00={a00:.6g}, Noise={noise:.2f}, Turnoff={maxDM:.2f}")
            yplot = ExpApprox(np.array(ptc.rawMeans[ampName]), gain, a00, noise)
            plt.subplot(2,4,2*n+1)
            plt.title(f"Detector - {det} - {ampName}\n{linNames[n]}")
            plt.scatter(ptc.rawMeans[ampName], ptc.rawVars[ampName], marker='+', label="Raw Data")
            plt.plot(ptc.rawMeans[ampName], yplot, ls = '--', color = 'red', label = 'ExpApprox')
            plt.plot([maxDM, maxDM], [0, 50000], ls = '--', color='black', label = "PTC Turnoff")
            plt.legend()
            plt.xlim(0, 100000)
            plt.xticks([0,25000,50000,75000,100000])
            plt.xlabel("Flux (ADU)")
            plt.ylabel("Variance (ADU)")
            #plt.ylim(30000, 40000)
            plt.subplot(2,4,2*n+2)
            plt.title(f"Detector - {det} - {ampName} PTC Residual\n{linNames[n]}")
            plt.scatter(ptc.rawMeans[ampName], yplot - ptc.rawVars[ampName], marker='+', label="Raw")
            plt.plot([maxDM, maxDM], [-1000, 1000], ls = '--', color='black', label = "PTC Turnoff")
            plt.xlim(0,100000)
            plt.xticks([0,25000,50000,75000,100000])
            plt.ylim(-1000,1000)
            plt.xlabel("Flux (ADU)")
            plt.ylabel("PTC Residual (ADU)")

        # Now get and plot the linearizer fit
        kk = 1
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
        # Get the spline coordinates from the stored linearizer
        binCenters, values = np.split(lin.linearityCoeffs[amp.getName()], 2)

        interp = afwMath.makeInterpolate(binCenters.tolist(), values.tolist(),
                                         afwMath.stringToInterpStyle("AKIMA_SPLINE"))
        modelOrdinate = linearOrdinate + interp.interpolate(linearOrdinate)
        plt.gcf().text(0.4, 0.94 - kk * 0.46, f"Detector - {det} - {ampName}", fontsize = 18)
        plt.subplot(2,4,1+kk*4)
        plt.title("Spline fit to ExpTime data")
        plt.plot(inputAbscissa, (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
        plt.scatter(binCenters / linearFit[1], values, marker = 'x', s = 200, color='red', label="Spline knots")
        plt.scatter(inputAbscissa, (inputOrdinate - linearOrdinate), label="Input data")
        plt.xlabel("Exposure Time (sec)")
        plt.ylabel("Deviation from Linearity(ADU)")
        plt.ylim(-200,200)

        plt.legend()    
        plt.subplot(2,4,2 + kk * 4)
        plt.title("Spline fit residual")
        plt.scatter(inputAbscissa, (modelOrdinate - inputOrdinate) / modelOrdinate * 100.0)
        plt.xlabel("Exposure Time (sec)")
        plt.ylabel("Residual (%)") 
        plt.ylim(-0.1,0.1)
        plt.subplot(2,4,3 + kk * 4)
        plt.title("Spline fit to ExpTime data")
        plt.plot(linearOrdinate, (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
        plt.scatter(binCenters, values, marker = 'x', s = 200, color='red', label="Spline knots")
        plt.scatter(linearOrdinate, (inputOrdinate - linearOrdinate), label="Input data")
        plt.plot([maxDM, maxDM], [-50, 200], ls = '--', color='black', label = "PTC Turnoff")
        plt.ylim(-200,200)
        plt.xlabel("Flux (ADU)")
        plt.ylabel("Deviation from Linearity(ADU)")
        plt.xlim(0, 100000)
        plt.xticks([0,25000,50000,75000,100000])
        plt.legend()        
        plt.subplot(2,4,4 + kk * 4)
        plt.title("Spline fit residual")
        plt.scatter(linearOrdinate, (modelOrdinate - inputOrdinate) / modelOrdinate * 100.0)
        plt.plot([maxDM, maxDM], [-0.1, 0.1], ls = '--', color='black', label = "PTC Turnoff")
        plt.xlabel("Flux (ADU)")
        plt.ylabel("Residual (%)")
        plt.ylim(-0.1,0.1)
        plt.xlim(0, 100000)
        plt.xticks([0,25000,50000,75000,100000])


    pdf.savefig(fig)
    plt.close(fig)
    print(f"Finished {det} {ampName}")
pdf.close()



```

```python

```
