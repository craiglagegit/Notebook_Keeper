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

```python
from lsstDebug import getDebugFrame
```

```python
def debugFit(stepname, xVector, yVector, yModel, mask, ampName):
    """Debug method for linearity fitting.
    Parameters
    ----------
    stepname : `str`
        A label to use to check if we care to debug at a given
        line of code.
    xVector : `numpy.array`, (N,)
        The values to use as the independent variable in the
        linearity fit.
    yVector : `numpy.array`, (N,)
        The values to use as the dependent variable in the
        linearity fit.
    yModel : `numpy.array`, (N,)
        The values to use as the linearized result.
    mask : `numpy.array` [`bool`], (N,) , optional
        A mask to indicate which entries of ``xVector`` and
        ``yVector`` to keep.
    ampName : `str`
        Amplifier name to lookup linearity correction values.
    """
    fig, axs = plt.subplots(2)
    plt.subplots_adjust(hspace=1.0)

    if mask is None:
        mask = np.ones_like(xVector, dtype=bool)

    fig.suptitle(f"{stepname} {ampName} 'Spline'")
    if stepname == 'linearFit':
        axs[0].set_xlabel("Input Abscissa (time or mondiode)")
        axs[0].set_ylabel("Input Ordinate (flux)")
        axs[1].set_xlabel("Linear Ordinate (linear flux)")
        axs[1].set_ylabel("Flux Difference: (input - linear)")
    elif stepname in ('polyFit', 'splineFit'):
        axs[0].set_xlabel("Linear Abscissa (linear flux)")
        axs[0].set_ylabel("Input Ordinate (flux)")
        axs[1].set_xlabel("Linear Ordinate (linear flux)")
        axs[1].set_ylabel("Flux Difference: (input - full model fit)")
    elif stepname == 'solution':
        axs[0].set_xlabel("Input Abscissa (time or mondiode)")
        axs[0].set_ylabel("Linear Ordinate (linear flux)")
        axs[1].set_xlabel("Model flux (linear flux)")
        axs[1].set_ylabel("Flux Difference: (linear - model)")

    axs[0].set_yscale('log')
    axs[0].set_xscale('log')
    axs[0].scatter(xVector, yVector)
    axs[0].scatter(xVector[~mask], yVector[~mask], c='red', marker='x')
    axs[1].set_xscale('log')

    axs[1].scatter(yModel, yVector[mask] - yModel)
    fig.show()

    #plt.close()

```

```python
expId=3021120700200
det = 74
lin = butler.get('linearizer', detector=det, exposure=expId, instrument='LSSTCam')
linPtc = linPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
nonlinPtc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')

```

```python
from lsst.cp.pipe.utils import (funcPolynomial, irlsFit)
```

```python
lin.linearityType['C00']
```

```python

```

```python
ampName = 'C11'
linearityCoeffs = lin.linearityCoeffs[ampName]
minLinearAdu = 2000.0
maxLinearAdu = 20000.0
#mask = np.array(nonlinPtc.expIdMask[ampName], dtype=bool)
mask = np.repeat(True, len(nonlinPtc.expIdMask[ampName])) # if ignorePtcMask=True

inputAbscissa = np.array(nonlinPtc.rawExpTimes[ampName])[mask]
inputOrdinate = np.array(nonlinPtc.rawMeans[ampName])[mask]
#inputAbscissa = np.array(modExpTimes)[mask]
#inputOrdinate = np.array(fluxes)
```

```python
fluxMask = inputOrdinate < maxLinearAdu
lowMask = inputOrdinate > minLinearAdu
fluxMask = fluxMask & lowMask
linearAbscissa = inputAbscissa[fluxMask]
linearOrdinate = inputOrdinate[fluxMask]

linearFit, linearFitErr, chiSq, weights = irlsFit([0.0, 100.0], linearAbscissa,
                                                  linearOrdinate, funcPolynomial)
```

```python
debugFit('linearFit', inputAbscissa, inputOrdinate, linearOrdinate, fluxMask, ampName)
```

```python
nSigmaClipLinear = 5.0

# Convert this proxy-to-flux fit into an expected linear flux
linearOrdinate = linearFit[0] + linearFit[1] * inputAbscissa

# Exclude low end outliers
threshold = nSigmaClipLinear * np.sqrt(linearOrdinate)
fluxMask = np.abs(inputOrdinate - linearOrdinate) < threshold
linearOrdinate = linearOrdinate[fluxMask]
fitOrdinate = inputOrdinate[fluxMask]
```

```python
import lsst.afw.math as afwMath

fitOrder = 16 # Number of spline knots

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
debugFit('splineFit', binCenters, np.abs(values), values, None, ampName)

```

```python
plt.scatter(binCenters, values)
```

```python
interp = afwMath.makeInterpolate(binCenters.tolist(), values.tolist(),
                                 afwMath.stringToInterpStyle("AKIMA_SPLINE"))
modelOrdinate = linearOrdinate + interp.interpolate(linearOrdinate)
debugFit('splineFit', linearOrdinate, fitOrdinate, modelOrdinate, None, ampName)
```

```python
for i in range(fitOrder):
    print(binCenters[i], linearityCoeffs[i], values[i], linearityCoeffs[fitOrder+i])
```

```python
plt.subplot(1,1,1)
plt.scatter(inputAbscissa[fluxMask], (inputOrdinate[fluxMask] - linearOrdinate))
plt.scatter(binCenters / linearFit[1], values, marker = 'x', s = 200, color='red')
plt.plot(inputAbscissa[fluxMask], (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red')
```

```python
#Using monDiode
plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
plt.title("Spline fit to exposure data")
plt.plot(linearOrdinate, (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
plt.scatter(binCenters, values, marker = 'x', s = 200, color='red', label="Spline knots")
plt.scatter(linearOrdinate, (inputOrdinate[fluxMask] - linearOrdinate), label="Input data")
plt.xlabel("Flux (ADU)")
plt.ylabel("Deviation from Linearity(ADU)")
plt.legend()
plt.subplot(1,2,2)
plt.title("Residuals")
plt.scatter(linearOrdinate, (modelOrdinate - inputOrdinate[fluxMask]))
plt.xlabel("Flux (ADU)")
plt.ylabel("Residual (ADU)")
plt.ylim(-100,100)

```

```python
#Using ExpTimes
plt.figure(figsize=(12,6))
plt.subplot(1,2,1)
plt.title("Spline fit to exposure data")
plt.plot(linearOrdinate, (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
plt.scatter(binCenters, values, marker = 'x', s = 200, color='red', label="Spline knots")
plt.scatter(linearOrdinate, (inputOrdinate[fluxMask] - linearOrdinate), label="Input data")
plt.xlabel("Flux (ADU)")
plt.ylabel("Deviation from Linearity(ADU)")
plt.legend()
plt.subplot(1,2,2)
plt.title("Residuals")
plt.scatter(linearOrdinate, (modelOrdinate - inputOrdinate[fluxMask]))
plt.xlabel("Flux (ADU)")
plt.ylabel("Residual (ADU)")
plt.ylim(-100,100)

```

```python
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

Below is mondiode stuff

```python
expId=3021120700200
nonlinPtc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
expId = nonlinPtc.inputExpIdPairs['C00'][0][0][0]
expTime  = nonlinPtc.rawExpTimes['C00'][0]
```

```python
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
expId=3021120700200
nonlinPtc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
pairs = nonlinPtc.inputExpIdPairs['C00']
rawExpTimes  = nonlinPtc.rawExpTimes['C00']
rawMeans = nonlinPtc.rawMeans['C00']
```

```python
len(rawMeans)
```

```python
len(rawExpTimes)
```

```python
expTimes = []
monDiodes = []
fluxes = []
modExpTimes = []
for i, pair in enumerate(pairs):
    pair = pair[0]
    expTime = rawExpTimes[i]
    modExpTime = 0.0
    nExps = 0
    for j in range(2):
        expId = pair[j]
        try:
            monDiode = calcMondiode(expId)
            modExpTime += monDiode
            nExps += 1
            expTimes.append(expTime)
            monDiodes.append(monDiode)
            fluxes.append(rawMeans[i])
        except:
            continue
    if nExps > 0:
        # The 5E8 factor bring the modExpTimes back to about the same order as the expTimes
        modExpTime = 5.0E8 * modExpTime / nExps
    else:
        modExpTime = 0.0
    modExpTimes.append(modExpTime)
    #break
```

```python tags=[]
max(modExpTimes)
```

```python
len(modExpTimes)
```

```python
plt.scatter(expTimes, monDiodes)
```

```python
plt.scatter(expTimes, fluxes)
```

```python
plt.scatter(monDiodes, fluxes)
```

```python
plt.scatter(modExpTimes, fluxes)
```

```python

```
