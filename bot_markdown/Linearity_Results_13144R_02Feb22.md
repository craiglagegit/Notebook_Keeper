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
                                                    "u/cslage/calib/13144/calib.20220107"])
camera = butler.get('camera', instrument='LSSTCam')
```

```python
linPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144R"])
```

```python
nonlinPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144M"])
```

```python
linButler = Butler("/repo/main", collections=["u/cslage/linearizer_28jan22"])
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

```python tags=[]
expId=3021120700200
pdf = PdfPages("/repo/main/u/cslage/bps_13144R/plots/Linearity_Results_13144R_Corrected_02Feb22.pdf")

names = ["E2V", "ITL"]
linNames = ["Not Linearized", "Linearized"]

for i, det in enumerate([55, 74]):
    linPtc = linPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    nonlinPtc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    lin = linButler.get('linearizer', detector=det, exposure=expId, instrument='LSSTCam')
    for amp in camera[0].getAmplifiers():
        ampName = amp.getName()
        if [det, ampName] not in [[55, 'C17'], [74, 'C01']]:
            continue
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
                print(f"{names[i]}-{det}-{ampName}-{linNames[n]} Gain={gain:.4f}, A00={a00:.6g}, Noise={noise:.2f}, Turnoff={maxDM:.2f}")
                yplot = ExpApprox(np.array(ptc.rawMeans[ampName]), gain, a00, noise)
                plt.subplot(2,4,2*n+1)
                plt.title(f"{names[i]} - {det} - {ampName}\n{linNames[n]} - {correct}")
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
                plt.title(f"{names[i]} - {det} - {ampName} PTC Residual\n{linNames[n]}")
                plt.scatter(ptc.rawMeans[ampName], yplot - ptc.rawVars[ampName], marker='+', label="Raw")
                plt.plot([maxDM, maxDM], [-1000, 1000], ls = '--', color='black', label = "PTC Turnoff")
                plt.xlim(0,100000)
                plt.xticks([0,25000,50000,75000,100000])
                plt.ylim(-1000,1000)
                plt.xlabel("Flux (ADU)")
                plt.ylabel("PTC Residual (ADU)")

        # Now get and plot the linearizer fit
        # This code is copied from cp_pipe/linearity.py
        filename = '/project/cslage/BOT_LSSTCam/linearizer/corrections_13144M_06jan22.pkl'
        infile = open(filename,'rb')
        abscissaCorrections = pkl.load(infile)
        for kk, correct in enumerate(["Uncorrected", "Corrected"]):
            modExpTimes = []
            for ii, pair in enumerate(linPtc.inputExpIdPairs[ampName]):
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


                # Get the photodiode correction                                                                                                                        
                try:
                    correction = np.nanmedian(abscissaCorrections[str(pair)])
                except:
                    correction = 0.0
                if correct == "Corrected":
                    modExpTimes.append(modExpTime + correction)
                else:
                    modExpTimes.append(modExpTime)                    
            inputAbscissa = np.array(modExpTimes)[mask]

            #inputAbscissa = np.array(nonlinPtc.rawExpTimes[ampName])[mask]

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
            plt.gcf().text(0.4, 0.94 - kk * 0.46, f"{names[i]} - {det} - {ampName} - {correct}", fontsize = 18)
            plt.subplot(2,4,1+kk*4)
            plt.title("Spline fit to MonDiode data")
            plt.plot(inputAbscissa, (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
            plt.scatter(binCenters / linearFit[1], values, marker = 'x', s = 200, color='red', label="Spline knots")
            plt.scatter(inputAbscissa, (inputOrdinate - linearOrdinate), label="Input data")
            plt.xlabel("Modified Exposure Time (sec)")
            plt.ylabel("Deviation from Linearity(ADU)")
            if names[i] == 'E2V':
                plt.ylim(-500,50)
            elif names[i] == 'ITL':
                plt.ylim(-100,200)

            plt.legend()    
            plt.subplot(2,4,2 + kk * 4)
            plt.title("Spline fit residual")
            plt.scatter(inputAbscissa, (modelOrdinate - inputOrdinate) / modelOrdinate * 100.0)
            plt.xlabel("Modified Exposure Time (sec)")
            plt.ylabel("Residual (%)") 
            plt.ylim(-0.1,0.1)
            plt.subplot(2,4,3 + kk * 4)
            plt.title("Spline fit to MonDiode data")
            plt.plot(linearOrdinate, (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
            plt.scatter(binCenters, values, marker = 'x', s = 200, color='red', label="Spline knots")
            plt.scatter(linearOrdinate, (inputOrdinate - linearOrdinate), label="Input data")
            if names[i] == 'E2V':
                plt.plot([maxDM, maxDM], [-400, 0], ls = '--', color='black', label = "PTC Turnoff")
                plt.ylim(-500,50)
            elif names[i] == 'ITL':
                plt.plot([maxDM, maxDM], [-50, 200], ls = '--', color='black', label = "PTC Turnoff")
                plt.ylim(-100,200)
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
        #print(f"Finished {det} {ampName}")
pdf.close()



```

```python

```
