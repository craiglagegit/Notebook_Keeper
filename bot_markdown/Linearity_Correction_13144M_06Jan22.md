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
```

```python
linPtcButler = Butler("/repo/main", collections=["u/cslage/bps_13144N"])
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
    raftCol = int(list(raft)[2])
    sensorRow = int(list(sensor)[1])
    sensorCol = int(list(sensor)[2])
    detectorNum = (rows[raftRow] + (raftCol - startingCol[raftRow])) * 9
    detectorNum += 3 * sensorRow + sensorCol
    return vendor, detectorNum, 4 - raftRow, raftCol

rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

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

abscissaCorrections = {}
linNames = ["Not Linearized", "Linearized"]
for RAFT in ['R22']:#rafts:
    for SENSOR in ['S11']:#sensors:
        VENDOR, det, raftRow, raftCol = detector(RAFT, SENSOR)
        linPtc = linPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
        nonlinPtc = nonlinPtcButler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
        lin = linPtcButler.get('linearizer', detector=det, exposure=expId, instrument='LSSTCam')
        #pdf = PdfPages(f"/repo/main/u/cslage/bps_13144N/plots/Linearity_det{det:03d}.pdf")

        for ampName in nonlinPtc.gain.keys():
            fig = plt.figure(figsize=(16,8))
            plt.subplots_adjust(wspace = 0.5, hspace = 0.5)
            for n, ptc in enumerate([nonlinPtc, linPtc]):
                gain = ptc.gain[ampName]
                a00 = ptc.ptcFitPars[ampName][0]
                noise = ptc.noise[ampName]
                mask = np.array(ptc.expIdMask[ampName], dtype=bool)
                maxDM = np.max(np.array(ptc.rawMeans[ampName])[mask])
                yplot = ExpApprox(np.array(ptc.rawMeans[ampName]), gain, a00, noise)
                plt.subplot(2,4,2*n+1)
                plt.title(f"{VENDOR} - {det} - {ampName}\n{linNames[n]}")
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
                plt.title(f"{VENDOR} - {det} - {ampName} PTC Residual\n{linNames[n]}")
                plt.scatter(ptc.rawMeans[ampName], yplot - ptc.rawVars[ampName], marker='+', label="Raw")
                plt.plot([maxDM, maxDM], [-1000, 1000], ls = '--', color='black', label = "PTC Turnoff")
                plt.xlim(0,100000)
                plt.xticks([0,25000,50000,75000,100000])
                plt.ylim(-1000,1000)
                plt.xlabel("Flux (ADU)")
                plt.ylabel("PTC Residual (ADU)")

            # Now get and plot the linearizer fit
            # This code is copied from cp_pipe/linearity.py

            modExpTimes = []
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
            binCenters, values = np.split(lin.linearityCoeffs[ampName], 2)

            interp = afwMath.makeInterpolate(binCenters.tolist(), values.tolist(),
                                             afwMath.stringToInterpStyle("AKIMA_SPLINE"))
            modelOrdinate = linearOrdinate + interp.interpolate(linearOrdinate)
            # Now calculate and store the correction
            fluxResidual = (inputOrdinate - modelOrdinate)
            abscissaCorrection = fluxResidual / linearFit[1]
            maskCounter = 0
            for ii, pair in enumerate(ptc.inputExpIdPairs[ampName]):
                key = str(pair[0])
                
                try:
                    if mask[ii]:
                        abscissaCorrections[key].append(abscissaCorrection[ii - maskCounter])
                    else:
                        abscissaCorrections[key].append(np.nan)
                        maskCounter += 1
                except KeyError:
                    abscissaCorrections[key] = []
                    if mask[ii]:
                        abscissaCorrections[key].append(abscissaCorrection[ii - maskCounter])
                    else:
                        abscissaCorrections[key].append(np.nan)
                        maskCounter += 1

            plt.subplot(2,4,5)
            plt.title("Spline fit to MonDiode data")
            plt.plot(inputAbscissa, (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
            plt.scatter(binCenters / linearFit[1], values, marker = 'x', s = 200, color='red', label="Spline knots")
            plt.scatter(inputAbscissa, (inputOrdinate - linearOrdinate), label="Input data")
            plt.xlabel("Exposure Time (sec)")
            plt.ylabel("Deviation from Linearity(ADU)")
            plt.legend()    
            plt.subplot(2,4,6)
            plt.title("Spline fit residual")
            plt.scatter(inputAbscissa, (modelOrdinate - inputOrdinate) / modelOrdinate * 100.0)
            plt.xlabel("Exposure Time (sec)")
            plt.ylabel("Residual (%)") 
            plt.ylim(-0.2,0.2)
            plt.subplot(2,4,7)
            plt.title("Spline fit to MonDiode data")
            plt.plot(linearOrdinate, (modelOrdinate-linearOrdinate), ls = '--', lw=3, color='red', label="Spline fit")
            plt.scatter(binCenters, values, marker = 'x', s = 200, color='red', label="Spline knots")
            plt.scatter(linearOrdinate, (inputOrdinate - linearOrdinate), label="Input data")
            plt.plot([maxDM, maxDM], [-200, 0], ls = '--', color='black', label = "PTC Turnoff")
            plt.xlabel("Flux (ADU)")
            plt.ylabel("Deviation from Linearity(ADU)")
            plt.xlim(0, 100000)
            plt.xticks([0,25000,50000,75000,100000])
            plt.legend()        
            plt.subplot(2,4,8)
            plt.title("Spline fit residual")
            plt.scatter(linearOrdinate, (modelOrdinate - inputOrdinate) / modelOrdinate * 100.0)
            plt.plot([maxDM, maxDM], [-200, 0], ls = '--', color='black', label = "PTC Turnoff")
            plt.xlabel("Flux (ADU)")
            plt.ylabel("Residual (%)")
            plt.ylim(-0.2,0.2)
            plt.xlim(0, 100000)
            plt.xticks([0,25000,50000,75000,100000])


            #pdf.savefig(fig)
            plt.close(fig)
            #print(f"Finished {det} {ampName}")
    #pdf.close()



```

```python
print(len(mask), len(ptc.inputExpIdPairs[ampName]), len(fluxResidual), len(modExpTimes), len(abscissaCorrections))
```

```python jupyter={"outputs_hidden": true} tags=[]
abscissaCorrections
```

```python
# Exploring corrections to the Monitor diode
```

```python
fluxResidual = (inputOrdinate - modelOrdinate)
```

```python
plt.plot(fluxResidual)
```

```python
plt.figure(figsize = (16,4))
plt.scatter(inputAbscissa, fluxResidual)
```

```python
linearFit[1]
```

```python
# This is how much to correct the modified exposure times.
# Correction is ~ .05/100 = .05%
plt.figure(figsize = (16,4))
plt.scatter(inputAbscissa, fluxResidual / linearFit[1])
```

```python
len(inputAbscissa)
```

```python
# Get the stored corrections
filename = '/project/cslage/BOT_LSSTCam/linearizer/corrections_13144M_06jan22.pkl'
infile = open(filename,'rb')
abscissaCorrections = pkl.load(infile)
infile.close()
```

```python
len(abscissaCorrections)
```

```python
maskedCorrections = []
for ii, pair in enumerate(ptc.inputExpIdPairs[ampName]):
    key = str(pair[0])
    if mask[ii]:
        maskedCorrections.append(abscissaCorrections[key])

```

```python
len(maskedCorrections)
```

```python
num = 304
plt.hist(maskedCorrections[num], bins=50, range=(-0.10, 0.10))
plt.text(-0.10, 500, f"Median = {np.nanmedian(maskedCorrections[num]):.4f}")
```

```python
print(RAFT, SENSOR, det, ampName)
for i in range(302, 308):
    print(inputAbscissa[i], (fluxResidual / linearFit[1])[i], np.nanmedian(maskedCorrections[i]))

```

```python
#correctedAbscissa = inputAbscissa + fluxResidual / linearFit[1]
```

```python
maskedCorrections = []
for ii, pair in enumerate(ptc.inputExpIdPairs[ampName]):
    key = str(pair[0])
    if mask[ii]:
        maskedCorrections.append(np.nanmedian(abscissaCorrections[key]))
maskedCorrections = np.array(maskedCorrections)
```

```python
correctedAbscissa = inputAbscissa + maskedCorrections
```

```python
plt.figure(figsize = (16,4))
plt.scatter(correctedAbscissa, inputOrdinate)
```

```python
#fitOrder = 30
fluxMask = inputOrdinate < maxLinearAdu
lowMask = inputOrdinate > minLinearAdu
fluxMask = fluxMask & lowMask

linearAbscissa = correctedAbscissa[fluxMask]
linearOrdinate = inputOrdinate[fluxMask]
linearFit, linearFitErr, chiSq, weights = irlsFit([0.0, 100.0], linearAbscissa,
                                                  linearOrdinate, funcPolynomial)
# Convert this proxy-to-flux fit into an expected linear flux
linearOrdinate = linearFit[0] + linearFit[1] * correctedAbscissa
threshold = nSigmaClipLinear * np.sqrt(linearOrdinate)
fluxMask = np.abs(inputOrdinate - linearOrdinate) < threshold
linearOrdinate = linearOrdinate[fluxMask]
numPerBin, binEdges = np.histogram(linearOrdinate, bins=fitOrder)
values = np.histogram(linearOrdinate, bins=fitOrder,
                      weights=(inputOrdinate[fluxMask] - linearOrdinate))[0]/numPerBin


binCenters = np.histogram(linearOrdinate, bins=fitOrder,
                          weights=linearOrdinate)[0]/numPerBin
values = values[numPerBin > 0]
binCenters = binCenters[numPerBin > 0]

interp = afwMath.makeInterpolate(binCenters.tolist(), values.tolist(),
                             afwMath.stringToInterpStyle("AKIMA_SPLINE"))

modelOrdinate = linearOrdinate + interp.interpolate(linearOrdinate)
newFluxResidual = (inputOrdinate[fluxMask] - modelOrdinate)
plotMask = newFluxResidual < 0.1
plt.figure(figsize = (16,4))
plt.scatter(correctedAbscissa, newFluxResidual)
```

```python
# This reduced the residuals by about a factor of 5.
```
