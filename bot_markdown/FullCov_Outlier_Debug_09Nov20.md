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

# Notebook for plotting FullCov results from ptc.py

Initially written 09 Nov 2019 by Craig Lage.\


```python
! eups list -s | grep lsst_distrib
! eups list -s ip_isr 
! eups list -s cp_pipe
```

```python
import matplotlib.pyplot as plt
from scipy.stats import norm
import numpy as np
import pickle as pkl
from scipy.stats import median_abs_deviation as mad
```

```python
# Get the pickled results from ptc.py
amp = 'C01'
counter = 1
nSigma = 20
maxFlux = 80
filename = "/project/cslage/BOT_LSSTCam/fullcov_tests/dump_94_12673_new_%s.pkl"%(amp)
file = open(filename, 'rb')
data = pkl.load(file)
file.close()
mu = data['mu']
cov = data['cov']
covModel = data['covModel']
sqrtW = data['sqrtW']
```

```python
for n in range(cov.shape[0]):
    print(n, mu[n], cov[n,0,0])#, covModel[n,0,0], sqrtW[n,0,0])
```

```python
wres = (covModel-cov)*sqrtW
flatwres = wres.flatten()
sig = mad(flatwres[flatwres != 0], scale='normal')
mask = (np.abs(flatwres) > (nSigma*sig))

nOutliers = mask.sum()
print(nOutliers)
print(mask.shape)
(nF, nx, ny) = wres.shape
```

```python
xmin = -10.0; xmax = 10.0
x = np.linspace(xmin, xmax, 100)

n, bins, patches = plt.hist(flatwres, bins=20, range=(xmin,xmax))
y = norm.pdf(x, 0.0, sig)*len(flatwres)
l = plt.plot(x, y, 'r--', linewidth=2)
print(sig)
print(np.std(flatwres))
print(cov[13,0,0], covModel[13,0,0], sqrtW[13,0,0])
print(wres[13,3,5], flatwres[13*nx*ny + 3*nx + 5])
```

```python
fig = plt.figure(figsize = (16,16))
plt.suptitle("FULLCOVARIANCE, Amp %s, Iteration %d"%(amp,counter), fontsize = 18)
plotNum = 0
for i in range(3):
    for j in range(3):
        plotNum += 1
        plt.subplot(3,3,plotNum)
        plt.title("Covariance, PixX = %d, PixY = %d"%(i,j))
        covData = []
        covFit = []    
        mus = []
        goodMus = []
        badMus = []
        maskedCovFit = []    

        for n in range(cov.shape[0]):
            mus.append(mu[n])
            covData.append(cov[n,i,j])
            if mask[n*nx*ny + i*nx + j]:
                maskedCovFit.append(covModel[n,i,j])
                badMus.append(mu[n])
            else:
                covFit.append(covModel[n,i,j])
                goodMus.append(mu[n])
        plt.scatter(mus, covData, marker='o', color='blue', label = 'Data')
        plt.scatter(goodMus, covFit, marker='x', s=200, color = 'green', label = 'Fit - kept')
        plt.scatter(badMus, maskedCovFit, marker='+', s=200, color = 'red', label = 'Fit, discarded')
        plt.xlabel("Flux(ADU)")
        plt.ylabel("Covariance(ADU^2)")
        plt.legend()
plt.savefig('/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_FullCov_12673/plots/Discard_Debug_94_%dK_S%d_%s_%d_09Nov20.pdf'\
            %(maxFlux, nSigma, amp, counter))
```

```python
from lsst.ip.isr import PhotonTransferCurveDataset
#datasetFile = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_FullCov_12606/calibrations/ptc/ptcDataset-det183.fits'
datasetFile = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_Test4_12606/calibrations/ptc/ptcDataset-det183.fits'
ptcDataset = PhotonTransferCurveDataset.readFits(datasetFile)
```

```python
cov = ptcDataset.covariances
covSqrtW = ptcDataset.covariancesSqrtWeights
print(cov['C13'][12][0][0])
print(covSqrtW['C13'][:][0][0])
for n, arr in enumerate(covSqrtW['C13']):
    print(n, cov['C13'][n][0][0], arr[0][0], sqrtW[n,0,0], mask[n*nx*ny + 0*nx + 0])

```

```python

```
