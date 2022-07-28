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

# Notebook for interrogating FULLCOV data.

Initially written 06 Nov 2020 by Craig Lage

```python
! eups list -s | grep lsst_distrib
! eups list -s cp_pipe
```

```python
import sys, os, glob, subprocess
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
```

```python
run = '12673'
DATA_DIR = '/project/shared/BOT/'
REPO_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_FullCov_12673T'
DETECTOR = 1
#butler = Butler(REPO_DIR)
#ptcDataset = butler.get('photonTransferCurveDataset', dataId={'detector': DETECTOR})
from lsst.ip.isr import PhotonTransferCurveDataset
datasetFile = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_FullCov_12673T/calibrations/ptc/ptcDataset-det001.fits'
ptcDataset = PhotonTransferCurveDataset.readFits(datasetFile)
```

```python
print(ptcDataset.badAmps)
print(ptcDataset.expIdMask['C01'])
```

```python

cov = ptcDataset.covariances
covSqrtW = ptcDataset.covariancesSqrtWeights
finalVars = ptcDataset.finalVars
finalMeans = ptcDataset.finalMeans
```

```python
amp = 'C07'

print(cov[amp][12][0][0])
print(covSqrtW[amp][:][0][0])
for n, arr in enumerate(covSqrtW[amp]):
    print(n, finalMeans[amp][n], finalVars[amp][n], cov[amp][n][0][0])
#arr = np.array(covSqrtW[amp])
#print(arr.shape)
#print(arr[:,0,0])
#print(np.nanmax(finalMeans[amp]))
```

```python
print(len(finalVars[amp]), len(cov[amp][:][0][0]))
```

```python
plt.scatter(finalMeans[amp], finalVars[amp], marker='x', color='green')
plt.scatter(finalMeans[amp], np.array(cov[amp])[:,0,0], marker='o', color='red')
```

```python

```
