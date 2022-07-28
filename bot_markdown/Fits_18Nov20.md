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
from lsst.ip.isr import PhotonTransferCurveDataset
```

```python
MAIN_DIR = '/project/shared/BOT/rerun/cslage/'
dirs = ['PTC_LSSTCAM_12673T3/','PTC_LSSTCAM_12673T2/','PTC_LSSTCAM_12673T1/','PTC_LSSTCAM_FullCov_12673T/']
names = ['LeastSq', 'Cauchy ', 'Fair   ', 'FullCov']
DETECTOR = 183

```

```python
slacAmps = {'C10':'AMP01','C11':'AMP02','C12':'AMP03','C13':'AMP04',\
           'C14':'AMP05','C15':'AMP06','C16':'AMP07','C17':'AMP08',\
           'C07':'AMP09','C06':'AMP10','C05':'AMP11','C04':'AMP12',\
           'C03':'AMP13','C02':'AMP14','C01':'AMP15','C00':'AMP16'}
```

```python
gains =[]
noises = []
a00s = []
for amp in slacAmps.keys():
    print(amp)
    for i, thisDir in enumerate(dirs):
        name = names[i]
        datasetFile = MAIN_DIR + thisDir + '/calibrations/ptc/ptcDataset-det%03d.fits'%DETECTOR
        ptcDataset = PhotonTransferCurveDataset.readFits(datasetFile)
        gain = ptcDataset.gain[amp]
        noise = ptcDataset.noise[amp]
        if ptcDataset.ptcFitType == 'EXPAPPROXIMATION':
            a00 = ptcDataset.ptcFitPars[amp][0]
        if ptcDataset.ptcFitType == 'FULLCOVARIANCE':
            a00 = ptcDataset.aMatrix[amp][0][0]
        print("%s: Gain = %.6f, Noise = %.2f, a00 = %.6g"%(name,gain,noise,a00))
     
```

```python
print(dir(ptcDataset))
print(ptcDataset.aMatrix['C00'][0][0])
```

```python

```
