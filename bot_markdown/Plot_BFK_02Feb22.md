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

# Notebook for plotting BF Kernels generated with Gen 3.

Initially written 22 Jun 2021 by Craig Lage.

```python
import os, sys, time, datetime, glob, subprocess
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pf
from scipy import stats
from lsst.daf.butler import Butler
```

```python
defaultButler = Butler('/repo/main', collections=
                ['LSSTCam/raw/all', 'LSSTCam/calib', 'u/cslage/tests/bfkB_01feb22'])
zeroButler = Butler('/repo/main', collections=
                ['LSSTCam/raw/all', 'LSSTCam/calib', 'u/cslage/tests/bfkD_01feb22'])
quadButler = Butler('/repo/main', collections=
                ['LSSTCam/raw/all', 'LSSTCam/calib', 'u/cslage/tests/bfkE_01feb22'])

names = ['Default', 'ForceZeroSum', 'Quad']
butlers = [defaultButler, zeroButler, quadButler]
```

```python
expId = 3021120600576 # I think any exposure within the set of flat pairs will work.
DETECTOR = 55
dataId={'instrument':'LSSTCam', 'detector':DETECTOR, 'exposure':expId}
```

```python
# Now plot the correlations and the kernel. 
# So much variation in the (0,0) covariance!

amp = 'C14'
fig = plt.figure(figsize=(16,16))
plt.subplots_adjust(hspace=0.2)

for ii, butler in enumerate(butlers):
    print(ii, names[ii])
    ptc_dataset = butler.get('ptc', dataId=dataId)
    bf_kernel = butler.get('bfk', dataId=dataId)
    gains = bf_kernel.gain
    means = bf_kernel.means # Mean flux of flat pairs in electrons
    rawMeans = ptc_dataset.rawMeans # Mean flux of flat pairs in electrons
    rawVars = ptc_dataset.rawVars # Mean flux of flat pairs in electrons
    rawXcorrs = bf_kernel.rawXcorrs # Raw extracted covariances in ADU^2. [0,0] is the variance}
    meanXcorrs = bf_kernel.meanXcorrs # Extracted covariances used to extract kernel. These are per e-.
    kernels = bf_kernel.ampKernels # ampwise kernel
    ptcResults = ptc_dataset.ptcFitPars        
    plt.suptitle("COVARIANCES(*1E7)                          KERNEL(*1E7)", fontsize=24)
    plt.subplot(3,4,4*ii+1)
    plt.imshow(np.log10(abs(np.array(meanXcorrs[amp]))))
    plt.subplot(3,4,4*ii+2)
    plt.title("     Amp %s   %s"%(amp, names[ii]), fontsize=18)
    plt.plot([0,16],[0,0], ls='--', color='black')
    plt.plot(-meanXcorrs[amp][:,8]*1E7, color='blue', drawstyle='steps-mid')
    plt.plot(-meanXcorrs[amp][8,:]*1E7, linestyle='--', color='red', drawstyle='steps-mid')
    plt.ylim(-40,10)
    plt.subplot(3,4,4*ii+3)
    plt.imshow(kernels[amp])
    plt.subplot(3,4,4*ii+4)  
    plt.plot([0,16],[0,0], ls='--', color='black')
    plt.plot(kernels[amp][:,8]*1E7, color='blue', drawstyle='steps-mid')
    plt.plot(kernels[amp][8,:]*1E7, linestyle='--', color='red', drawstyle='steps-mid')
    plt.ylim(-20,2)
plt.savefig("/repo/main/u/cslage/tests/bfkB_01feb22/plots/BF_Kernel_Tests_02Feb22.pdf")
```

```python

```

```python

```
