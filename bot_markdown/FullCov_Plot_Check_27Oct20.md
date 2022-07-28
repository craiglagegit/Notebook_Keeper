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

# Notebook for plotting extracted covariances on BOT data

Initially written 16 Oct 2019 by Craig Lage.\


```python
! eups list -s | grep lsst_distrib
! eups list -s obs_lsst 
! eups list -s cp_pipe
```

```python
import matplotlib.pyplot as plt
import numpy as np
from lsst.daf.persistence import Butler
```

```python
butler = Butler('/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_FullCov_12606')
ptcDataset = butler.get('photonTransferCurveDataset', dataId={'detector': 94})
print(ptcDataset.ptcFitType)
rawMeans = ptcDataset.rawMeans['C15']
rawVars = ptcDataset.rawVars['C15']
finalMeans = ptcDataset.finalMeans['C15']
finalVars = ptcDataset.finalVars['C15']
plt.plot(rawMeans,rawVars, marker='o')
plt.plot(finalMeans,finalVars, marker='+')
plt.xscale('log')
plt.yscale('log')
```

```python
# To remember, an alternate way to get the ptc data:
from lsst.ip.isr import PhotonTransferCurveDataset
datasetFile = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_FullCov_12606/calibrations/ptc/ptcDataset-det094.fits'
ptcDataset = PhotonTransferCurveDataset.readFits(datasetFile)
means = ptcDataset.finalMeans['C15']
vars = ptcDataset.finalVars['C15']
for i, mean in enumerate(means):
    print(i, mean, vars[i])
    means = ptcDataset.finalMeans['C15']
#rawMeans = ptcDataset.rawMeans['C15']
#for i, mean in enumerate(rawMeans):
#    print(i, mean)
```

```python
# From within plotPtc.py
means  = [   225.4469599,  329.21562534, 481.21701855, 703.85519784,
   1028.49397225,   2200.01217804,   3216.80205311,   4702.43936309,
   6875.21612947 ,  4972.59095203 ,  7266.06818337,  10622.48769831,
  15533.26499826 , 22708.70748936 , 33197.63705894 , 70740.36807742,
 103553.33074419]
vars = [  241.74304921 ,  334.78784098 ,  473.76934613,   673.71161401,
   968.91323279 , 2055.53412341 , 2950.21384776,  4298.62188172,
  4488.53009676 , 6210.96495104,  6531.54504735 , 9386.82295156,
 13450.7147755 , 19155.31113836 ,27034.00942317, 51094.59809799,
 67935.72149933]
for i, mean in enumerate(means):
    print(i, mean, vars[i])

plt.plot(means,vars, marker='x')
plt.xscale('log')
plt.yscale('log')
```

```python
gain = ptcDataset.gain['C15']
print(gain)
```

```python
print("i \t rawMean \t rawVar \t plotMean \t plotVar \t finalMean*g \t finalVar*g^2")

for i in range(len(rawMeans)):
    if i < 5:
        print("%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f"\
              %(i, rawMeans[i], rawVars[i], means[i], vars[i], finalMeans[i]*gain, finalVars[i] * gain**2))
    elif i == 5:
        print("%d \t %.2f \t %.2f "\
              %(i, rawMeans[i], rawVars[i]))
    elif i < 16:
        print("%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f"\
              %(i, rawMeans[i], rawVars[i], means[i-1], vars[i-1], finalMeans[i-1]*gain, finalVars[i-1] * gain**2))
    elif i == 16:
        print("%d \t %.2f \t %.2f "\
              %(i, rawMeans[i], rawVars[i]))
    else:
        print("%d \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f \t %.2f"\
              %(i, rawMeans[i], rawVars[i], means[i-2], vars[i-2], finalMeans[i-2]*gain, finalVars[i-2] * gain**2))

```

```python
# It's the covariances that are swapped.
from lsst.ip.isr import PhotonTransferCurveDataset
datasetFile = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_FullCov_12606/calibrations/ptc/ptcDataset-det183.fits'
ptcDataset = PhotonTransferCurveDataset.readFits(datasetFile)
rawMeans = ptcDataset.rawMeans['C15']
rawVars = ptcDataset.rawVars['C15']
Cov00 = np.array(ptcDataset.covariances['C15'])[:,0,0]
print("i \t rawMean \t rawVar \t Cov[0,0]")
for i in range(len(rawMeans)):
    if i==9 or i==10:
        print("%d \t %.2f \t %.2f \t %.2f **Swapped**"%(i, rawMeans[i], rawVars[i], Cov00[i]))
    else:
        print("%d \t %.2f \t %.2f \t %.2f "%(i, rawMeans[i], rawVars[i], Cov00[i]))
```

```python

```
