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

# Notebook for plotting ComCam gains.

Initially written 14 Aug 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
! eups list -s | grep cp_pipe
```

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
```

```python
DIR_1 = '/home/cslage/ComCam/20191113/'
DIR_2 = '/project/shared/comCam/rerun/cslage/PTC_2020-08-13/'
RAFT = 'R22'

dirs = [DIR_1, DIR_2]
names = ['Tucson-2019-11-13', 'LaSerena-2020-08-13']
markers = ['o', 'x', '+', '*', '^', 'v']
```

```python
# First, get the old gain data.
#It is too old to get with the Butler, so we have to go directly to the pkl files
gain_pickle_file = DIR_1+'calibrations/ptc/ptcDataGainAndNoise-det000.pkl'
gain_file = open(gain_pickle_file, 'rb')
gain = pkl.load(gain_file)
gain_file.close()
gains_Tuc = gain['gain']
print(gains_Tuc)
```

```python
# Next get the new data
butler = Butler(DIR_2)
ptcDataset = butler.get('photonTransferCurveDataset', raftName=RAFT, detectorName='S11')
gains_LaSerena = ptcDataset.gain
print(gains_LaSerena['C10'])
```

```python
tucPlot = []
lasPlot = []
for amp in gains_Tuc.keys():
    tucPlot.append(gains_Tuc[amp][0])
    lasPlot.append(gains_LaSerena[amp])
print(tucPlot)
print(lasPlot)
```

```python
minGain = 0.8
maxGain = 1.2

plt.figure(figsize=(16,16))
plt.title("Gain Comparison COMCAM", fontsize = 18)
plt.scatter(tucPlot, lasPlot, marker='o', s = 100)
plt.xlim(minGain, maxGain)
plt.ylim(minGain, maxGain)
plt.xlabel("Tucson gains-2019-11-13", fontsize = 18)
plt.ylabel("La Serena gains-2020-08-13", fontsize = 18)
xplot = np.linspace(minGain, maxGain,100)
plt.plot(xplot, xplot, ls = '--', color='black', label = '1.00')
plt.savefig(DIR_2+"plots/Gain_Tucson_LaSerena_14Aug20.pdf")
```

```python

```
