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

# Notebook for taking gains extracted from PTC and writing them to a yaml file.

Initially written 09 Mar 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
! eups list -s obs_lsst
```

```python
from lsst.daf.persistence import Butler
import sys, os, glob
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
```

```python jupyter={"outputs_hidden": true}
REPO_DIR = '/project/shared/auxTel/'
GAIN_DIR = '/project/shared/auxTel/rerun/cslage/PTC_Defect_2021-02-18/'
butler = Butler(REPO_DIR)
dataId = {'dayObs':"2021-02-18", 'detector':0}
```

```python
numAmps = 16
# Get the yaml file
file = open('/home/cslage/alternate_branches/obs_lsst/policy/latiss/RXX.yaml', 'r')
lines = file.readlines()
file.close()

file = open('/home/cslage/alternate_branches/obs_lsst/policy/latiss/test.yaml', 'w')
# First, copy the header lines from the old file
for i in range(9):
    file.write(lines[i])

# Get the gain/noise data
gain_pickle_file = GAIN_DIR+'calibrations/ptc/ptcDataset-det000.pkl'
gain_file = open(gain_pickle_file, 'rb')
data = pkl.load(gain_file)
raw = butler.get('raw', dataId=dataId)
ccd = raw.getDetector()

for amp in ccd:
    ampName = amp.getName()
    newGain = data.gain[ampName]
    newNoise = data.noise[ampName]
    newLine = '      %s : { gain : %.4f, readNoise : %.1f }\n'%(ampName, newGain, newNoise)
    file.write(newLine)
file.close()

```

```python

```
