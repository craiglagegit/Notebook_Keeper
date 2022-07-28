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

Initially written 09 Mar 2020 by Craig Lage.\
Updated 18 Mar 20 to use empirical read noise

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
from lsst.ip.isr import PhotonTransferCurveDataset
```

```python jupyter={"outputs_hidden": true}
REPO_DIR = '/project/shared/auxTel/'
GAIN_DIR = '/project/shared/auxTel/rerun/cslage/PTC_Defect_2021-02-18/'
noise_file = '/project/cslage/AuxTel/noise/empirical_read_noise_17mar21.out'
butler = Butler(GAIN_DIR)
dataId = {'dayObs':'2021-02-18','detector':0, 'expId':2021021800100}
```

```python
numAmps = 16
# Get the empirical read noise file
file = open(noise_file, 'r')
noise_lines = file.readlines()
file.close()

# Get the yaml file
file = open('/home/cslage/alternate_branches/obs_lsst/policy/latiss/RXX.yaml', 'r')
lines = file.readlines()
file.close()

file = open('/home/cslage/alternate_branches/obs_lsst/policy/latiss/test.yaml', 'w')
# First, copy the header lines from the old file
for i in range(9):
    file.write(lines[i])

# Get the gain/noise data
datasetFile = GAIN_DIR+'/calibrations/ptc/ptcDataset-det000.fits'
ptc_data = PhotonTransferCurveDataset.readFits(datasetFile)
gain = ptc_data.gain
oldNoises = ptc_data.noise

raw = butler.get('raw', dataId=dataId)
ccd = raw.getDetector()

for amp in ccd:
    ampName = amp.getName()
    newGain = gain[ampName]
    oldNoise = oldNoises[ampName]
    noise = 0.0
    numNoises = 0
    for line in noise_lines:
        items = line.split(' ')
        thisAmpName = items[-2].strip(':')
        if ampName == thisAmpName:
            noise += float(items[-1].rstrip().strip('.'))
            numNoises += 1
    newNoise = noise / float(numNoises)
    print(ampName, newNoise, oldNoise)
    newLine = '      %s : { gain : %.4f, readNoise : %.1f }\n'%(ampName, newGain, newNoise)
    file.write(newLine)
file.close()

```

```python

```
