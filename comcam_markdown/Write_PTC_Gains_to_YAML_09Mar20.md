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

```python
REPO_DIR = '/project/cslage/ComCam/20200303'
GAIN_DIR = '/home/cslage/ComCam/20191113/'
raftName = 'R22'
butler = Butler(REPO_DIR)
visit = 3020030300034

```

```python
numCCDs = 9
numAmps = 16
# Get the yaml file
file = open('/home/cslage/alternate_branches/obs_lsst/policy/comCam/R22.yaml', 'r')
lines = file.readlines()
file.close()


file = open('/home/cslage/alternate_branches/obs_lsst/policy/comCam/test.yaml', 'w')
# First, copy the header lines from the old file
for i in range(15):
    file.write(lines[i])

# Now loop through the detectors, correcting the gain and noise
for detector in range(numCCDs):
    # Get the gain/noise data
    gain_pickle_file = GAIN_DIR+'calibrations/ptc/ptcDataGainAndNoise-det%03d.pkl'%detector
    gain_file = open(gain_pickle_file, 'rb')
    gain_data = pkl.load(gain_file)
    raw = butler.get('raw', detector=detector, visit=visit)
    ccd = raw.getDetector()
    ccdName = ccd.getName()
    newLine = '    %s :\n'%ccdName.split('_')[1]
    file.write(newLine)
    for amp in ccd:
        ampName = amp.getName()
        newGain = gain_data['gain'][ampName][0]
        newNoise = gain_data['noise'][ampName][0]
        newLine = '      %s : { gain : %.4f, readNoise : %.1f }\n'%(ampName, newGain, newNoise)
        file.write(newLine)
file.close()

```

```python

```
