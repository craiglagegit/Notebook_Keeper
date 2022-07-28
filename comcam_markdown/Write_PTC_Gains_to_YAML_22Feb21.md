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
REPO_DIR = '/project/shared/comCam-CCS'
GAIN_DIR = '/project/shared/comCam-CCS/rerun/cslage/PTC_2020-12-29/'
NOISE_DIR = '/project/cslage/ComCam/noise/'
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
    # Get the empirical read noise file
    noise_filename = NOISE_DIR + 'empirical_read_noise_%d_18mar21.out'%detector
    noise_file = open(noise_filename, 'r')
    noise_lines = noise_file.readlines()
    noise_file.close()

    # Get the gain/noise data
    datasetFile = GAIN_DIR+'/calibrations/ptc/ptcDataset-det%03d.fits'%detector
    #ptc_data = PhotonTransferCurveDataset.readFits(datasetFile)
    hdulist = pf.open(datasetFile, mode='readonly', do_not_scale_image_data=True)
    data=hdulist[1].data

    gain_data = data['GAIN']
    old_noise_data = data['NOISE']
    ampName_data = data['AMPLIFIER_NAME']
    raw = butler.get('raw', detector=detector, visit=visit)
    ccd = raw.getDetector()
    ccdName = ccd.getName()
    newLine = '    %s :\n'%ccdName.split('_')[1]
    file.write(newLine)
    for amp in ccd:
        ampName = amp.getName()
        noise = 0.0
        numNoises = 0
        for line in noise_lines:
            items = line.split(' ')
            thisAmpName = items[-2].strip(':')
            if ampName == thisAmpName:
                noise += float(items[-1].rstrip().strip('.'))
                numNoises += 1
        newNoise = noise / float(numNoises)

        for i in range(16):
            fitsAmpName = ampName_data[i]
            newGain = gain_data[i]
            oldNoise = old_noise_data[i]
            if fitsAmpName == ampName:
                print(detector, ampName, newNoise, oldNoise)
                newLine = '      %s : { gain : %.4f, readNoise : %.1f }\n'%(ampName, newGain, newNoise)
                file.write(newLine)
file.close()

```

```python

```

```python

```
