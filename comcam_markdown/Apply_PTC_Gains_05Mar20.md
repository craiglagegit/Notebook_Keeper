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

# Notebook for taking gains extracted from PTC and applying them to a new repo.

Initially written 05 Mar 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
```

```python
import eups
from lsst.daf.persistence import Butler
from lsst.cp.pipe.ptc import MeasurePhotonTransferCurveTask
import sys, os, glob
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
```

```python
REPO_DIR = '/project/shared/comCam/'
GAIN_DIR = '/home/cslage/ComCam/20191113/'
OUTPUT_DIR = '/project/cslage/ComCam/20200303/'
raftName = 'R22'
```

```python
butler = Butler(OUTPUT_DIR)
visit = 3020030300034
dataRef = butler.dataRef('brighterFatterGain', detector=4, visit=visit)
```

```python
raw = butler.get('raw', detector=4, visit=visit)
```

```python
for key in dataRef.dataId.keys():
    print(key)

```

```python
print(dir(dataRef.put()))
```

```python
ccd = raw.getDetector()
print(dir(ccd))
print(dir(raw))
for amp in ccd:
    print(dir(amp))
    sys.exit()
```

```python
#exposure = dataRef.exposure
ccd = raw.getDetector()
builder = ccd.rebuild()
for amp in builder:
    #amp = amp.rebuild()
    amp.setGain(1.33)
    #amp.finish()
    print(amp.getName(),amp.getGain())
raw.setDetector(builder.finish())
ccd = raw.getDetector()
for amp in ccd:
    print(amp.getName(),amp.getGain())

```

```python
print(dir(dataRef))
```

```python
test = dataRef.get
print(dir(test))
```

```python
test = dataRef.put
```

```python
for detector in range(9):
    gain_pickle_file = GAIN_DIR+'calibrations/ptc/ptcDataGainAndNoise-det%03d.pkl'%detector
    gain_file = open(gain_pickle_file, 'rb')
    gain_data = pkl.load(gain_file)
    amps = gain_data['gain'].keys()
    gains = []
    names = []
    for amp in amps:
        print(detector, amp, gain_data['gain'][amp][0])


```

```python

```
