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

# Notebook for exploring the spotCatalog

Initially written 05 Feb 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
```

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import astropy.io.fits as pf

import eups
from lsst.daf.persistence import Butler
import lsst.afw.image as afwImage
import lsst.geom as geom
from lsst.daf.persistence import Butler
from lsst.ip.isr.isrTask import IsrTask, IsrTaskConfig
from lsst.ip.isr.isrFunctions import brighterFatterCorrection
from lsst.meas.algorithms import SourceDetectionTask
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
from lsst.geom import Point2I, Box2I
```

```python
DATA_DIR = '/project/shared/comCam/'
REPO_DIR = '/project/cslage/ComCam/20191230/'
OUTPUT_DIR = '/project/cslage/ComCam/20191230/'
DETECTOR = 4
raftName = 'R22'
```

```python
# Now set up the isrConfig and charConfig 
# The master bias, flat, and dark images have already been created and ingested.
butler = Butler(OUTPUT_DIR)

isrConfig = IsrTaskConfig()
isrConfig.doLinearize = False
isrConfig.doBias = True
isrConfig.doFlat = False
isrConfig.doDark = False
isrConfig.doFringe = False
isrConfig.doDefect = False
isrConfig.doAddDistortionModel = False
isrConfig.doWrite = False
isrConfig.doBrighterFatter = False
isrTask = IsrTask(config=isrConfig)

charConfig = CharacterizeImageConfig()
charConfig.installSimplePsf.fwhm = 1.0
charConfig.doMeasurePsf = False
charConfig.doApCorr = False
charConfig.doDeblend = False
charConfig.repair.doCosmicRay = True
charConfig.repair.doInterpolate = False   
charConfig.detection.background.binSize = 32
charConfig.detection.minPixels = 30
charTask = CharacterizeImageTask(config=charConfig)
```

```python
# First just try a single image with medium brightness
spot_visit=3019123000031
rawSpotDataRef = butler.dataRef('raw', detector=DETECTOR, visit=spot_visit)
postIsrSpot = isrTask.runDataRef(rawSpotDataRef).exposure
charResult = charTask.run(postIsrSpot)
spotCatalog = charResult.sourceCat

```

```python
# now let's look at several things to explore the spotCatalog
# type tells you what kind of object it is
print(type(spotCatalog))
```

```python
# After you know this, you can go to github and look at the code.  Sometimes this helps me, sometimes it doesn't.  
# For example, this code would be in https://github.com/lsst/afw This time, it doesn't help much
```

```python
# In the notebook environment, you can also use ? marks to explore variables:
spotCatalog?
```

```python
# dir tells you what it understands:
print(dir(spotCatalog))
```

```python
# You can then drill deeper, but it takes time.  For example:
print(dir(spotCatalog.schema))
print(dir(spotCatalog.table))
```

```python
# To see what types of measurements are available in the spotCatalog, try this:
for name in spotCatalog.schema.getOrderedNames():
    print(name)
```

```python

```
