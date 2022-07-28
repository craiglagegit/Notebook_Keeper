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

# Notebook for ingesting ComCam pinhole images.

Initially written 04 Mar 2020 by Craig Lage\
This ingests the images into my own repo, \
and does a test assembly

```python
! eups list -s | grep lsst_distrib
```

```python
import eups
import sys, os, glob
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf

from lsst.daf.persistence import Butler
from lsst.ip.isr.isrTask import IsrTask
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python
DATA_DIR = '/lsstdata/offline/teststand/comcam/CCS/storage/'
REPO_DIR = '/project/cslage/ComCam/20200303/'
OUTPUT_DIR = '/project/cslage/ComCam/20200303/'
! mkdir -p {'/project/cslage/ComCam/20200303/images'}
DETECTOR = 4
raftName = 'R22'
```

```python jupyter={"outputs_hidden": true}
# First check the exposure times
filedir = DATA_DIR+'20200303/'
files = glob.glob(filedir+'CC_C_20200303_00????/CC_C_20200303_00????_R22_S11.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)

for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
        
    phdr=hdulist[0].header
    filename = file.split('/')[8][14:]#phdr['FILENAME']
    exptime = phdr['EXPTIME']
    imgtype = phdr['IMGTYPE'] 
    print("%s\t%s\t%f"%(filename, imgtype, exptime))

```

```python jupyter={"outputs_hidden": true}
# Now ingest the images
! mkdir -p {REPO_DIR}
! echo "lsst.obs.lsst.comCam.LsstComCamMapper" > {REPO_DIR+"_mapper"}
args = REPO_DIR + " " + DATA_DIR + "20200303/*/*_R22_S??.fits" + " " + "--mode=link"
! ingestImages.py {args}
```

```python jupyter={"outputs_hidden": true}
# Now create a master bias
start=3019121800005
end=3019121800010
CALIB_DIR = REPO_DIR + "CALIB"
RERUN_DIR = REPO_DIR + "calib_construction"
! mkdir -p {CALIB_DIR}
args = REPO_DIR + " --calib " + CALIB_DIR + " --rerun " + RERUN_DIR + " --id visit=%d..%d"%(start,end) + \
" --batch-type=None" + " -c isr.doCrosstalk=False" + " --clobber-config"
! constructBias.py {args}
```

```python
# Now ingest the master bias image
args = REPO_DIR + " " + RERUN_DIR + "/bias/*/*.fits" + " --validity 9999" + " --calib " + CALIB_DIR + " --mode=link"
! ingestCalibs.py {args} 
```

```python jupyter={"outputs_hidden": true}
# Now create a master dark
# It failed with the default number of cosmic ray pixels = 10000
# Increased this to 100,000 and then it ran.
visits = []
starting_visit = 3019121800022
ending_visit = 3019121800022
visit = starting_visit
while visit < ending_visit + 1:
    visits.append(visit)
    visit += 2
print(len(visits))
CALIB_DIR = REPO_DIR + "CALIB"
RERUN_DIR = REPO_DIR + "calib_construction"

args = REPO_DIR + " --calib " + CALIB_DIR + " --rerun " + RERUN_DIR + " --id visit="
for visit in visits:
    if visit != starting_visit:
        args += "^"
    args += str(visit)

args += " --batch-type=None" + " -c isr.doCrosstalk=False repair.cosmicray.nCrPixelMax=100000" + " --clobber-config"
print(args)
! constructDark.py {args}
```

```python
# Now ingest the master dark image
args = REPO_DIR + " " + RERUN_DIR + "/dark/*/*.fits" + " --validity 9999" + " --calib " + CALIB_DIR + " --mode=link"
! ingestCalibs.py {args} 
```

```python
# Set up the ISR task
# For now, no ISR, just assembly
isrConfig = IsrTask.ConfigClass()
isrConfig.doLinearize = False
isrConfig.doBias = False
isrConfig.doFlat = False
isrConfig.doDark = False
isrConfig.doFringe = False
isrConfig.doDefect = False
isrConfig.doAddDistortionModel = False
isrConfig.doWrite = False

```

```python
# First just look at the center CCD
spot_visit = 3020030300054
butler = Butler(REPO_DIR)
rawSpotDataRef = butler.dataRef('raw', detector=4, visit=spot_visit)
isrTask = IsrTask(config=isrConfig)
postIsrSpot = isrTask.runDataRef(rawSpotDataRef).exposure
plt.figure(figsize=(16,16))    
plt.subplot(1,1,1)
plt.title("Raw Image")
plt.imshow(np.log10(postIsrSpot.image.array[0:4000,0:4000]),vmin=2.5, vmax=5.0)
plt.colorbar()
plt.savefig(OUTPUT_DIR+"images/Image_054_S11_Log.png")
```

```python
# Now assemble all 9
plt.figure(figsize=(16,16))
xs = [0.0,0.333,0.667,0.0,0.333,0.667,0.0,0.333,0.667]
ys = [0.0,0.0,0.0,0.333,0.333,0.333,0.667,0.667,0.667]

for ccd in range(9):
    rawSpotDataRef = butler.dataRef('raw', detector=ccd, visit=spot_visit)
    postIsrSpot = isrTask.runDataRef(rawSpotDataRef).exposure
    ax=plt.axes([xs[ccd],ys[ccd],0.333,0.333],aspect=1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.imshow(np.log10(postIsrSpot.image.array[0:4000,0:4000]),vmin=2.5, vmax=5.0)

plt.savefig(OUTPUT_DIR+"images/Image_054_All_Log.png")
```

```python

```
