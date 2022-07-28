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
and does assembly.  ISR is just CCD assembly, and bias subtraction.\
I applied the gains from last year's PTC measurements.\
The assembly of the raft is still manual at this point.

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
# Now create a master bias set
start=3020030300005
end=3020030300015
CALIB_DIR = REPO_DIR + "CALIB"
RERUN_DIR = REPO_DIR + "calib_construction"
! mkdir -p {CALIB_DIR}
args = REPO_DIR + " --calib " + CALIB_DIR + " --rerun " + RERUN_DIR + " --id visit=%d..%d"%(start,end) + \
" --batch-type=None" + " -c isr.doCrosstalk=False" + " --clobber-config"
! constructBias.py {args}
```

```python jupyter={"outputs_hidden": true}
# Now ingest the master bias images
args = REPO_DIR + " " + RERUN_DIR + "/bias/*/*.fits" + " --validity 9999" + " --calib " + CALIB_DIR + " --mode=link"
! ingestCalibs.py {args} 
```

```python jupyter={"outputs_hidden": true}
# Now create a master dark set
visits = []
starting_visit = 3020030300022
ending_visit = 3020030300026
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

args += " --batch-type=None" + " -c isr.doCrosstalk=False" + " --clobber-config"
print(args)
! constructDark.py {args}
```

```python jupyter={"outputs_hidden": true}
# Now ingest the master dark images
args = REPO_DIR + " " + RERUN_DIR + "/dark/*/*.fits" + " --validity 9999" + " --calib " + CALIB_DIR + " --mode=link"
! ingestCalibs.py {args} 
```

```python
# Set up the ISR task
# For now, this is just applying the bias and the gains
# For some reason, the darks are not working
isrConfig = IsrTask.ConfigClass()
isrConfig.doLinearize = False
isrConfig.doBias = True
isrConfig.doApplyGains = True
isrConfig.doFlat = False
isrConfig.doDark = False
isrConfig.doFringe = False
isrConfig.doDefect = False
isrConfig.doAddDistortionModel = False
isrConfig.doWrite = False
isrTask = IsrTask(config=isrConfig)
butler = Butler(REPO_DIR)
```

```python
# Apply gains from PTC data from last year
# and assemble all 9 CCDs
# Run all three pinhole images in two colors with and without gaps.
# The gaps are realistic estimates of the gaps between CCDs
GAIN_DIR = '/home/cslage/ComCam/20191113/'

aSize = 0.32 # Size of the pixel array
for vis in [34, 60, 64]:
    # 34 = Brian selfie, 60 = Rubin blouse, 64 = Rubin face
    if vis == 34:
        vmax = 8000
    else:
        vmax = 150000
    visit = 3020030300000 + vis
    for gap in ['Gap', 'NoGap']:
        for color in ['gray', 'viridis']:
            if gap == 'Gap':
                xGap = 0.0072 # Gap between imaging arrays
                yGap = 0.0032 # Gap between imaging arrays
            else:
                xGap = 0.0
                yGap = 0.0
            plt.figure(figsize=(16,16))
            xs = [0.0,aSize+xGap,2*(aSize+xGap),0.0,aSize+xGap,2*(aSize+xGap),0.0,aSize+xGap,2*(aSize+xGap)]
            ys = [0.0,0.0,0.0,aSize+yGap,aSize+yGap,aSize+yGap,2*(aSize+yGap),2*(aSize+yGap),2*(aSize+yGap)]
            for detector in range(9):
                gain_pickle_file = GAIN_DIR+'calibrations/ptc/ptcDataGainAndNoise-det%03d.pkl'%detector
                gain_file = open(gain_pickle_file, 'rb')
                gain_data = pkl.load(gain_file)
                raw = butler.get('raw', detector=detector, visit=visit)
                dataRef = butler.dataRef('raw', detector=detector, visit=visit)
                ccd = raw.getDetector()
                for amp in ccd:
                    amp = amp.rebuild()
                    amp.setGain(gain_data['gain'][amp.getName()][0])
                    amp.finish()
                    #print(detector, amp.getName(), amp.getGain())
                postIsr = isrTask.runDataRef(dataRef).exposure
                ax=plt.axes([xs[detector],ys[detector],aSize*(509.0/500.0),aSize],aspect=1)
                ax.set_xticks([])
                ax.set_yticks([])
                #ax.imshow(np.log10(postIsr.image.array[0:4072,0:4000]),vmin=2.5, vmax=vmax, cmap=color)
                ax.imshow(postIsr.image.array[0:4072,0:4000],vmin=10, vmax=vmax, cmap=color)
            plt.savefig(OUTPUT_DIR+"images/Image_Lin_%d_%s_%s.png"%(vis,gap,color))
```

```python

```
