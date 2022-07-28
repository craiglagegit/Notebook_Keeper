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

# Notebook for re-ingesting ComCam data.

Initially written 15 Nov 2019 by Craig Lage\
This ingests the images into my own repo, \
creates the master bias images, and ingests them.\
02-Dec-19 - I've added the master flats and master darks.
19-Dec-19  Adding code to add missing header keywords\
I'm setting up for BF correction, so this requires the following:\
obs_base: tickets/DM-18683\
obs_lsst: tickets/DM-18683\
cp_pipe: tickets/DM-18683\
ip_isr: tickets/DM-22659

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
DATA_DIR = '/project/shared/comCam/'
REPO_DIR = '/project/cslage/ComCam/20191218/'
OUTPUT_DIR = '/project/cslage/ComCam/20191218/'
DETECTOR = 4
raftName = 'R22'
```

```python
# First check the exposure times
filedir = DATA_DIR+'raw/20191218/'
files = glob.glob(filedir+'CC_C_20191218_00????/CC_C_20191218_00????_R22_S11.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)

for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
        
    phdr=hdulist[0].header
    filename = file.split('/')[6][14:]#phdr['FILENAME']
    exptime = phdr['EXPTIME']
    imgtype = phdr['IMGTYPE'] 
    print("%s\t%s\t%f"%(filename, imgtype, exptime))


```

```python jupyter={"outputs_hidden": true}
# Now we create the correction files to add the missing header keywords
# Correction file name format is comCam-CC_C_20191218_000001.yaml
filedir = DATA_DIR+'raw/20191218/'
files = glob.glob(filedir+'CC_C_20191218_00????/CC_C_20191218_00????_R22_S11.fits')
files = np.sort(files)
lines = 'LSST_NUM: ITL-3800C-206\n'
lines+= 'RAFTNAME: R22\n'
for file in files:
    date = file.split('/')[-1][5:13]
    seqnum =  file.split('/')[-1][14:20]
    #print(date, seqnum)
    #break
    filename = '/home/cslage/alternate_branches/obs_lsst/corrections/comCam-CC_C_%s_%s.yaml'%(date, seqnum)
    print(filename)
    correctionFile = open(filename, 'w')
    correctionFile.write(lines)
    correctionFile.close()
```

```python jupyter={"outputs_hidden": true}
# Now ingest the images
! mkdir -p {REPO_DIR}
! echo "lsst.obs.lsst.comCam.LsstComCamMapper" > {REPO_DIR+"_mapper"}
args = REPO_DIR + " " + DATA_DIR + "raw/20191218/*/*_R22_S11.fits" + " " + "--mode=link"
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

```python jupyter={"outputs_hidden": true}
# Skipping master flats - no flats in this dataset
# Now create a master flat
visits = []
starting_visit = 2019111300044
ending_visit = 2019111300062
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

args += " batch-type=None" + " -c isr.doCrosstalk=False" + " --clobber-config"
print(args)
! constructFlat.py {args}
```

```python jupyter={"outputs_hidden": true}
# Now ingest the master flat images
args = REPO_DIR + " " + RERUN_DIR + "/flat/*/*/*.fits" + " --validity 9999" + " --calib " + CALIB_DIR + " --mode=link"
! ingestCalibs.py {args} 
```

```python
# Now let's try running the ISR on a spot image. This is a medium exposure.
spot_visit = 3019121800081
butler = Butler(REPO_DIR)
rawSpot = butler.get('raw', detector=4, visit=spot_visit)
# this is the dataRef for running isr
rawSpotDataRef = butler.dataRef('raw', detector=4, visit=spot_visit)

isrConfig = IsrTask.ConfigClass()
isrConfig.doLinearize = False
isrConfig.doBias = True
isrConfig.doFlat = False
isrConfig.doDark = True
isrConfig.doFringe = False
isrConfig.doDefect = False
isrConfig.doAddDistortionModel = False
isrConfig.doWrite = False
isrTask = IsrTask(config=isrConfig)
# run the task and take the exposure
postIsrSpot = isrTask.runDataRef(rawSpotDataRef).exposure

charConfig = CharacterizeImageConfig()
charConfig.installSimplePsf.fwhm = 1.0
charConfig.doMeasurePsf = False
charConfig.doApCorr = False
charConfig.doDeblend = False
charConfig.repair.doCosmicRay = True
charConfig.repair.cosmicray.nCrPixelMax=100000
charConfig.repair.doInterpolate = False   
charConfig.detection.background.binSize = 32
charConfig.detection.minPixels = 10
charTask = CharacterizeImageTask(config=charConfig)
```

```python
# First just try a single image with medium brightness
rawSpotDataRef = butler.dataRef('raw', detector=DETECTOR, visit=spot_visit)
postIsrSpot = isrTask.runDataRef(rawSpotDataRef).exposure
charResult = charTask.run(postIsrSpot)
spotCatalog = charResult.sourceCat
# These two liines weed out some bad data
badSelect = spotCatalog['base_SdssShape_instFlux'] < 800000
spotCatalog = spotCatalog.subset(badSelect)

maxFlux = np.nanmax(spotCatalog['base_SdssShape_instFlux'])
print(maxFlux)
minFluxRatio = 0.80
select = spotCatalog['base_SdssShape_instFlux'] > maxFlux * minFluxRatio
numCat = len(spotCatalog)
plt.figure(figsize=(16,8))
plt.subplot(1,2,1, aspect=1)
plt.title('X/Y locations of detections - numCat = %d'%numCat)
color = spotCatalog['base_SdssShape_instFlux'] / maxFlux * 100.0
plt.scatter(spotCatalog['base_SdssCentroid_x'],spotCatalog['base_SdssCentroid_y'],c=color, cmap=plt.cm.jet, s=5)
plt.colorbar()
plt.xlim(-100,4200)
plt.ylim(-100,4200)

spotCatalog  = spotCatalog.subset(select)
numCat = len(spotCatalog)
plt.subplot(1,2,2, aspect=1)
plt.title('X/Y locations of detections - numCat = %d'%numCat)
color = spotCatalog['base_SdssShape_instFlux'] / maxFlux * 100.0
norm = plt.Normalize(vmin=0.0, vmax=100.0)
plt.scatter(spotCatalog['base_SdssCentroid_x'],spotCatalog['base_SdssCentroid_y'],c=color, norm=norm, cmap=plt.cm.jet, s=5)
plt.colorbar()
plt.xlim(-100,4200)
plt.ylim(-100,4200)
plt.savefig(OUTPUT_DIR+"plots/Spot_Intensities_19Dec19.pdf")
```

```python jupyter={"outputs_hidden": true, "source_hidden": true}
# This finds the bad spots
plt.figure(figsize=(16,16))
plt.subplot(1,1,1)
plt.title('Flux')
plt.hist(charResult.sourceCat['base_SdssShape_instFlux'])

plt.show()
```

```python
smallSelect = ((spotCatalog['base_SdssShape_instFlux'] > maxFlux * 0.98) & \
               (spotCatalog['base_SdssShape_instFlux'] < maxFlux * (1.0 - 1.0E-6)))
smallSpotCatalog = spotCatalog.subset(smallSelect)
maxSelect = spotCatalog['base_SdssShape_instFlux'] > maxFlux * (1.0 - 1.0E-6)
maxSpotCatalog = spotCatalog.subset(maxSelect)
print(len(maxSpotCatalog), len(smallSpotCatalog))

sep = np.sqrt(np.square(smallSpotCatalog['base_SdssCentroid_x'] - maxSpotCatalog['base_SdssCentroid_x'][0]) + \
              np.square(smallSpotCatalog['base_SdssCentroid_y'] - maxSpotCatalog['base_SdssCentroid_y'][0]))
print("Separation = %.2f pixels"%np.nanmin(sep))
```

```python
# Now just look at one of the spots
# A much cleaner spot profile than in the ComCam images
# This is certainly impacting the spot size algorithm
spotCatalog = charResult.sourceCat
xs = spotCatalog['base_SdssCentroid_x']
ys = spotCatalog['base_SdssCentroid_y']
spotNum = 277
deltaX = deltaY = 25
plotX = range(deltaX)
xCen = int(round(xs[spotNum]))
yCen = int(round(ys[spotNum]))
print(spotNum, xCen, yCen)
xMin = int(xCen - (deltaX-1)/2)
xMax = xMin + deltaX
yMin = int(yCen - (deltaY-1)/2)
yMax = yMin + deltaY
plt.figure(figsize=(16,8))    
plt.subplot(1,2,1)
plt.title("Spot %d X=%d; Y=%d"%(spotNum,xCen,yCen))
plt.imshow(postIsrSpot.image.array[yMin:yMax,xMin:xMax])
plt.colorbar()
plt.subplot(1,2,2)
plt.plot(plotX, postIsrSpot.image.array[yCen,xMin:xMax], label = "Y=%d"%yCen)
plt.plot(plotX, postIsrSpot.image.array[yMin:yMax,xCen], label = "X=%d"%xCen)
plt.xlabel("Pixels")
plt.ylabel("ADU")
plt.legend()
plt.savefig(OUTPUT_DIR+"plots/Spot_Profile_%d_%d_19Dec19.pdf"%(spot_visit, spotNum))
```

```python

```
