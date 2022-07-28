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

## AuxTel AzEl offsets - 15-Apr-21

In this notebook, investigate az-el offsets from 11-Mar-21

```python
import sys, time, os, asyncio, glob

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pickle as pkl
import pandas as pd
import astropy.io.fits as pf
from astropy.time import Time, TimeDelta

from lsst.daf.butler import Butler as gen3Butler
from lsst.daf.persistence import Butler as gen2Butler
from lsst_efd_client import EfdClient
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python jupyter={"outputs_hidden": true}
# Gen3 butler
dayObs = 20210311
firstExpId = 2021031100346
lastExpId = 2021031100424
REPO_DIR = '/repo/main'
butler = gen3Butler(REPO_DIR, collections="LATISS/raw/all")

exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
    exposureList.append(record.id)
exposureList.sort()
myVisits = []
for exposure in exposureList:
    mData = butler.get('raw.metadata', detector=0, exposure=exposure)
    expTime = mData['EXPTIME']
    imgType = mData['IMGTYPE']
    obj = mData['OBJECT']
    filter = mData['FILTER']
    rotpa = mData['ROTPA']
    date_beg = mData['DATE-BEG']
    elstart = mData['ELSTART']
    azstart = mData['AZSTART']
    rastart = mData['RASTART']
    decstart = mData['DECSTART']
    dummy=0.0
    # Need to use DATE-BEG to get the right timestamp
    visit = (exposure, expTime, imgType, obj, filter, date_beg, rotpa, dummy, dummy, azstart, elstart, rastart, decstart)
    if (int(visit[0]) >= firstExpId) and (int(visit[0]) <= lastExpId):
        myVisits.append(visit)
        print(visit)
```

```python jupyter={"outputs_hidden": true}
mData = butler.get('raw.metadata', detector=0, exposure=2021031100422)
for key in mData.keys():
    print(key, mData[key])
```

```python
for visit in myVisits[0:2]:
    print(visit)
```

```python
# Get EFD client
client = EfdClient('summit_efd')
```

```python
# These are for finding the timestamps of the offset events
backUp = 5 # seconds before first image to get initial offset
start = Time(myVisits[0][5],scale='tai') - TimeDelta(backUp, format='sec')
end = Time(myVisits[-1][5],scale='tai') - TimeDelta(0, format='sec')
timestamp = f"time >= {start} AND time <= {end}"
```

```python
# Now get the offsets applied
offsets = await client.select_time_series("lsst.sal.ATPtg.command_offsetAzEl", ['*'],
                                          start, end)
```

```python
print(len(offsets))
```

```python
offsets.columns
```

```python
offsets.head(1)
```

```python
# Plot the first few to check the interleaving of the offsets with the exposures
# Blue are the times of setting the offsets, and red are the start of the exposure
startPlot = Time('2021-03-12T04:09:00') #this is UTC
endPlot = Time('2021-03-12T04:12:00')

fig = plt.figure()
plt.suptitle(f"Offsets - 11-Mar-21", fontsize = 18)
# Azimuth axis
plt.subplot(1,1,1)
ax1 = offsets['num'].plot(color='red', label='azimuth')
ax1.set_ylim(0,1.0)
for i in range(5):
    t1 = Time(offsets.index[i]).tai.isot
    ax1.axvline(t1, ymin=0.5, ymax=0.9, color="blue")
    t2 = Time(myVisits[i][5]).tai.isot
    ax1.axvline(t2, ymin=0.1, ymax=0.5, color="red")
ax1.set_xlim(startPlot.tai.isot,endPlot.tai.isot)
```

```python
# Plot a few more to check the interleaving of the offsets with the exposures
# Blue are the times of setting the offsets, and red are the start of the exposure
startPlot = Time('2021-03-12T04:23:00') #this is UTC
endPlot = Time('2021-03-12T04:32:00')

fig = plt.figure()
plt.suptitle(f"Offsets - 11-Mar-21", fontsize = 18)
# Azimuth axis
plt.subplot(1,1,1)
ax1 = offsets['num'].plot(color='red', label='azimuth')
ax1.set_ylim(0,1.0)
for i in range(77):
    t1 = Time(offsets.index[i]).tai.isot
    ax1.axvline(t1, ymin=0.5, ymax=0.9, color="blue")
    t2 = Time(myVisits[i][5]).tai.isot
    ax1.axvline(t2, ymin=0.1, ymax=0.5, color="red")
ax1.set_xlim(startPlot.tai.isot,endPlot.tai.isot)
```

```python jupyter={"outputs_hidden": true}
# Now append the applied offsets to the list of visits
# A few drop out because the offsets are not clear
backUp = 10
fullVisits = []
for i, visit in enumerate(myVisits):
    if i == 0:
        startTime = Time(myVisits[i][5],scale='tai',precision=0) - TimeDelta(backUp, format='sec')
        startTime = startTime.tai.isot
    else:
        startTime = Time(myVisits[i - 1][5],scale='tai',precision=0).tai.isot
    endTime = Time(myVisits[i][5],scale='tai',precision=0).tai.isot
    print(startTime, endTime)
    try:
        offset = offsets.loc[startTime:endTime].values
        if len(offset) == 1:
            newList = list(visit)
            newList[7] = offset[0][0]
            newList[8] = offset[0][1]
            fullVisits.append(newList)
        else:
            print("Not = 1", len(offset))
            continue
    except:
        print("Failed the try")
        continue
```

```python
len(fullVisits)
```

```python
# The last two values are the applied offsets in az and el in arcseconds
for fullVisit in fullVisits:
    print(fullVisit[0],fullVisit[7], fullVisit[8])
```

```python jupyter={"outputs_hidden": true}
# Get the raw quickLook data.  Only Gen2 works
REPO_DIR = '/project/shared/auxTel/rerun/quickLook'
gen2_butler = gen2Butler(REPO_DIR)
dayObs = '2021-03-11'
```

```python
charConfig = CharacterizeImageConfig()
charConfig.doMeasurePsf = False#True
charConfig.doApCorr = False
charConfig.doDeblend = False
charConfig.repair.doCosmicRay = True
charConfig.repair.doInterpolate = True   
charConfig.detection.minPixels = 500
charTask = CharacterizeImageTask(config=charConfig)
```

```python jupyter={"outputs_hidden": true}
# Try doing them all
charVisits = []
for fullVisit in fullVisits:
    expId = fullVisit[0]
    try:
        charVisit = {}
        charVisit['Visit'] = fullVisit
        exp = gen2_butler.get('quickLookExp', detector=0, expId=expId)
        charResult = charTask.run(exp)
        sourceCatalog = charResult.sourceCat

        maxFlux = np.nanmax(sourceCatalog['base_SdssShape_instFlux'])
        selectBrightestSource = sourceCatalog['base_SdssShape_instFlux'] > maxFlux * 0.99
        brightestSource = sourceCatalog.subset(selectBrightestSource)
        brightestCentroid = (brightestSource['base_SdssCentroid_x'][0], \
                             brightestSource['base_SdssCentroid_y'][0])
        brightCatalog = sourceCatalog.subset(sourceCatalog['base_SdssShape_instFlux'] > maxFlux * 0.001)
        print(f"expId:{expId}. Found {len(sourceCatalog)} sources, {len(brightCatalog)} bright sources")
        print(f"Brightest centroid at {brightestCentroid}")
        charVisit['exp'] = exp
        charVisit['brightestCentroid'] = brightestCentroid
        charVisit['brightCatalog'] = brightCatalog
        charVisits.append(charVisit)
    except:
        print(f"Skipping expId {expId}.")
        continue
```

```python jupyter={"outputs_hidden": true}
for charVisit in charVisits:
    print(charVisit['Visit'][6], charVisit['Visit'][7], charVisit['Visit'][8])
    print(charVisit['brightestCentroid'])
```

```python
outfile = open('/project/cslage/AuxTel/offsets/offsets_16apr21.pkl','wb')

pkl.dump(charVisits,outfile)
outfile.close()
```

```python
infile = open('/project/cslage/AuxTel/offsets/offsets_16apr21.pkl','rb')
charVisits = pkl.load(infile)
infile.close()
```

```python
for charVisit in charVisits:
    print(charVisit['Visit'][0], charVisit['Visit'][6], charVisit['Visit'][7], charVisit['Visit'][8])
    print(charVisit['brightestCentroid'])
```

```python
#%matplotlib inline
# Look at the data with matplotlib
def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar


#Plot these four
expIds = [2021031100406, 2021031100407, 2021031100410, 2021031100411]

plt.figure(figsize=(16,16))
plotCounter = 1
for charVisit in charVisits:
    expId = charVisit['Visit'][0]
    offset = (charVisit['Visit'][7],charVisit['Visit'][8])
    if expId not in expIds:
        continue
    plt.subplot(2,2,plotCounter)
    plotCounter += 1
    plt.title(f"Image - {expId}",fontsize=18)
    arr = charVisit['exp'].image.array
    arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
    img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')
    cat = charVisit['brightCatalog']
    plt.scatter(cat['base_SdssCentroid_x'],cat['base_SdssCentroid_y']\
                ,color='red', marker='x', label=f"offset={offset}")

    colorbar(img)
    plt.legend()
plt.tight_layout(h_pad=1)
#plt.savefig(f"/project/cslage/AuxTel/offsets/Offsets_{expIds[0]}_{expIds[1]}_{expIds[2]}_{expIds[3]}_16Apr21.pdf")
```

```python
dat = []
expIds = []
errors = []
angErrors = []
for charVisit in charVisits:
    expId = charVisit['Visit'][0]
    if expId == 2021031100424:
        break
    nextExpId = expId + 1
    for nextCharVisit in charVisits:
        thisExpId = nextCharVisit['Visit'][0]
        if thisExpId == nextExpId:
            break
    rotpa = charVisit['Visit'][6]
    off_x = charVisit['Visit'][7] - nextCharVisit['Visit'][7]
    off_y = charVisit['Visit'][8] - nextCharVisit['Visit'][8]
    off = np.sqrt(off_x**2 + off_y**2)  # This is the distance in arcseconds that we commanded
    t1 = np.arctan2(off_x, off_y)*180.0/np.pi # Angle relative to "up" in image
    
    shift_x = charVisit['brightestCentroid'][0] - nextCharVisit['brightestCentroid'][0]
    shift_y = charVisit['brightestCentroid'][1] - nextCharVisit['brightestCentroid'][1]
    shift_x *= exp.getWcs().getPixelScale().asArcseconds()
    shift_y *= exp.getWcs().getPixelScale().asArcseconds()
    shift = np.sqrt(shift_x**2 + shift_y**2)  # This is the distance in arcseconds that we measured
    t2 = np.arctan2(shift_x, shift_y)*180.0/np.pi # Angle relative to "up" in image
    dat.append([expId, rotpa, off, shift, t1, t2])
    error = off - shift
    if abs(error) < 10.0:
        expIds.append(expId - 2021031100000)
        errors.append(off - shift)
        angErrors.append(t1 - t2 + rotpa)
```

```python

```

```python
rmsError = np.sqrt(np.mean(np.array(errors) * np.array(errors)))

fig = plt.figure(figsize = (8,8))
plt.subplots_adjust(hspace = 0.7)
plt.subplot(2,1,1)
plt.title("Centroid Error \n Commanded shift - Measured shift", fontsize=18)
plt.scatter(expIds, errors, marker = 'x')
plt.text(360, -3.5, f"RMS error = {rmsError:.2f} arcseconds", fontsize=12)
plt.ylim(-5.0, 5.0)
plt.xlabel("Sequence number", fontsize=12)
plt.ylabel("Centroid error (arcseconds)", fontsize=12)
plt.subplot(2,1,2)
plt.title("Angular Error \n Commanded angle - Measured angle + ROTPA", fontsize=18)
plt.scatter(expIds, angErrors, marker = 'x')
plt.ylim(-12.0, 12.0)
plt.xlabel("Sequence number", fontsize=12)
plt.ylabel("Angular error (degrees)", fontsize=12)
plt.savefig(f"/project/cslage/AuxTel/offsets/Offset_errors_19Apr21.pdf")
```

```python
# Look at the data with matplotlib
def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar
```

```python

myExpId = 2021031100371
for charVisit in charVisits:
    expId = charVisit['Visit'][0]
    if expId == myExpId:
        break
nextExpId = myExpId + 1
for nextCharVisit in charVisits:
    thisExpId = nextCharVisit['Visit'][0]
    if thisExpId == nextExpId:
        break
cat = charVisit['brightCatalog']
nextCat = nextCharVisit['brightCatalog']
# These are the measured shifts between the two catalogs
shift_x = nextCharVisit['brightestCentroid'][0] - charVisit['brightestCentroid'][0]
shift_y = nextCharVisit['brightestCentroid'][1] - charVisit['brightestCentroid'][1] 
exp = charVisit['exp']
nextExp = nextCharVisit['exp']
rotpa = charVisit['Visit'][6]
# These are the commanded offsets in Az, El
off_az = nextCharVisit['Visit'][7] - charVisit['Visit'][7]
off_el = nextCharVisit['Visit'][8] - charVisit['Visit'][8]

# Now put off_az and off_el in pixels, and rotate them using rotpa
off_az /= exp.getWcs().getPixelScale().asArcseconds()
off_el /= exp.getWcs().getPixelScale().asArcseconds()

off = np.array([off_az, off_el])
theta = np.radians(rotpa)
c, s = np.cos(theta), np.sin(theta)
# This is the rotation matrix that puts the commanded offsets into the detector coordinates
R = np.array(((c, s), (-s, c))) 
rotated_off = R.dot(off)

# Now plot it all
plt.figure(figsize=(16,8))

plt.subplot(1,2,1)
plt.title(f"Image - {myExpId}",fontsize=18)
arr = exp.image.array
arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')
plt.scatter(cat['base_SdssCentroid_x'],cat['base_SdssCentroid_y']\
            ,color='red', marker='x', label="Measured")
plt.arrow(charVisit['brightestCentroid'][0],charVisit['brightestCentroid'][1], rotated_off[0], rotated_off[1],\
            color='green', width = 20, label='Commanded offset')
plt.arrow(charVisit['brightestCentroid'][0],charVisit['brightestCentroid'][1], shift_x, shift_y,\
            color='red', width=20, label='Measured offset')
plt.xlim(0,4000)
plt.ylim(4000,0)
colorbar(img)
plt.legend()

plt.subplot(1,2,2)
plt.title(f"Image - {nextExpId}",fontsize=18)
nextArr = nextExp.image.array
nextArr = np.clip(nextArr, 1, 100000) # This image has some negative values, and this removes them
img = plt.imshow(nextArr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')
plt.scatter(nextCat['base_SdssCentroid_x'],nextCat['base_SdssCentroid_y']\
            ,color='red', marker='x', label="Measured")
plt.scatter(cat['base_SdssCentroid_x'] + rotated_off[0],cat['base_SdssCentroid_y'] + rotated_off[1]\
            ,color='green', marker='+', s=200, label="Expected")
plt.xlim(0,4000)
plt.ylim(4000,0)
colorbar(img)
plt.legend()

plt.tight_layout(h_pad=1)
plt.savefig(f"/project/cslage/AuxTel/offsets/Offsets_Meas_vs_Expected_{myExpId}_19Apr21.pdf")
```

```python
# Looking at just one for debug.
```

```python
off = np.array([off_x, off_y])
theta = np.radians(rotpa)
c, s = np.cos(theta), np.sin(theta)
R = np.array(((c, -s), (s, c)))

rotated_off = R.dot(off)

```

```python
rotated_off
```

```python
charVisits = []
for fullVisit in fullVisits[0:1]:
    expId = fullVisit[0]
    #if expId == 2021031100348:
    #    continue
    print(expId)
    charVisit = {}
    charVisit['Visit'] = fullVisit
    exp = gen2_butler.get('quickLookExp', detector=0, expId=expId)
    charResult = charTask.run(exp)
    sourceCatalog = charResult.sourceCat

    maxFlux = np.nanmax(sourceCatalog['base_SdssShape_instFlux'])
    selectBrightestSource = sourceCatalog['base_SdssShape_instFlux'] > maxFlux * 0.99
    brightestSource = sourceCatalog.subset(selectBrightestSource)
    brightestCentroid = (brightestSource['base_SdssCentroid_x'][0], \
                         brightestSource['base_SdssCentroid_y'][0])
    brightCatalog = sourceCatalog.subset(sourceCatalog['base_SdssShape_instFlux'] > maxFlux * 0.001)
    print(f"Found {len(sourceCatalog)} sources, {len(brightCatalog)} bright sources")
    print(f"Brightest centroid at {brightestCentroid}")
    charVisit['exp'] = exp
    charVisit['brightestCentroid'] = brightestCentroid
    charVisit['brightCatalog'] = brightCatalog
    charVisits.append(charVisit)
```

```python
#%matplotlib inline
# Look at the data with matplotlib
def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

plt.figure(figsize=(8,8))
plt.subplot(1,1,1)
plt.title(f"Image - {expId}",fontsize=18)
arr = exp.image.array
arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')
cat = brightCatalog
plt.scatter(cat['base_SdssCentroid_x'],cat['base_SdssCentroid_y']\
            ,color='red', marker='x', label=f"offset={offset}")

colorbar(img)
plt.legend()
plt.tight_layout(h_pad=1)
#plt.savefig(f"/project/cslage/AuxTel/offsets/Offsets_{expIds[0]}_{expIds[1]}_{expIds[2]}_{expIds[3]}_16Apr21.pdf")
```
