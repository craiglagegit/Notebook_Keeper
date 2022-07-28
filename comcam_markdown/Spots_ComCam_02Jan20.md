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

# Notebook for running characterizing spots on ComCam data.

Initially written 02 jan 2020 by Craig Lage.\
Looking at spot profiles.\
I'm setting up for BF correction, so this requires the following:\
obs_base: tickets/DM-18683\
obs_lsst: tickets/DM-18683\
cp_pipe: tickets/DM-18683\
ip_isr: tickets/DM-22659

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
# This sets up the visits.  
starting_visit = 3019123000031
ending_visit =   3019123000079

visits = []
visit_1 = starting_visit
while visit_1 < ending_visit+1:
    visits.append(visit_1)
    visit_1 += 2 # Skipping the bias frames in between
print(visits)
print(len(visits))
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
plt.savefig(OUTPUT_DIR+"plots/Spot_Intensities_30Dec19.pdf")
```

```python
smallSelect = ((spotCatalog['base_SdssShape_instFlux'] > maxFlux * 0.90) & \
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
plt.figure(figsize=(16,16))
plt.subplot(2,2,1)
plt.title('X second moment')
plt.hist(charResult.sourceCat['base_SdssShape_xx'], range=(0.0,50.0))
plt.subplot(2,2,2)
plt.title('Y second moment')
plt.hist(charResult.sourceCat['base_SdssShape_yy'], range=(0.0,50.0))

plt.figure(figsize=(16,8))
plt.subplot(2,2,3)
plt.title('X second moment')
plt.hist(charResult.sourceCat['base_SdssShape_xx'], range=(2.0,7.0))
plt.subplot(2,2,4)
plt.title('Y second moment')
plt.hist(charResult.sourceCat['base_SdssShape_yy'], range=(2.0,7.0))

plt.show()
```

```python
# Now just look at one of the spots
xs = spotCatalog['base_SdssCentroid_x']
ys = spotCatalog['base_SdssCentroid_y']
spotNum = 53
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
plt.savefig(OUTPUT_DIR+"plots/Spot_Profile_%d_%d_30Dec19.pdf"%(spot_visit, spotNum))
```

```python jupyter={"outputs_hidden": true}
# Now we try running on all of the spot images

minSizeX = 2.0
maxSizeX = 7.0
minSizeY = 2.0
maxSizeY = 7.0
minX = minY = 500
maxX = maxY = 3500
minFluxRatio = 0.80

byamp_results = []
byamp_corrected_results = []
for i, spot_visit in enumerate(visits):
    numVisits = 1
    for do_bf_corr in [False]:
        isrConfig.doBrighterFatter = do_bf_corr
        isrTask = IsrTask(config=isrConfig)
        for j in range(numVisits):
            sub_visit = spot_visit + j
            #print("Getting exposure # %d"%sub_visit)
            rawSpotDataRef = butler.dataRef('raw', detector=DETECTOR, visit=sub_visit)

            postIsrSpot = isrTask.runDataRef(rawSpotDataRef).exposure
            charResult = charTask.run(postIsrSpot)
            spotCatalog = charResult.sourceCat
            maxFlux = np.nanmax(spotCatalog['base_SdssShape_instFlux'])
            select = spotCatalog['base_SdssShape_instFlux'] > maxFlux * minFluxRatio
            spotCatalog  = spotCatalog.subset(select)
            select = ((spotCatalog['base_SdssShape_xx'] >= minSizeX) & (spotCatalog['base_SdssShape_xx'] <= maxSizeX) & 
                    (spotCatalog['base_SdssShape_yy'] >= minSizeY) & (spotCatalog['base_SdssShape_yy'] <= maxSizeY) &
                    (spotCatalog['base_SdssCentroid_x'] >= minX) & (spotCatalog['base_SdssCentroid_x'] <= maxX) &
                    (spotCatalog['base_SdssCentroid_y'] >= minY) & (spotCatalog['base_SdssCentroid_y'] <= maxY))

            spotCatalog  = spotCatalog.subset(select)
            if j == 0:
                print("Correction = %s, Exposure # %d, %d spots"%(str(do_bf_corr),sub_visit, len(spotCatalog['base_SdssShape_instFlux'])))
                x2 = spotCatalog['base_SdssShape_xx']
                y2 = spotCatalog['base_SdssShape_yy']
                flux = spotCatalog['base_SdssShape_instFlux']
            else:
                print("Correction = %s, Exposure # %d, %d spots"%(str(do_bf_corr),sub_visit, len(spotCatalog['base_SdssShape_instFlux'])))
                x2 = np.concatenate((x2, spotCatalog['base_SdssShape_xx']))
                y2 = np.concatenate((y2, spotCatalog['base_SdssShape_yy']))
                flux = np.concatenate((flux, spotCatalog['base_SdssShape_instFlux']))

        numspots = len(flux)
        print("Correction = %s, Detected %d objects, Flux = %f, X2 = %.3f +/- %.3f, Y2 = %.3f +/- %.3f"%(str(do_bf_corr),numspots, \
                                np.nanmean(flux),np.nanmean(x2),np.nanstd(x2),np.nanmean(y2),np.nanstd(y2)))
        sys.stdout.flush()                                
        if do_bf_corr:
            byamp_corrected_results.append([numspots, np.nanmean(flux), np.nanstd(flux), np.nanmean(x2), np.nanstd(x2),
                                   np.nanmean(y2), np.nanstd(y2)])
        else:
            byamp_results.append([numspots, np.nanmean(flux), np.nanstd(flux), np.nanmean(x2), np.nanstd(x2),
                                   np.nanmean(y2), np.nanstd(y2)])
spots_pickle = {'results':byamp_results, 'corrected_results': byamp_corrected_results}
filename = OUTPUT_DIR+"/spots_results.pkl"
with open(filename, 'wb') as f:
    pkl.dump(spots_pickle, f)


```

```python
# Now plot the result
from scipy import stats
plotCorrection=False
syst_fraction = 0.25
min_slope_index = 2
max_slope_index = len(byamp_results) - 1
max_flux_index = len(byamp_results)
minSpot = 2.5
maxSpot = 4.5
#with open(filename, 'rb') as f:
#    spots_pickle= pkl.load(f)
#byamp_results = spots_pickle['results']
textDelta = (maxSpot - minSpot) / 10
# These next are in case not all fluxes produced good results

try:
    results = np.array([byamp_results[i] for i in range(max_flux_index)])
    max_slope_ind = max_slope_index
except:
    results = np.array(byamp_results)
    max_slope_ind = min(len(results) - 4, max_slope_index)
xerror = results[:,2]/np.sqrt(results[:,0])
xyerror = results[:,4] * (syst_fraction + (1 - syst_fraction) / np.sqrt(results[:,0]))
yyerror = results[:,6] * (syst_fraction + (1 - syst_fraction) / np.sqrt(results[:,0]))

if plotCorrection:
    #byamp_corrected_results = spots_pickle['corrected_results']
    try:
        corrected_results = np.array([byamp_corrected_results[i] for i in range(max_flux_index)])
        max_slope_ind_corr = max_slope_index
    except:
        corrected_results = np.array(byamp_corrected_results)           
        max_slope_ind_corr = min(len(corrected_results) - 4, max_slope_index)

    corrected_xerror = corrected_results[:,2]/np.sqrt(corrected_results[:,0])
    corrected_xyerror = corrected_results[:,4] * (syst_fraction + (1 - syst_fraction) / np.sqrt(corrected_results[:,0]))
    corrected_yyerror = corrected_results[:,6] * (syst_fraction + (1 - syst_fraction) / np.sqrt(corrected_results[:,0]))

plt.figure(figsize=(16,8))
plt.title("Brighter-Fatter - ComCam R22:S11", fontsize = 36)
# First plot the uncorrected data
plt.errorbar(results[:,1], results[:,3], xerr = xerror, 
             yerr = xyerror, color = 'green', lw = 2, label = 'X2', ls='', marker='x')
plt.errorbar(results[:,1], results[:,5], xerr = xerror, 
             yerr = yyerror, color = 'red', lw = 2, label = 'Y2', ls='',marker='x')
slope, intercept, r_value, p_value, std_err = stats.linregress(results[min_slope_index:max_slope_ind,1], results[min_slope_index:max_slope_ind,3])
xplot=np.linspace(-5000.0,3000000.0,100)
yplot = slope * xplot + intercept
#plt.plot(xplot, yplot, color='green', lw = 2, ls = '--')
tslope = slope * 100.0 * 200000.0
#plt.text(10000.0,maxSpot-textDelta,"X Slope = %.2f %% per 50K e-"%tslope, fontsize=24)

slope, intercept, r_value, p_value, std_err = stats.linregress(results[min_slope_index:max_slope_ind,1], results[min_slope_index:max_slope_ind,5])
xplot=np.linspace(-5000.0,3000000.0,100)
yplot = slope * xplot + intercept
#plt.plot(xplot, yplot, color='red', lw = 2, ls = '--')
tslope = slope * 100.0 * 200000.0
#plt.text(10000.0,maxSpot-2*textDelta,"Y Slope = %.2f %% per 50K e-"%tslope, fontsize=24)

if plotCorrection:
    # Now plot the corrected data
    plt.errorbar(corrected_results[:,1], corrected_results[:,3], xerr = corrected_xerror, 
                yerr = corrected_xyerror, color = 'cyan', lw = 2, ls='', marker='x', label = 'Corrected X2')
    plt.errorbar(corrected_results[:,1], corrected_results[:,5], xerr = corrected_xerror,
                yerr = corrected_yyerror, color = 'magenta', lw = 2, ls='', marker='x', label = 'Corrected Y2')
    slope, intercept, r_value, p_value, std_err = stats.linregress(corrected_results[min_slope_index:max_slope_ind_corr,1], corrected_results[min_slope_index:max_slope_ind_corr,3])
    xplot=np.linspace(-5000.0,3200000.0,100)
    yplot = slope * xplot + intercept
    plt.plot(xplot, yplot, color='cyan', lw = 2, ls = '--')
    tslope = slope * 100.0 * 200000.0
    plt.text(10000.0,maxSpot-3*textDelta,"Corrected X Slope = %.2f %% per 50K e-"%tslope, fontsize=24)

    slope, intercept, r_value, p_value, std_err = stats.linregress(corrected_results[min_slope_index:max_slope_ind_corr,1], corrected_results[min_slope_index:max_slope_ind_corr,5])
    xplot=np.linspace(-5000.0,3000000.0,100)
    yplot = slope * xplot + intercept
    plt.plot(xplot, yplot, color='magenta', lw = 2, ls = '--')
    tslope = slope * 100.0 * 200000.0
    plt.text(10000.0,maxSpot-4*textDelta,"Corrected Y Slope = %.2f %% per 50K e-"%tslope, fontsize=24)

plt.xlim(0.0,2000000.0)
plt.xticks([0,1000000,2000000])
plt.ylim(minSpot, maxSpot)
plt.xlabel('Spot Flux(electrons)',fontsize=24)
plt.ylabel('Second Moment (Pixels^2)',fontsize=24)
plt.legend(loc= 'upper right',fontsize = 18)
plt.savefig(OUTPUT_DIR+"plots/BF_Slopes_30Dec19.pdf")

```

```python

```

```python
# Looking at the brightest image
spot_visit=3019123000079
rawSpotDataRef = butler.dataRef('raw', detector=DETECTOR, visit=spot_visit)
postIsrSpot = isrTask.runDataRef(rawSpotDataRef).exposure
charResult = charTask.run(postIsrSpot)
spotCatalog = charResult.sourceCat

# These lines weed out some bad data
#medianFlux = np.nanmedian(spotCatalog['base_SdssShape_instFlux'])
#badSelect = (spotCatalog['base_SdssShape_instFlux'] > 0.25 * medianFlux) & \
#            (spotCatalog['base_SdssShape_instFlux'] < 4.0 * medianFlux)
#spotCatalog = spotCatalog.subset(badSelect)

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
#plt.savefig(OUTPUT_DIR+"plots/Spot_Intensities_30Dec19.pdf")
```

```python
# Now just look at one of the spots
# A much cleaner spot profile than in the ComCam images
# This is certainly impacting the spot size algorithm
#spotCatalog = charResult.sourceCat
xs = spotCatalog['base_SdssCentroid_x']
ys = spotCatalog['base_SdssCentroid_y']
spotNum = 53
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
plt.savefig(OUTPUT_DIR+"plots/Spot_Profile_%d_%d_30Dec19.pdf"%(spot_visit, spotNum))
```

```python

```
