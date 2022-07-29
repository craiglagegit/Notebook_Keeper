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

## AuxTel Astrometry - 15Jun-21

In this notebook, investigate astrometry for Merlin's solutions.

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
from astropy.coordinates import AltAz, ICRS, EarthLocation, Angle, FK5
import astropy.units as u
from astropy.table import Table
from astroquery.astrometry_net import AstrometryNet

from lsst.daf.butler import Butler as gen3Butler
from lsst.daf.persistence import Butler as gen2Butler
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python
# Set Cerro Pachon location
location = EarthLocation.from_geodetic(lon=-70.749417*u.deg,
                                       lat=-30.244639*u.deg,
                                       height=2663.0*u.m)

```

```python
# Now get the image data

expId = 2020031500163
# Get the raw quickLook data.  Only Gen2 works for now
REPO_DIR = '/project/shared/auxTel/rerun/quickLook'

gen2_butler = gen2Butler(REPO_DIR)
charConfig = CharacterizeImageConfig()
charConfig.doMeasurePsf = False#True
charConfig.doApCorr = False
charConfig.doDeblend = False
charConfig.repair.doCosmicRay = True
charConfig.repair.doInterpolate = True   
charConfig.detection.minPixels = 500
charTask = CharacterizeImageTask(config=charConfig)


exp = gen2_butler.get('quickLookExp', detector=0, expId=expId)
mData = exp.getMetadata()
charResult = charTask.run(exp)
sourceCatalog = charResult.sourceCat
maxFlux = np.nanmax(sourceCatalog['base_CircularApertureFlux_3_0_instFlux'])
selectBrightestSource = sourceCatalog['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.99
brightestSource = sourceCatalog.subset(selectBrightestSource)
brightestCentroid = (brightestSource['base_SdssCentroid_x'][0], \
                     brightestSource['base_SdssCentroid_y'][0])
brightCatalog = sourceCatalog.subset(sourceCatalog['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.001)
print(f"expId:{expId}. Found {len(sourceCatalog)} sources, {len(brightCatalog)} bright sources")
print(f"Brightest centroid at {brightestCentroid}")

```

```python
cat = brightCatalog
```

```python
Ncenter = (700, 900)
Nlength = 500.0
NcenterAzEl = (3200, 700)
Nlabel = 650.0
yShift = 150.0

plt.figure(figsize=(16,16))
arr = exp.image.array
arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')

plt.scatter(cat['base_SdssCentroid_x'],cat['base_SdssCentroid_y']\
            ,color='red', marker='x')
rotpa = Angle(mData['ROTPA']*u.deg)
plt.arrow(Ncenter[0],Ncenter[1], -Nlength*np.sin(rotpa), Nlength*np.cos(rotpa),\
    color='green', width = 20)
plt.text(Ncenter[0]-Nlabel*np.sin(rotpa),Ncenter[1]+Nlabel*np.cos(rotpa), 'N', \
    color='green', fontsize=12, weight='bold')
plt.arrow(Ncenter[0],Ncenter[1], Nlength*np.cos(rotpa), Nlength*np.sin(rotpa),\
    color='green', width = 20)
plt.text(Ncenter[0]+Nlabel*np.cos(rotpa),Ncenter[1]+Nlabel*np.sin(rotpa), 'E', \
    color='green', fontsize=12, weight='bold')

plt.ylim(0,4000)
plt.savefig(f'/project/cslage/AuxTel/offsets/{expId}.png')
```

```python
ast = AstrometryNet()
ast.api_key = 'xxawwhvleirxcswx'

```

```python
AstrometryNet.show_allowed_settings()
```

```python
sources = sourceCatalog.asAstropy()
sources.keep_columns(['base_SdssCentroid_x', 'base_SdssCentroid_y', 'base_CircularApertureFlux_3_0_instFlux'])
sources.sort('base_CircularApertureFlux_3_0_instFlux', reverse=True)
# Keep just the 17 bright sources
#sources.remove_rows(slice(17, 53, 1))
    
image_width = 4072
image_height = 4000
scale_units = 'arcsecperpix'
scale_type='ev' # ev means submit estimate and % error
scale_est = 0.095
scale_err = 2.0
center_ra = mData['RA']
center_dec = mData['DEC']
radius = 0.5
wcs_header = ast.solve_from_source_list(sources['base_SdssCentroid_x'], sources['base_SdssCentroid_y'],
                                        image_width, image_height, scale_units=scale_units,
                                        scale_type=scale_type, scale_est=scale_est, scale_err=scale_err,
                                        center_ra=center_ra, center_dec=center_dec, radius=radius,
                                        solve_timeout=240)
```

```python

```

```python
# Now get the image data

expId = 2021061000202
# Get the raw quickLook data.  Only Gen2 works for now

exp = gen2_butler.get('quickLookExp', detector=0, expId=expId)
mData = exp.getMetadata()
charResult = charTask.run(exp)
sourceCatalog = charResult.sourceCat
maxFlux = np.nanmax(sourceCatalog['base_CircularApertureFlux_3_0_instFlux'])
selectBrightestSource = sourceCatalog['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.99
brightestSource = sourceCatalog.subset(selectBrightestSource)
brightestCentroid = (brightestSource['base_SdssCentroid_x'][0], \
                     brightestSource['base_SdssCentroid_y'][0])
brightCatalog = sourceCatalog.subset(sourceCatalog['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.001)
print(f"expId:{expId}. Found {len(sourceCatalog)} sources, {len(brightCatalog)} bright sources")
print(f"Brightest centroid at {brightestCentroid}")

```

```python
cat = brightCatalog
```

```python
Ncenter = (700, 900)
Nlength = 500.0
NcenterAzEl = (3200, 700)
Nlabel = 650.0
yShift = 150.0

plt.figure(figsize=(16,16))
arr = exp.image.array
arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')

plt.scatter(cat['base_SdssCentroid_x'],cat['base_SdssCentroid_y']\
            ,color='red', marker='x')
rotpa = Angle(mData['ROTPA']*u.deg)
plt.arrow(Ncenter[0],Ncenter[1], -Nlength*np.sin(rotpa), Nlength*np.cos(rotpa),\
    color='green', width = 20)
plt.text(Ncenter[0]-Nlabel*np.sin(rotpa),Ncenter[1]+Nlabel*np.cos(rotpa), 'N', \
    color='green', fontsize=12, weight='bold')
plt.arrow(Ncenter[0],Ncenter[1], Nlength*np.cos(rotpa), Nlength*np.sin(rotpa),\
    color='green', width = 20)
plt.text(Ncenter[0]+Nlabel*np.cos(rotpa),Ncenter[1]+Nlabel*np.sin(rotpa), 'E', \
    color='green', fontsize=12, weight='bold')

plt.ylim(0,4000)
plt.savefig(f'/project/cslage/AuxTel/offsets/{expId}.png')
```

```python
sources = sourceCatalog.asAstropy()
sources.keep_columns(['base_SdssCentroid_x', 'base_SdssCentroid_y', 'base_CircularApertureFlux_3_0_instFlux'])
sources.sort('base_CircularApertureFlux_3_0_instFlux', reverse=True)
# Keep just the 17 bright sources
#sources.remove_rows(slice(17, 53, 1))
    
image_width = 4072
image_height = 4000
scale_units = 'arcsecperpix'
scale_type='ev' # ev means submit estimate and % error
scale_est = 0.095
scale_err = 2.0
center_ra = mData['RA']
center_dec = mData['DEC']
radius = 0.5
wcs_header = ast.solve_from_source_list(sources['base_SdssCentroid_x'], sources['base_SdssCentroid_y'],
                                        image_width, image_height, scale_units=scale_units,
                                        scale_type=scale_type, scale_est=scale_est, scale_err=scale_err,
                                        center_ra=center_ra, center_dec=center_dec, radius=radius,
                                        solve_timeout=240)
```

```python
# Brought the position data back in with:
# curl http://nova.astrometry.net/corr_file/5339307 --output corr_2021061000202.fits
# Need to automate this
```

```python
corr_file = pf.open(f"/project/cslage/AuxTel/offsets/corr_{expId}.fits")
corr_data = corr_file[1].data
```

```python
for dat in corr_data:
    print(dat[4], dat[5])
```

```python
Ncenter = (700, 900)
Nlength = 500.0
NcenterAzEl = (3200, 700)
Nlabel = 650.0
yShift = 150.0

plt.figure(figsize=(16,16))
arr = exp.image.array
arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')


plt.scatter(cat['base_SdssCentroid_x'],cat['base_SdssCentroid_y']\
            ,facecolors='none', edgecolors='r', s=200)

for dat in corr_data:
    plt.scatter(dat[4],dat[5]\
            ,facecolors='none', edgecolors='g', s=200, lw=2)

rotpa = Angle(mData['ROTPA']*u.deg)
plt.arrow(Ncenter[0],Ncenter[1], -Nlength*np.sin(rotpa), Nlength*np.cos(rotpa),\
    color='green', width = 20)
plt.text(Ncenter[0]-Nlabel*np.sin(rotpa),Ncenter[1]+Nlabel*np.cos(rotpa), 'N', \
    color='green', fontsize=12, weight='bold')
plt.arrow(Ncenter[0],Ncenter[1], Nlength*np.cos(rotpa), Nlength*np.sin(rotpa),\
    color='green', width = 20)
plt.text(Ncenter[0]+Nlabel*np.cos(rotpa),Ncenter[1]+Nlabel*np.sin(rotpa), 'E', \
    color='green', fontsize=12, weight='bold')

plt.ylim(0,4000)
plt.savefig(f'/project/cslage/AuxTel/offsets/{expId}_astropy.png')
```

```python

```
