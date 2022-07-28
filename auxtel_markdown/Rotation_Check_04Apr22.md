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

## AuxTel Image elongation due to azimuth oscillation

Craig Lage 17-Mar-22

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
from astropy.coordinates import AltAz, ICRS, EarthLocation, Angle, FK5, SkyCoord
import astropy.units as u
from lsst.obs.lsst.translators.latiss import AUXTEL_LOCATION
import lsst.geom as geom
from lsst.daf.butler import Butler
from lsst_efd_client import EfdClient
from lsst.geom import Box2I, Point2I, SpherePoint
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python
def RotatedMoments(Ixx, Iyy, Ixy, theta):
    # Rotates the moments about an angle theta.
    # Formulae are fron the Sextractor documentation
    # https://sextractor.readthedocs.io/en/latest/Position.html\
    # ?highlight=shape#basic-shape-parameters-a-b-theta
    # theta is assumed to be an astropy.coordinates.Angle instance
    c = np.cos(theta)
    s = np.sin(theta)
    IxxRot = c * c * Ixx + s * s * Iyy - 2.0 * c * s * Ixy
    IyyRot = s * s * Ixx + c * c * Iyy + 2.0 * c * s * Ixy
    IxyRot = c * s * (Ixx - Iyy) + (c * c - s * s) * Ixy
    return [IxxRot, IyyRot, IxyRot] 
```

```python
butler = Butler('/repo/main', collections="LATISS/runs/quickLook")
```

```python
# Now get the image data and the metadata

expId = 2022031700454
exp = butler.get('quickLookExp', detector=0, exposure=expId)
mData = exp.getMetadata()
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
time = date_beg.utc
print(expId, time)
```

```python
# Get the galsim generated data
filename = '/project/cslage/AuxTel/offsets/ellipses.fits'
hdulist = pf.open(filename, mode='readonly', do_not_scale_image_data=True)
hdr=hdulist[0].header
ellipse_data = hdulist[0].data
print(type(ellipse_data))
print(ellipse_data.shape)
```

```python
# Replace the real data with the GalSim ellipses
for i in range(4000):
    for j in range(4000):
        exp.image.array[i,j] = ellipse_data[i,j]
```

```python
charConfig = CharacterizeImageConfig()
charConfig.doMeasurePsf = False
charConfig.doApCorr = False
charConfig.doDeblend = False
charConfig.repair.doCosmicRay = True
charConfig.repair.doInterpolate = True   
charConfig.detection.minPixels = 500
charTask = CharacterizeImageTask(config=charConfig)
#crop = Box2I(minimum=Point2I(250,250), maximum=Point2I(3750,3750))
#cropped_exp = exp.subset(bbox=crop)
charResult = charTask.run(exp)
#charResult = charTask.run(cropped_exp)
sourceCatalog = charResult.sourceCat
print(f"expId:{expId}. Found {len(sourceCatalog)} sources,")

```

```python
# Used for arrow locations
Ncenter = (700, 900)
Nlength = 500.0
NcenterAzEl = (3200, 700)
NcenterPA = (3000, 1400)
Nlabel = 650.0
yShift = 150.0
NcenterFW = (3200, 3000)
FW_shift = 100
labelShift = 40

el = Angle(mData['ELSTART'] * u.deg)
az = Angle(mData['AZSTART'] * u.deg)
dec = Angle(mData['DECSTART'] * u.deg)

plt.figure(figsize=(16,16))
arr = exp.image.array
arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')

# Get Az, El shift from wcs:
wcs = exp.getWcs()
origin = wcs.getSkyOrigin()
vi = exp.getInfo().getVisitInfo()

skyLocation = SkyCoord(origin.getRa().asRadians(), origin.getDec().asRadians(), unit=u.rad)
altAz = AltAz(obstime=vi.date.toPython(), location=AUXTEL_LOCATION)
obsAltAz = skyLocation.transform_to(altAz)
shiftInArcsec = 50.0
altShifted = SkyCoord(obsAltAz.az, obsAltAz.alt + Angle(shiftInArcsec * u.arcsec), frame=altAz)
altShiftedSpherePoint = SpherePoint(altShifted.icrs.ra.deg*geom.degrees, altShifted.icrs.dec.deg*geom.degrees)
azShifted = SkyCoord(obsAltAz.az + Angle(shiftInArcsec / np.cos(obsAltAz.alt) * u.arcsec), obsAltAz.alt, frame=altAz)
azShiftedSpherePoint = SpherePoint(azShifted.icrs.ra.deg*geom.degrees, azShifted.icrs.dec.deg*geom.degrees)
originPixel = wcs.skyToPixel(origin)
altShiftedPixel = wcs.skyToPixel(altShiftedSpherePoint)
altShift = altShiftedPixel - originPixel
azShiftedPixel = wcs.skyToPixel(azShiftedSpherePoint)
azShift = azShiftedPixel - originPixel

# Now plot the Az El arrows as determined by wcs
plt.arrow(originPixel.x, originPixel.y, azShift.x, azShift.y,\
    color='orange', width = 20)
plt.text(azShiftedPixel.x + labelShift, azShiftedPixel.y + labelShift  , 'AZ - WCS', fontsize=12, weight='bold', color='orange')
plt.arrow(originPixel.x, originPixel.y, altShift.x, altShift.y,\
    color='orange', width = 20)    
plt.text(altShiftedPixel.x + labelShift  , altShiftedPixel.y + labelShift  , 'EL - WCS', fontsize=12, weight='bold', color='orange')


# Get N, E shift from wcs:
decShifted = SkyCoord(skyLocation.ra, skyLocation.dec + Angle(shiftInArcsec * u.arcsec))
decShiftedSpherePoint = SpherePoint(decShifted.ra.deg*geom.degrees, decShifted.dec.deg*geom.degrees)
raShifted = SkyCoord(skyLocation.ra + Angle((shiftInArcsec / np.cos(skyLocation.dec)) * u.arcsec), skyLocation.dec)
raShiftedSpherePoint = SpherePoint(raShifted.ra.deg*geom.degrees, raShifted.dec.deg*geom.degrees)
originPixel = wcs.skyToPixel(origin)
decShiftedPixel = wcs.skyToPixel(decShiftedSpherePoint)
decShift = decShiftedPixel - originPixel
raShiftedPixel = wcs.skyToPixel(raShiftedSpherePoint)
raShift = raShiftedPixel - originPixel

# Now plot the N E arrows as determined by wcs
plt.arrow(originPixel.x, originPixel.y, decShift.x, decShift.y,\
    color='magenta', width = 20)
plt.text(decShiftedPixel.x + labelShift, decShiftedPixel.y + labelShift  , 'N - WCS', fontsize=12, weight='bold', color='magenta')
plt.arrow(originPixel.x, originPixel.y, raShift.x, raShift.y,\
    color='magenta', width = 20)    
plt.text(raShiftedPixel.x + labelShift  , raShiftedPixel.y + labelShift  , 'E - WCS', fontsize=12, weight='bold', color='magenta')

# Now plot the N, E arrows as determined by ROTPA
rotpa = Angle(mData['ROTPA']*u.deg)
plt.arrow(Ncenter[0],Ncenter[1], -Nlength*np.sin(rotpa), Nlength*np.cos(rotpa),\
    color='lightgreen', width = 20)
plt.text(Ncenter[0]-Nlabel*np.sin(rotpa),Ncenter[1]+Nlabel*np.cos(rotpa), 'N', \
    color='lightgreen', fontsize=12, weight='bold')
plt.arrow(Ncenter[0],Ncenter[1], Nlength*np.cos(rotpa), Nlength*np.sin(rotpa),\
    color='lightgreen', width = 20)
plt.text(Ncenter[0]+Nlabel*np.cos(rotpa),Ncenter[1]+Nlabel*np.sin(rotpa), 'E', \
    color='lightgreen', fontsize=12, weight='bold')

# Now calculate and plot the Az, El arrows from my algorithm
# Now calculate and plot the Az, El arrows from my algorithm
lat = AUXTEL_LOCATION.lat
# This calculates the angle theta between (X,Y) and (Az,El)
sinTheta =  np.cos(lat) / np.cos(dec) * np.sin(az)
cosTheta = (np.sin(el) * np.sin(dec) - np.sin(lat)) / (np.cos(el) * np.cos(dec))
theta = Angle(np.arcsin(sinTheta))
# The following removes the ambiguity in arcsin(theta)
if cosTheta > 0:
    rotAzEl = rotpa - theta
else:    
    rotAzEl = rotpa - (Angle(180.0 * u.deg) - theta)

plt.arrow(NcenterAzEl[0],NcenterAzEl[1], Nlength*np.sin(rotAzEl),-Nlength*np.cos(rotAzEl),\
    color='cyan', width = 20)
plt.text(NcenterAzEl[0] + Nlabel*np.sin(rotAzEl),NcenterAzEl[1] - Nlabel*np.cos(rotAzEl), 'EL', \
    color='cyan', fontsize=12, weight='bold')
plt.arrow(NcenterAzEl[0],NcenterAzEl[1], Nlength*np.cos(rotAzEl), Nlength*np.sin(rotAzEl),\
    color='cyan', width = 20)
plt.text(NcenterAzEl[0]+Nlabel*np.cos(rotAzEl),NcenterAzEl[1]+Nlabel*np.sin(rotAzEl), 'AZ', \
    color='cyan', fontsize=12, weight='bold')

# Now determine and plot the ellipticities
Ixx = sourceCatalog.getIxx()
Iyy = sourceCatalog.getIyy()
Ixy = sourceCatalog.getIxy()
Ip = (Ixx + Iyy) / 2.0
Im = (Ixx - Iyy) / 2.0
phi = np.arctan2(Ixy, Im) / 2.0
median_phi = Angle(np.median(phi) * u.rad)
A2 = Ip + np.sqrt(Im**2 + Ixy**2)
B2 = Ip - np.sqrt(Im**2 + Ixy**2)
[Iaa, Iee, Iae] = RotatedMoments(Ixx, Iyy, Ixy, rotAzEl)
FWHM_x = 2.35 * np.median(np.sqrt(Ixx))
FWHM_y = 2.35 * np.median(np.sqrt(Iyy))
FWHM_az = 2.35 * np.median(np.sqrt(Iaa)) 
FWHM_el = 2.35 * np.median(np.sqrt(Iee)) 
FWHM_a = 2.35 * np.median(np.sqrt(A2)) 
FWHM_b = 2.35 * np.median(np.sqrt(B2)) 

print(theta.deg, rotpa.deg, rotAzEl.deg, median_phi.deg)
print(FWHM_x, FWHM_y, FWHM_az, FWHM_el, FWHM_a, FWHM_b)

plt.arrow(NcenterPA[0],NcenterPA[1], Nlength*np.cos(median_phi), Nlength*np.sin(median_phi),\
    color='yellow', width = 20)
plt.text(NcenterPA[0]+Nlabel*np.cos(median_phi),NcenterPA[1]+Nlabel*np.sin(median_phi), 'Object_PA', \
    color='yellow', fontsize=12, weight='bold')
names = ['FWHM_x', 'FWHM_y', 'FWHM_az', 'FWHM_el', 'FWHM_a', 'FWHM_b']
for ii, FW in enumerate([FWHM_x, FWHM_y, FWHM_az, FWHM_el, FWHM_a, FWHM_b]):
    plt.text(NcenterFW[0], NcenterFW[1] - FW_shift * ii, names[ii] + f" = {FW:0.2f}", fontsize=12, \
             weight='bold', color='yellow')



plt.ylim(0,4000)
#plt.savefig(f'/project/cslage/AuxTel/offsets/Mount_Motions_{expId}.png')
plt.savefig(f'/project/cslage/AuxTel/offsets/GalSim_Test2_05Apr22.png')
```

```python

```
