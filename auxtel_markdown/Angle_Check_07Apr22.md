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

expId = 2022021700311
exp = butler.get('quickLookExp', detector=0, exposure=expId)
mData = exp.getMetadata()
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
time = date_beg.utc
print(expId, time)
```

```python
# Used for arrow locations
Ncenter = (700, 900)
Nlength = 500.0
NcenterAzEl = (3200, 700)
NcenterPA = (3000, 2000)
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
sinTheta =  np.cos(AUXTEL_LOCATION.lat) / np.cos(dec) * np.sin(az)
cosTheta = (np.sin(el) * np.sin(dec) - np.sin(AUXTEL_LOCATION.lat)) / (np.cos(el) * np.cos(dec))

theta = Angle(np.arcsin(sinTheta))
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
"""
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
"""
print(az.deg, el.deg, rotpa.deg, rotAzEl.deg, theta.deg)
print(mData['AZSTART'], mData['ELSTART'], mData['ROTPA'])
print(mData['ROTCOORD'])
#print(FWHM_x, FWHM_y, FWHM_az, FWHM_el, FWHM_a, FWHM_b)
"""
plt.arrow(NcenterPA[0],NcenterPA[1], Nlength*np.cos(median_phi), Nlength*np.sin(median_phi),\
    color='yellow', width = 20)
plt.text(NcenterPA[0]+Nlabel*np.cos(median_phi),NcenterPA[1]+Nlabel*np.sin(median_phi), 'Object_PA', \
    color='yellow', fontsize=12, weight='bold')
names = ['FWHM_x', 'FWHM_y', 'FWHM_az', 'FWHM_el', 'FWHM_a', 'FWHM_b']
for ii, FW in enumerate([FWHM_x, FWHM_y, FWHM_az, FWHM_el, FWHM_a, FWHM_b]):
    plt.text(NcenterFW[0], NcenterFW[1] - FW_shift * ii, names[ii] + f" = {FW:0.2f}", fontsize=12, \
             weight='bold', color='yellow')

"""

plt.ylim(0,4000)
#plt.savefig(f'/project/cslage/AuxTel/offsets/Mount_Motions_{expId}.png')

```

```python
plt.close()
```

```python
expIdList = [2022031700454, 2022021600721, 2022031600645, 2022040600457, 2022040600654, 2022040600660, \
             2022040600723, 2022040601000, 2022021700321, 2022040500400, 2022021700311]
for expId in expIdList:
    exp = butler.get('quickLookExp', detector=0, exposure=expId)
    mData = exp.getMetadata()
    ni = Angle(90.0 * u.deg)
    on = Angle(180.0 * u.deg)
    el = Angle(mData['ELSTART'] * u.deg)
    az = Angle(mData['AZSTART'] * u.deg)
    dec = Angle(mData['DECSTART'] * u.deg)
    lat = AUXTEL_LOCATION.lat
    # Now calculate and plot the Az, El arrows from my algorithm
    sinTheta =  np.cos(lat) / np.cos(dec) * np.sin(az)
    cosTheta = (np.sin(el) * np.sin(dec) - np.sin(lat)) / (np.cos(el) * np.cos(dec))

    #sinTheta1 = np.sin(on-az) / np.sin(ni - dec) * np.sin(ni + lat)
    cosTheta1 = (np.cos(ni+lat) - np.cos(ni-el) * np.cos(ni+dec)) / (np.sin(ni-el) * np.sin(ni+dec))
    #print(expId, sinTheta, cosTheta)
    #print(expId, sinTheta, cosTheta1)
    print(cosTheta, cosTheta1)

```

```python

```
