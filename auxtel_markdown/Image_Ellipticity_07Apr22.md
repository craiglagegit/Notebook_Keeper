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
from matplotlib.backends.backend_pdf import PdfPages
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

def PlotArrows(center, theta, xlabel, ylabel, x=[500.0,0], y=[0,500.0], flipX=False, flipY=False, color='lightgreen'):
    lshift = 100
    r = np.array(( (np.cos(theta), -np.sin(theta)),
               (np.sin(theta),  np.cos(theta)) ))
    x = np.array(x)
    if flipX: 
        x *= -1.0
    y = np.array(y)
    if flipY: 
        y *= -1.0
    rx = r.dot(x)
    ry = r.dot(y)
    tipx = center + rx
    tipy = center + ry
    plt.arrow(center[0],center[1], rx[0], rx[1], color=color, width = 20)
    plt.arrow(center[0],center[1], ry[0], ry[1], color=color, width = 20)
    plt.text(tipx[0]+lshift, tipx[1]+lshift, xlabel, \
        color=color, fontsize=12, weight='bold')
    plt.text(tipy[0]+lshift, tipy[1]+lshift, ylabel, \
        color=color, fontsize=12, weight='bold')
    return
    
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
filename = '/project/cslage/AuxTel/offsets/ellipses_267.fits'
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

```python tags=[]
charConfig = CharacterizeImageConfig()
charConfig.doMeasurePsf = False
charConfig.doApCorr = False
charConfig.doDeblend = False
charConfig.repair.doCosmicRay = True
charConfig.repair.doInterpolate = True   
charConfig.detection.minPixels = 500
charTask = CharacterizeImageTask(config=charConfig)
crop = Box2I(minimum=Point2I(250,250), maximum=Point2I(3750,3750))
cropped_exp = exp.subset(bbox=crop)
#charResult = charTask.run(cropped_exp)
charResult = charTask.run(exp)
sourceCatalog = charResult.sourceCat
print(f"expId:{expId}. Found {len(sourceCatalog)} sources,")

```

```python
plt.clf()

flipY = True
el = Angle(mData['ELSTART'] * u.deg)
az = Angle(mData['AZSTART'] * u.deg)
dec = Angle(mData['DECSTART'] * u.deg)
lat = AUXTEL_LOCATION.lat
rotpa = Angle(mData['ROTPA']*u.deg)

# This calculates the angle theta between (X,Y) and (Az,El)
sinTheta =  np.cos(lat) / np.cos(dec) * np.sin(az)
cosTheta = (np.sin(el) * np.sin(dec) - np.sin(lat)) / (np.cos(el) * np.cos(dec))
theta = Angle(np.arcsin(sinTheta))
# The following removes the ambiguity in arcsin(theta)
if cosTheta > 0:
    rotAzEl = rotpa - theta
else:    
    rotAzEl = rotpa - (Angle(180.0 * u.deg) - theta)

    # There are two posibilities for rotAzEl:
if flipY:
    flipX = False
else:
    rotAzEl += Angle(180.0 * u.deg)
    flipX = True

if rotAzEl.deg >= 360.0:
    rotAzEl -= Angle(360.0 * u.deg)
if rotAzEl.deg <= 0.0:
    rotAzEl += Angle(360.0 * u.deg)

plt.figure(figsize=(16,16))
arr = exp.image.array
arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')

# Get Az, El shift from wcs:
wcs = exp.getWcs()
origin = wcs.getSkyOrigin()
vi = exp.getInfo().getVisitInfo()

# Get Az, El shift from wcs:
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
PlotArrows(originPixel, 0.0, 'AZ-WCS', 'EL-WCS', x=azShift, y=altShift, color='lightgreen')

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
PlotArrows(originPixel, 0.0, 'N-WCS', 'E-WCS', x=decShift, y=raShift, color='orange')

# Plot N,E from rotpa
paOrigin = originPixel + np.array([500,500])
PlotArrows(paOrigin, rotpa, 'E-rotpa', 'N-rotpa', color='magenta')

# plot Az,El from rot AzEl
azelOrigin = originPixel + np.array([-500, -500])
PlotArrows(azelOrigin, rotAzEl, 'AZ-rotAzEl', 'EL-rotAzEl', flipY=True, color='cyan')

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
Nlength = 500.0
NcenterPA = (3000, 2000)
Nlabel = 650.0
NcenterFW = (2800, 3500)
FW_shift = 100
plt.arrow(NcenterPA[0],NcenterPA[1], Nlength*np.cos(median_phi), Nlength*np.sin(median_phi),\
    color='yellow', width = 20)
plt.text(NcenterPA[0]+Nlabel*np.cos(median_phi),NcenterPA[1]+Nlabel*np.sin(median_phi), 'Object_PA', \
    color='yellow', fontsize=12, weight='bold')
names = ['FWHM_x', 'FWHM_y', 'FWHM_az', 'FWHM_el', 'FWHM_a', 'FWHM_b']
for ii, FW in enumerate([FWHM_x, FWHM_y, FWHM_az, FWHM_el, FWHM_a, FWHM_b]):
    plt.text(NcenterFW[0], NcenterFW[1] - FW_shift * ii, names[ii] + f" = {FW:0.2f}", fontsize=24, \
             weight='bold', color='yellow')



plt.ylim(0,4000)
#plt.savefig(f'/project/cslage/AuxTel/offsets/Mount_Motions_{expId}.png')
plt.savefig(f'/project/cslage/AuxTel/offsets/GalSim_Test4_20Apr22.png')
```

```python
# Now check with multiple images
```

```python
flipY = True

expIdList = [2022031700454, 2022021600721, 2022031600645, 2022040600457, 2022040600654, 2022040600660, \
             2022040600723, 2022040601000, 2022021700321, 2022040500400, 2022021700311]

pdf = PdfPages(f"/project/cslage/AuxTel/offsets/Rotation_Check_{flipY}_20Apr22.pdf")
fig = plt.figure(figsize=(8,8))
for expId in expIdList:
    plt.title(f"Directions {expId}")
    exp = butler.get('quickLookExp', detector=0, exposure=expId)
    mData = exp.getMetadata()
    el = Angle(mData['ELSTART'] * u.deg)
    az = Angle(mData['AZSTART'] * u.deg)
    dec = Angle(mData['DECSTART'] * u.deg)
    lat = AUXTEL_LOCATION.lat
    rotpa = Angle(mData['ROTPA']*u.deg)

    # This calculates the angle theta between (X,Y) and (Az,El)
    sinTheta =  np.cos(lat) / np.cos(dec) * np.sin(az)
    cosTheta = (np.sin(el) * np.sin(dec) - np.sin(lat)) / (np.cos(el) * np.cos(dec))
    theta = Angle(np.arcsin(sinTheta))
    # The following removes the ambiguity in arcsin(theta)
    if cosTheta > 0:
        rotAzEl = rotpa - theta
    else:    
        rotAzEl = rotpa - (Angle(180.0 * u.deg) - theta)
        
        
    # There are two posibilities for rotAzEl:
    if flipY:
        flipX = False
    else:
        rotAzEl += Angle(180.0 * u.deg)
        flipX = True
        
    if rotAzEl.deg >= 360.0:
        rotAzEl -= Angle(360.0 * u.deg)
    if rotAzEl.deg <= 0.0:
        rotAzEl += Angle(360.0 * u.deg)
        
    print(f"{expId}, ROTPA={rotpa}, ROTAZEL={rotAzEl}")
    arr = exp.image.array
    arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
    plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', origin='lower', cmap='gray')

    # Get Az, El shift from wcs:
    wcs = exp.getWcs()
    origin = wcs.getSkyOrigin()
    vi = exp.getInfo().getVisitInfo()

    # Get Az, El shift from wcs:
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
    PlotArrows(originPixel, 0.0, 'AZ-WCS', 'EL-WCS', x=azShift, y=altShift, color='lightgreen')

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
    PlotArrows(originPixel, 0.0, 'N-WCS', 'E-WCS', x=decShift, y=raShift, color='orange')

    # Plot N,E from rotpa
    paOrigin = originPixel + np.array([500,500])
    PlotArrows(paOrigin, rotpa, 'E-rotpa', 'N-rotpa', color='magenta')

    # plot Az,El from rot AzEl
    azelOrigin = originPixel + np.array([-500, -500])
    #PlotArrows(azelOrigin, rotAzEl, 'AZ-rotAzEl', 'EL-rotAzEl', flipY=True, color='cyan')
    PlotArrows(azelOrigin, rotAzEl, 'AZ-rotAzEl', 'EL-rotAzEl', flipX=flipX, flipY=flipY, color='cyan')
    pdf.savefig(fig)  # saves the current figure into a pdf page
    plt.clf()
pdf.close()


```

```python

```
