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

## AuxTel AzEl offset coordinates - 22-Apr-21

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
from astropy.coordinates import SkyCoord, AltAz, ICRS, GCRS, CIRS, EarthLocation, Angle, FK5, SkyOffsetFrame
import astropy.units as u

from lsst.daf.butler import Butler as gen3Butler
from lsst.daf.persistence import Butler as gen2Butler
from lsst_efd_client import EfdClient
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python
# Gen3 butler
dayObs = 20210311
expId = 2021031100422
REPO_DIR = '/repo/main'
butler = gen3Butler(REPO_DIR, collections="LATISS/raw/all")
mData = butler.get('raw.metadata', detector=0, exposure=expId)
for key in mData.keys():
    print(key, mData[key])

```

```python
# Set Cerro Pachon location and observation time
location = EarthLocation.from_geodetic(lon=mData['OBS-LONG']*u.deg,
                                       lat=mData['OBS-LAT']*u.deg,
                                       height=mData['OBS-ELEV']*u.m)

utcoffset = -3*u.hour  
time = Time(mData['DATE-BEG']) + utcoffset
time.format = 'iso'
```

```python
print(location.lat.deg, time)
```

```python
objectCoords = SkyCoord.from_name(mData['OBJECT'])
```

```python
objectCoords.dec.deg
```

```python
objectAltaz = objectCoords.transform_to(AltAz(obstime=time,location=location, \
                                              pressure=mData['PRESSURE']*u.torr, temperature=mData['AIRTEMP']*u.deg_C, \
                                             relative_humidity=mData['HUMIDITY']*u.percent))
```

```python
objectAltaz
```

```python
sinThetaObject =  np.cos(location.lat.rad) / np.cos(objectCoords.dec.rad) * np.sin(objectAltaz.az.rad)
print(sinThetaObject,np.cos(location.lat.rad),np.cos(objectCoords.dec.rad),np.sin(objectAltaz.az.rad))
```

```python
pointingCoords_1 = SkyCoord(ra=mData['RA']*u.degree, dec=mData['DEC']*u.degree)
```

```python
pointingCoords_1.dec.deg
```

```python
pointingAltaz_1 = pointingCoords_1.transform_to(AltAz(obstime=time,location=location, \
                                              pressure=mData['PRESSURE']*u.torr, temperature=mData['AIRTEMP']*u.deg_C, \
                                             relative_humidity=mData['HUMIDITY']*u.percent))
```

```python
pointingAltaz_1
```

```python
sinThetaPointing_1 =  np.cos(location.lat.rad) / np.cos(pointingCoords_1.dec.rad) * np.sin(pointingAltaz_1.az.rad)
print(sinThetaPointing_1,np.cos(location.lat.rad),np.cos(pointingCoords_1.dec.rad),np.sin(pointingAltaz_1.az.rad))
```

```python
pointingCoords_2 = SkyCoord(ra=mData['RASTART']*u.degree, dec=mData['DECSTART']*u.degree)
```

```python
pointingAltaz_2 = pointingCoords_2.transform_to(AltAz(obstime=time,location=location, \
                                              pressure=mData['PRESSURE']*u.torr, temperature=mData['AIRTEMP']*u.deg_C, \
                                             relative_humidity=mData['HUMIDITY']*u.percent))
```

```python
pointingAltaz_2
```

```python
sinThetaPointing_2 =  np.cos(location.lat.rad) / np.cos(pointingCoords_2.dec.rad) * np.sin(pointingAltaz_2.az.rad)
print(sinThetaPointing_2,np.cos(location.lat.rad),np.cos(pointingCoords_2.dec.rad),np.sin(pointingAltaz_2.az.rad))
```

```python
pointingCoords_3 = SkyCoord(ra=mData['RAEND']*u.degree, dec=mData['DECEND']*u.degree)
```

```python
pointingAltaz_3 = pointingCoords_3.transform_to(AltAz(obstime=time,location=location, \
                                              pressure=mData['PRESSURE']*u.torr, temperature=mData['AIRTEMP']*u.deg_C, \
                                             relative_humidity=mData['HUMIDITY']*u.percent))
```

```python
pointingAltaz_3
```

```python
sinThetaPointing_3 =  np.cos(location.lat.rad) / np.cos(pointingCoords_3.dec.rad) * np.sin(pointingAltaz_3.az.rad)
print(sinThetaPointing_3,np.cos(location.lat.rad),np.cos(pointingCoords_3.dec.rad),np.sin(pointingAltaz_3.az.rad))
```

```python
pointingCoords_header = SkyCoord(ra=mData['RASTART']*u.degree, dec=mData['DECSTART']*u.degree)
```

```python
altazCoords_header = AltAz(alt=mData['ELSTART']*u.deg, az=mData['AZSTART']*u.deg, obstime=time, location=location)
```

```python
sinThetaPointing_header =  np.cos(location.lat.rad) / np.cos(pointingCoords_header.dec.rad) * np.sin(altazCoords_header.az.rad)
print(sinThetaPointing_header,np.cos(location.lat.rad),np.cos(pointingCoords_header.dec.rad),np.sin(altazCoords_header.az.rad))
```

```python
pointingCoords_header_2 = SkyCoord(ra=mData['RAEND']*u.degree, dec=mData['DECEND']*u.degree)
```

```python
altazCoords_header_2 = AltAz(alt=mData['ELEND']*u.deg, az=mData['AZEND']*u.deg, obstime=time, location=location)
```

```python
sinThetaPointing_header_2 =  np.cos(location.lat.rad) / np.cos(pointingCoords_header_2.dec.rad) * np.sin(altazCoords_header_2.az.rad)
print(sinThetaPointing_header_2,np.cos(location.lat.rad),np.cos(pointingCoords_header_2.dec.rad),np.sin(altazCoords_header_2.az.rad))
```

```python
# These are as extracted by astrometry.net in the notebook AzEl_Offset_Astrometry_26Apr21.ipynb
astroNetCoords = SkyCoord(ra=240.141563967*u.degree, dec=-89.2964243801*u.degree)
```

```python
astroNetCoords.dec.deg
```

```python
astroNetAltaz = astroNetCoords.transform_to(AltAz(obstime=time,location=location, \
                                              pressure=mData['PRESSURE']*u.torr, temperature=mData['AIRTEMP']*u.deg_C, \
                                             relative_humidity=mData['HUMIDITY']*u.percent))
```

```python
astroNetAltaz
```

```python
sinThetaAstroNet =  np.cos(location.lat.rad) / np.cos(astroNetCoords.dec.rad) * np.sin(astroNetAltaz.az.rad)
print(sinThetaAstroNet,np.cos(location.lat.rad),np.cos(astroNetCoords.dec.rad),np.sin(astroNetAltaz.az.rad))
```

```python
# These are from Stellarium
stellCoords = SkyCoord(ra=250.830*u.degree, dec=-89.3566*u.degree)
```

```python
stellCoords.dec.deg
```

```python
stellAltaz = stellCoords.transform_to(AltAz(obstime=time,location=location, \
                                              pressure=mData['PRESSURE']*u.torr, temperature=mData['AIRTEMP']*u.deg_C, \
                                             relative_humidity=mData['HUMIDITY']*u.percent))
```

```python
stellAltaz
```

```python
stellAltaz_2 = SkyCoord(AltAz(az=179.3781*u.deg, alt=30.600*u.deg))
```

```python
sinThetaStell =  np.cos(location.lat.rad) / np.cos(stellCoords.dec.rad) * np.sin(stellAltaz_2.az.rad)
print(sinThetaStell,np.cos(location.lat.rad),np.cos(stellCoords.dec.rad),np.sin(stellAltaz_2.az.rad))
print(90.0 - np.arcsin(sinThetaStell)*180.0/np.pi)
```

```python
print(f"Header discrepancy expId {expId}")
print()
print("Source\t\t\tRA\t\tDec\t\tAlt\t\tAz\t\tsin(theta)")
print(f"Object+Conversion\t{objectCoords.ra.deg:.6f}\t{objectCoords.dec.deg:.6f}"  
 + f"\t{objectAltaz.alt.deg:.6f}\t{objectAltaz.az.deg:.6f}\t{sinThetaObject:.4f}")
print(f"RA/DEC+Conversion\t{pointingCoords_1.ra.deg:.6f}\t{pointingCoords_1.dec.deg:.6f}"  
 + f"\t{pointingAltaz_1.alt.deg:.6f}\t{pointingAltaz_1.az.deg:.6f}\t{sinThetaPointing_1:.4f}")
print(f"RASTART/DECSTART+Conv.\t{pointingCoords_2.ra.deg:.6f}\t{pointingCoords_2.dec.deg:.6f}"  
 + f"\t{pointingAltaz_2.alt.deg:.6f}\t{pointingAltaz_2.az.deg:.6f}\t{sinThetaPointing_2:.4f}")
print(f"RAEND/DECEND+Conv.\t{pointingCoords_3.ra.deg:.6f}\t{pointingCoords_3.dec.deg:.6f}"  
 + f"\t{pointingAltaz_3.alt.deg:.6f}\t{pointingAltaz_3.az.deg:.6f}\t{sinThetaPointing_3:.4f}")
print(f"Header START values\t{pointingCoords_header.ra.deg:.6f}\t{pointingCoords_header.dec.deg:.6f}"  
 + f"\t{altazCoords_header.alt.deg:.6f}\t{altazCoords_header.az.deg:.6f}\t{sinThetaPointing_header:.4f}")
print(f"Header END values\t{pointingCoords_header_2.ra.deg:.6f}\t{pointingCoords_header_2.dec.deg:.6f}"  
 + f"\t{altazCoords_header_2.alt.deg:.6f}\t{altazCoords_header_2.az.deg:.6f}\t{sinThetaPointing_header_2:.4f}")
print(f"Astrometry.net values\t{astroNetCoords.ra.deg:.6f}\t{astroNetCoords.dec.deg:.6f}"  
 + f"\t{astroNetAltaz.alt.deg:.6f}\t{astroNetAltaz.az.deg:.6f}\t{sinThetaAstroNet:.4f}")
print(f"Stellarium values\t{stellCoords.ra.deg:.6f}\t{stellCoords.dec.deg:.6f}"  
 + f"\t{stellAltaz_2.alt.deg:.6f}\t{stellAltaz_2.az.deg:.6f}\t{sinThetaStell:.4f}")
print(f"Stellarium values\t{stellCoords.ra.deg:.6f}\t{stellCoords.dec.deg:.6f}"  
 + f"\t{stellAltaz_2.alt.deg:.6f}\t{stellAltaz_2.az.deg:.6f}\t{sinThetaStell:.4f}")
```

```python
np.sin(Angle((90-36.8)*u.deg))
```

```python
def pixelLocation(center, object, rotpa):
    plateScale = 0.095695 # arcseconds/pixel
    centerPixel = np.array([2000.0, 2000.0])
    dd = (object.dec.deg - center.dec.deg) * 3600.0 / plateScale
    dr = (object.ra.deg - center.ra.deg) * np.cos(center.dec.rad) * 3600.0 / plateScale
    off = np.array([dr, dd])
    theta = np.radians(rotpa)
    c, s = np.cos(theta), np.sin(theta)
    # This is the rotation matrix that rotates (RA,Dec) into the detector coordinates
    R = np.array(((c, -s), (s, c))) 
    rotated_off = R.dot(off)
    location = centerPixel + rotated_off
    return location
```

```python
centerShiftRa = 800.0*u.arcsec
centerShiftDec = 0.0*u.arcsec
objectRa = 175.0*u.deg
objectDec = -85.0*u.deg
objectCoords = SkyCoord(ra=objectRa, dec=objectDec)
centerCoords = SkyCoord(ra=(objectRa + centerShiftRa), dec=(objectDec + centerShiftDec))
rotpa = Angle(mData['ROTPA']*u.deg)
shiftMag = 200.0
shift = np.array([[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]]) * shiftMag
dt = 300.0
deltaT = TimeDelta(dt, format='sec')
oldPixelCoords = pixelLocation(centerCoords, objectCoords, rotpa)
startPixelCoords = oldPixelCoords
currentTime = time
print(f"Start = {oldPixelCoords}")
fig = plt.figure(figsize=(16,8))
plt.subplot(1,2,1)
plt.title("Non-closure when doing shifts in AltAz", fontsize=18)
for step in range(4):
    centerAltAz = centerCoords.transform_to(AltAz(obstime=currentTime,location=location))
    shiftedAltAz = SkyCoord(AltAz(alt=centerAltAz.alt + shift[step][0]*u.arcsec, \
                                  az=centerAltAz.az + shift[step][1]*u.arcsec, \
                                  obstime=currentTime, location=location))
    centerCoords = shiftedAltAz.transform_to(ICRS)
    #print(centerCoords)
    newPixelCoords = pixelLocation(centerCoords, objectCoords, rotpa)
    print(f"Step {step} = {newPixelCoords}")
    currentTime = currentTime + deltaT
    delta = newPixelCoords - oldPixelCoords
    plt.arrow(x=oldPixelCoords[0], y=oldPixelCoords[1], dx=delta[0], dy=delta[1], \
              length_includes_head=True, head_width=75.0, color='red')
    oldPixelCoords = newPixelCoords
error = np.sqrt((newPixelCoords[0] - startPixelCoords[0])**2 + (newPixelCoords[1] - startPixelCoords[1])**2)
plt.xlim(0,4000)
plt.ylim(0,4000)
plt.text(200, 3800, f"Shift = {shiftMag} arcsec per step",fontsize=12)
plt.text(200, 3600, f"Time delay = {dt} seconds per step",fontsize=12)
plt.text(200, 3400, f"Final error = {error:.1f} pixels",fontsize=12)

print(f"Error = {error:.1f} pixels")

centerShiftRa = -1600.0*u.arcsec
centerShiftDec = -100.0*u.arcsec
objectCoords = SkyCoord(ra=objectRa, dec=objectDec)
centerCoords = SkyCoord(ra=(objectRa + centerShiftRa), dec=(objectDec + centerShiftDec))
shift = np.array([[10.0, 0.0], [0.0, 1.0], [-10.0, 0.0], [0.0, -1.0]]) * shiftMag
oldPixelCoords = pixelLocation(centerCoords, objectCoords, rotpa)
startPixelCoords = oldPixelCoords
currentTime = time
print(f"Start = {oldPixelCoords}")
plt.subplot(1,2,2)
plt.title("Closure when doing shifts in RaDec", fontsize=18)
for step in range(4):
    centerCoords = SkyCoord(ra=(centerCoords.ra + shift[step][0]*u.arcsec), \
                            dec=(centerCoords.dec + shift[step][1]*u.arcsec))
    #print(centerCoords)
    newPixelCoords = pixelLocation(centerCoords, objectCoords, rotpa)
    print(f"Step {step} = {newPixelCoords}")
    currentTime = currentTime + deltaT
    delta = newPixelCoords - oldPixelCoords
    plt.arrow(x=oldPixelCoords[0], y=oldPixelCoords[1], dx=delta[0], dy=delta[1], \
              length_includes_head=True, head_width=75.0, color='red')
    oldPixelCoords = newPixelCoords
error = np.sqrt((newPixelCoords[0] - startPixelCoords[0])**2 + (newPixelCoords[1] - startPixelCoords[1])**2)
plt.xlim(0,4000)
plt.ylim(0,4000)
plt.text(200, 3800, f"Shift = {shiftMag} arcsec per step",fontsize=12)
plt.text(200, 3600, f"Time delay = {dt} seconds per step",fontsize=12)
plt.text(200, 3400, f"Final error = {error:.1f} pixels",fontsize=12)

print(f"Error = {error:.1f} pixels")
plt.savefig(f"/project/cslage/AuxTel/offsets/Offset_Closure_AzEl_vs_RADec_23Apr21.pdf")
```

```python
altazCoords_header_2 = AltAz(alt=mData['ELEND']*u.deg, az=mData['AZEND']*u.deg, obstime=time, location=location)
```

```python
altazCoords_header_2.transform_to(GCRS)
```

```python
pointingCoords_2
```

```python
objectCoords
```

```python
objectCoords.transform_to(CIRS)
```

```python
objectCoords.transform_to(GCRS)
```

```python

```

```python
pointingCoords_2.transform_to(CIRS)
```

```python

```

```python
print(np.arcsin(0.7456)*180.0/np.pi, np.arcsin(0.756162127)*180.0/np.pi)
```

# Calculation of center pixel coordinates
# Need to do this more carefully if it matters

CRVAL1  =        240.141563967 / RA  of reference point                         
CRVAL2  =       -89.2964243801 / DEC of reference point                         
CRPIX1  =        1811.34913958 / X reference pixel                              
CRPIX2  =        1032.95088371 / Y reference pixel                              
CUNIT1  = 'deg     ' / X pixel scale units                                      
CUNIT2  = 'deg     ' / Y pixel scale units                                      
CD1_1   =   -1.96427829069E-05 / Transformation matrix                          
CD1_2   =    1.78943966674E-05 / no comment                                     
CD2_1   =   -1.79421818042E-05 / no comment                                     
CD2_2   =   -1.95647709172E-05 / no comment                                     
IMAGEW  =                 4072 / Image width,  in pixels.                       
IMAGEH  =                 4000 / Image height, in pixels.   
COMMENT scale: 0.0956118 arcsec/pix

```python
240.141563967 + (4072.0/2.0-1811.34913958)*0.0956118/3600.0
```

```python
-89.2964243801 - (4000.0/2.0-1032.95088371)*0.0956118/3600.0
```

```python

```
