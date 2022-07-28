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

## Analysis of multiple stuttered images

Craig Lage 10-May-22

```python
import sys, time, os, asyncio, glob
from datetime import datetime, timedelta
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import scipy.signal as signal
from astropy.time import Time, TimeDelta
from astropy.coordinates import AltAz, ICRS, EarthLocation, Angle, FK5, SkyCoord
import astropy.units as u
from lsst.obs.lsst.translators.latiss import AUXTEL_LOCATION
from lsst.daf.butler import Butler
```

```python
# This gives the azimuth tracking rate
def dAzdT(ALT, AZ):
    # I tried doing this anlatically, but then went to
    # this brute force method that is guaranteed to work.
    deltaT = 10.0
    time = datetime.now()
    shiftedTime = time + timedelta(seconds=deltaT)
    altAz = AltAz(obstime=time, location=AUXTEL_LOCATION)
    shiftedAltAz = AltAz(obstime=shiftedTime, location=AUXTEL_LOCATION)
    skyLocation = SkyCoord(AZ, ALT, frame=altAz)
    az1 = skyLocation.altaz.az.deg
    shiftedSkyLocation = skyLocation.transform_to(shiftedAltAz)
    az2 = shiftedSkyLocation.altaz.az.deg
    deltaAz = abs(az1 - az2)
    if deltaAz > 1.0:
        deltaAz = abs(deltaAz - 360.0)
    azSpeed = abs(deltaAz) / deltaT * 3600.0
    if azSpeed > 500.0:
        return 0.0
    else:
        return azSpeed
        
```

```python
butler = Butler('/repo/main', collections="LATISS/raw/all")
```

```python
# Plot all 5 images and calculate the centroid shifts
center = np.array([25.0,25.0])
shift = 50
xmin = 2655
yfirst = 1995
bgfirst = 1960
bgwidth = 30
bgxmin = 2550
bgxmax = 3000

xcentroid = []
ycentroid = []
altcentroid = []
azcentroid = []
times = []
expIds = [2022050500695, 2022050500696, 2022050500697, 2022050500698, 2022050500699]

plt.figure(figsize=(11, 8.5))
plt.suptitle(f"Stuttered Sequence - {expIds[0]} - {expIds[-1]}", fontsize = 18)
plotcounter = 0
for expId in expIds:
    exp = butler.get('raw', detector=0, exposure=expId)
    mData = exp.getMetadata()
    date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
    time = date_beg.utc
    el = Angle(mData['ELSTART'] * u.deg)
    az = Angle(mData['AZSTART'] * u.deg)
    dec = Angle(mData['DECSTART'] * u.deg)
    lat = AUXTEL_LOCATION.lat
    rotpa = Angle(mData['ROTPA']*u.deg)
    r = np.array(( (np.cos(rotpa), -np.sin(rotpa)),
               (np.sin(rotpa),  np.cos(rotpa)) ))

    # This calculates the angle theta between (X,Y) and (Az,El)
    sinTheta =  np.cos(lat) / np.cos(dec) * np.sin(az)
    cosTheta = (np.sin(el) * np.sin(dec) - np.sin(lat)) / (np.cos(el) * np.cos(dec))
    theta = Angle(np.arcsin(sinTheta))
    # The following removes the ambiguity in arcsin(theta)
    if cosTheta > 0:
        rotAzEl = rotpa - theta
    else:    
        rotAzEl = rotpa - (Angle(180.0 * u.deg) - theta)
    print(rotAzEl)
    
    for i in range(39):
        times.append(time.mjd * 86400.0)
        yfinish = yfirst - i * shift
        ystart = yfinish - shift
        bgstart = bgfirst - i * shift
        bgfinish = bgstart + bgwidth
        arr = exp.image.array[ystart:yfinish, xmin:xmin+shift]
        bg = exp.image.array[bgstart:bgfinish, bgxmin:bgxmax]
        background = np.nanmedian(bg)
        arr = arr - background
        arr = np.clip(arr, 0.1, 200000)
        nx = plotcounter % 20
        ny = int(plotcounter / 20)
        ax = plt.axes([0.05 + 0.045 * nx, 0.85 - 0.055 * ny, 0.050, 0.050], aspect = 1)
        ax.imshow(arr,   interpolation='Nearest', cmap='gray')
        ax.set_xticks([])
        ax.set_yticks([])
        plotcounter += 1


        xsum = 0
        ysum = 0
        imagesum = 0
        for ii in range(50):
            for jj in range(50):
                imagesum += arr[ii,jj]
                xsum += ii * arr[ii,jj]
                ysum += jj * arr[ii,jj]
        xsum /= imagesum
        ysum /= imagesum
        dx = np.array([xsum, ysum]) - center
        daltaz = r.dot(dx)
        altaz = daltaz + center
        xcentroid.append(xsum)
        ycentroid.append(ysum)
        altcentroid.append(altaz[0])
        azcentroid.append(altaz[1])
        time = time + TimeDelta(1.0, format='sec')
times = times - times[0]
plt.savefig(f"/project/cslage/AuxTel/stuttered/Stuttered_Multiple_{expIds[0]}_{expIds[-1]}.pdf")
```

```python
w = np.linspace(0.02, 0.20, 10)
coggingPeriod = 675.0 # arcseconds

plt.suptitle(f"Stuttered Images - {expIds[0]} - {expIds[-1]}")
plt.subplots_adjust(hspace=1.0, wspace=0.5)
plt.subplot(2,2,1)
plt.plot(times, altcentroid)
plt.xlabel("Time(sec)")
plt.ylabel("Alt Centroid(pixels)")
plt.subplot(2,2,3)
plt.title("Alt Periodogram")
# Periodogram
# Since the data is not all equally spaced, we need to use the 
# Lomb - Scargle periodogram instead of an FFT
pgram = signal.lombscargle(times, altcentroid, w, normalize=True)
azSpeed = dAzdT(el, az)
period = coggingPeriod / azSpeed
plt.plot(w, pgram)
plt.xlabel("Frequency (Hz)")
freq = 1.0/period
plt.plot([freq, freq], [0.0,0.2], color='red', ls = '--')
plt.text(0.10, 0.15, f"Cogging\nfrequency", color='red')
plt.subplot(2,2,2)
plt.plot(times, azcentroid)
plt.xlabel("Time(sec)")
plt.ylabel("Az Centroid(pixels)")
plt.subplot(2,2,4)
plt.title("Az Periodogram")
# Periodogram
# Since the data is not all equally spaced, we need to use the 
# Lomb - Scargle periodogram instead of an FFT
pgram = signal.lombscargle(times, azcentroid, w, normalize=True)
azSpeed = dAzdT(el, az)
period = coggingPeriod / azSpeed
plt.plot(w, pgram)
plt.xlabel("Frequency (Hz)")
freq = 1.0/period
plt.plot([freq, freq], [0.0,0.2], color='red', ls = '--')
plt.text(0.10, 0.15, f"Cogging\nfrequency", color='red')

plt.savefig(f"/project/cslage/AuxTel/stuttered/Stuttered_Centroids_AltAz.pdf")

```

```python

```
