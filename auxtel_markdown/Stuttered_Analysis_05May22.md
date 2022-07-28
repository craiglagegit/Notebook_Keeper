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

## Stuttered image analysis

Craig Lage 10-May-22

```python
import sys, time, os, asyncio, glob
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from astropy.time import Time, TimeDelta
from astropy.coordinates import AltAz, ICRS, EarthLocation, Angle, FK5, SkyCoord
import astropy.units as u
from lsst.obs.lsst.translators.latiss import AUXTEL_LOCATION
from lsst.daf.butler import Butler
```

```python
butler = Butler('/repo/main', collections="LATISS/raw/all")
```

```python
# Now get the image data and the metadata

expId = 2022050500699
exp = butler.get('raw', detector=0, exposure=expId)
mData = exp.getMetadata()
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
time = date_beg.utc
print(expId, time)
```

```python
plt.figure(figsize=(8,8))
arr = exp.image.array
img = plt.imshow(arr, norm=LogNorm(vmin=10000, vmax=25000),  interpolation='Nearest', cmap='gray')

```

```python
# Manually find a good bounding box
shift = 50
xmin = 2655
yfirst = 1995
yfinish = yfirst - 38 * shift
ystart = yfinish - shift
bgfirst = 1960
bgwidth = 30
bgstart = bgfirst - 38 * shift
bgfinish = bgstart + bgwidth
bgxmin = 2550
bgxmax = 3000


plt.figure(figsize=(4,4))
arr = exp.image.array[ystart:yfinish, xmin:xmin+shift]
bg = exp.image.array[bgstart:bgfinish, bgxmin:bgxmax]
background = np.nanmedian(bg)
print(f"Background = {background}")
arr = arr - background
print(f"Min = {arr.min()}, Max = {arr.max()}")
arr = np.clip(arr, 0.1, 200000)

img = plt.imshow(arr,   interpolation='Nearest', cmap='gray')
plt.colorbar()
```

```python
# Now run the whole sequence
shift = 50
xmin = 2655
yfirst = 1995
bgfirst = 1960
bgwidth = 30
bgxmin = 2550
bgxmax = 3000


plt.figure(figsize=(11, 8.5))
plt.suptitle(f"Stuttered Sequence - {expId}", fontsize = 18)
xcentroid = []
ycentroid = []
for i in range(39):
    yfinish = yfirst - i * shift
    ystart = yfinish - shift
    bgstart = bgfirst - i * shift
    bgfinish = bgstart + bgwidth
    arr = exp.image.array[ystart:yfinish, xmin:xmin+shift]
    bg = exp.image.array[bgstart:bgfinish, bgxmin:bgxmax]
    background = np.nanmedian(bg)
    #print(background)
    arr = arr - background
    arr = np.clip(arr, 0.1, 200000)
                           
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
    xcentroid.append(xsum)
    ycentroid.append(ysum)

    nx = i % 8
    ny = int(i / 8)
    ax = plt.axes([0.05 + 0.12 * nx, 0.80 - 0.14 * ny, 0.118, 0.118], aspect = 1)
    ax.imshow(arr,   interpolation='Nearest', cmap='gray')
    #plotcounter += 1
    ax.set_xticks([])
    ax.set_yticks([])
plt.savefig(f"/project/cslage/AuxTel/stuttered/Stuttered_BG_{expId}.pdf")
    
```

```python
plt.subplots_adjust(wspace=0.5)
plt.subplot(1,2,1)
plt.plot(xcentroid)
plt.xlabel("Time(sec)")
plt.ylabel("XCentroid(pixels)")
plt.subplot(1,2,2)
plt.plot(ycentroid)
plt.xlabel("Time(sec)")
plt.ylabel("YCentroid(pixels)")

plt.savefig(f"/project/cslage/AuxTel/stuttered/Stuttered_Centroids_BG_{expId}.pdf")

```

```python
exp.getWcs().getPixelOrigin()
```

```python
exp.getWcs?
```

```python
import lsst.afw.geom as geom
```

```python
exp.getWcs?
```

```python

```
