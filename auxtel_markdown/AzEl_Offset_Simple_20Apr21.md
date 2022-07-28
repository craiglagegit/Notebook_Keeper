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

## AuxTel AzEl offsets - 20-Apr-21

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
from astropy.coordinates import AltAz, ICRS, EarthLocation, Angle, FK5
import astropy.units as u

from lsst.daf.butler import Butler as gen3Butler
from lsst.daf.persistence import Butler as gen2Butler
from lsst_efd_client import EfdClient
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python
# Set Cerro Pachon location
location = EarthLocation.from_geodetic(lon=-70.749417*u.deg,
                                       lat=-30.244639*u.deg,
                                       height=2663.0*u.m)

```

```python
infile = open('/project/cslage/AuxTel/offsets/offsets_16apr21.pkl','rb')
charVisits = pkl.load(infile)
infile.close()
```

```python
# This helps make the plots more compact
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
def rotAngle(charVisits, myExpId, location=location):
    for charVisit in charVisits:
        visit = charVisit['Visit']
        expId = visit[0]
        if expId == myExpId:
            break
    rotpa = Angle(visit[6]*u.deg)
    az = Angle(visit[9]*u.deg)
    #az = Angle(179.2616*u.deg)
    el = Angle(visit[10]*u.deg)
    #el = Angle(30.3307*u.deg)
    ra = Angle(visit[11]*u.deg)
    dec = Angle(visit[12]*u.deg)
    #dec = Angle(-89.3566*u.deg)
    sinTheta = np.cos(location.lat) * np.sin(az) / np.cos(dec)
    print(sinTheta)
    print(Angle(np.arcsin(sinTheta)).deg)
    #print(az.deg, el.deg, ra.deg, dec.deg)
    cosTheta = (np.sin(location.lat) - np.sin(el) * np.sin(dec)) / np.cos(el * np.cos(dec))
    #print(cosTheta)
    theta = Angle(np.arccos(cosTheta))
    #print(theta.deg, rotpa.deg, (theta+rotpa).deg)
    return #theta + rotpa
```

```python
rotAngle(charVisits, 2021031100422)
```

```python
el = Angle(179.2616*u.deg)
print(np.sin(el))
el = Angle((180.0 - 179.2616)*u.deg)
print(np.sin(el))
```

```python

point_in_time = pd.to_datetime('2021-03-12T04:37:36.415').tz_localize('UTC')


```

```python
# not what you want
local_time = point_in_time.tz_convert("America/Santiago")
print(local_time)
```

```python
# Pick an expId, and compare this image with the next in the sequence.
myExpId = 2021031100422
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
#theta = rotAngle(charVisits, 2021031100422).rad
#theta = -np.radians(rotpa)
theta = np.radians(37)
c, s = np.cos(theta), np.sin(theta)
# This is the rotation matrix that puts the commanded offsets into the detector coordinates
R = np.array(((c, -s), (s, c))) 
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
#plt.savefig(f"/project/cslage/AuxTel/offsets/Offsets_Meas_vs_Expected_{myExpId}_19Apr21.pdf")
```

```python
# Plot the problem
expId = 2021031100422
for charVisit in charVisits:
    visit = charVisit['Visit']
    expId = visit[0]
    if expId == myExpId:
        break
az = visit[9]
dec = -89.30854#visit[12]

r = 90.0 + dec
az_angle = 180.0 - az
zenith = -90.0 - location.lat.deg
dx = (zenith + 1.0) * np.tan(az_angle*np.pi / 180.0)

fig = plt.figure(figsize = (16,8))
plt.suptitle(f"Header values in expId {expId} are inconsistent", fontsize = 24)
ax1 = plt.subplot(1,2,1,aspect='equal')
decCircle = plt.Circle((0,0), r, fill=False)
poleCircle = plt.Circle((0,0), 0.05, fill=True)
ax1.add_patch(decCircle)
ax1.add_patch(poleCircle)
ax1.plot([0.0,0.0], [zenith, 1.0])
ax1.plot([0.0,dx], [zenith, 1.0])
ax1.text(0.1,0.2, "SCP", fontsize=12)
ax1.text(0.0,1.2, "Meridian", fontsize=12)
ax1.text(-1.8,1.2, "Azimuth = %.5f"%az, fontsize=12)
ax1.text(0.8,0.0, "Dec = %.5f"%dec, fontsize=12)
ax1.set_xlim(-2.0,2.0)
ax1.set_ylim(-2.0,2.0)
ax2 = plt.subplot(1,2,2,aspect='equal')
decCircle = plt.Circle((0,0), r, fill=False)
poleCircle = plt.Circle((0,0), 0.05, fill=True)
ax2.add_patch(decCircle)
ax2.add_patch(poleCircle)
ax2.plot([0.0,0.0], [zenith, 1.0])
ax2.plot([0.0,dx], [zenith, 1.0])
ax2.text(2.0,-1.0, "SCP", fontsize=12)
ax2.text(2.0,zenith, "Zenith", fontsize=12)
ax2.set_xlim(-36.0,36.0)
ax2.set_ylim(-70.0,2.0)
plt.savefig(f"/project/cslage/AuxTel/offsets/Header_Issue_{expId}_20Apr21.pdf")
```

```python
# Plot the problem
expId = 2021031100422
for charVisit in charVisits:
    visit = charVisit['Visit']
    myExpId = visit[0]
    if expId == myExpId:
        break
az = visit[9]
dec = visit[12]

r = 90.0 + dec
az_angle = 180.0 - az
zenith = 90.0 + location.lat.deg
dx = -(zenith + 1.0) * np.tan(az_angle*np.pi / 180.0)
dy = -location.lat.deg
fig = plt.figure(figsize = (16,8))
plt.suptitle(f"Header values in expId {expId} are inconsistent", fontsize = 24)
ax1 = plt.subplot(1,2,1,aspect='equal')
decCircle = plt.Circle((0,dy), r, fill=False)
poleCircle = plt.Circle((0,dy), 0.05, fill=True)
ax1.add_patch(decCircle)
ax1.add_patch(poleCircle)
ax1.plot([0.0,0.0], [zenith+dy, -1.0+dy])
ax1.plot([0.0,dx], [zenith+dy, -1.0+dy])
ax1.text(0.1,-0.2+dy, "SCP", fontsize=12)
ax1.text(0.0,-1.2+dy, "Meridian", fontsize=12)
ax1.text(-1.8,-1.2+dy, "Azimuth = %.5f"%az, fontsize=12)
ax1.text(0.8,0.0+dy, "Dec = %.5f"%dec, fontsize=12)
ax1.set_xlim(-2.0,2.0)
ax1.set_ylim(2.0+dy,-2.0+dy)
ax1.set_ylabel("Altitude")
ax2 = plt.subplot(1,2,2,aspect='equal')
decCircle = plt.Circle((0,dy), r, fill=False)
poleCircle = plt.Circle((0,dy), 0.05, fill=True)
ax2.add_patch(decCircle)
ax2.add_patch(poleCircle)
ax2.plot([0.0,0.0], [zenith+dy, -1.0+dy])
ax2.plot([0.0,dx], [zenith+dy, -1.0+dy])
ax2.text(2.0,1.0+dy, "SCP", fontsize=12)
ax2.text(2.0,zenith+dy, "Zenith", fontsize=12)
ax2.set_xlim(-36.0,36.0)
ax2.set_ylim(65.0+dy,-2.0+dy)
ax2.set_ylabel("Altitude")
plt.savefig(f"/project/cslage/AuxTel/offsets/Header_Issue_{expId}_20Apr21.pdf")
```

```python

```
