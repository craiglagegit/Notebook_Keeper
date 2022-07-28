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

## AuxTel AzEl offsets - 19-Apr-21

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

```
