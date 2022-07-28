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
from astropy.table import Table
from astroquery.astrometry_net import AstrometryNet

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
# Pick an expId, and compare this image with the next in the sequence.
myExpId = 2021031100406
for charVisit in charVisits:
    expId = charVisit['Visit'][0]
    if expId == myExpId:
        break
```

```python
cat = charVisit['brightCatalog']
```

```python
cat['base_SdssCentroid_x']
```

```python
cat['base_SdssShape_instFlux']
```

```python
cat.sort?
```

```python
test = cat.asAstropy()
test.keep_columns(['base_SdssCentroid_x', 'base_SdssCentroid_y', 'base_SdssShape_instFlux'])
print(test)
```

```python
test.sort('base_SdssShape_instFlux', reverse=True)
#test.remove_row(6)
print(test)
```

```python
ast = AstrometryNet()
ast.api_key = 'xxawwhvleirxcswx'

```

```python
AstrometryNet.show_allowed_settings()
```

```python
sources = cat.asAstropy()
sources.keep_columns(['base_SdssCentroid_x', 'base_SdssCentroid_y', 'base_SdssShape_instFlux'])
sources.sort('base_SdssShape_instFlux', reverse=True)
#sources.remove_row(6) # Hot pixel?
image_width = 4072
image_height = 4000
scale_units = 'arcsecperpix'
scale_type='ev'
scale_est = 0.095
scale_err = 2.0
center_ra = 241.467693
center_dec = -89.3087277
radius = 0.5
wcs_header = ast.solve_from_source_list(sources['base_SdssCentroid_x'], sources['base_SdssCentroid_y'],
                                        image_width, image_height, scale_units=scale_units,
                                        scale_type=scale_type, scale_est=scale_est, scale_err=scale_err,
                                        center_ra=center_ra, center_dec=center_dec, radius=radius,
                                        solve_timeout=240)
```

```python
wcs_header
```

```python
dir(wcs_header)
```

```python
dir(ast)
```

```python
ast._session_id
```

```python

```
