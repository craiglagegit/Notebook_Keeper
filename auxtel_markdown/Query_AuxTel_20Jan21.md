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

```python
import sys, os, glob, subprocess
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
```

```python
DATA_DIR = '/project/shared/auxTel/'
filedir = DATA_DIR+'_parent/raw/'
files = glob.glob(filedir+'2021-02-16/*00209-det000.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles, files[-1])
```

```python jupyter={"outputs_hidden": true}
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    hdr=hdulist[0].header
    for key in hdr.keys():
        print(key, hdr[key])

```

```python
bad = [164,165]+list(range(208,222))+list(range(252,286))
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    hdr=hdulist[0].header
    obj = hdr['OBJECT']
    seq = hdr['OBSID'].split('_')[3]
    exptime = hdr['EXPTIME']
    imgtype = hdr['IMGTYPE']
    domeaz = hdr['DOMEAZ']
    if domeaz is None:
        domeaz= 0.0
    rotpa = hdr['ROTPA']
    delaz = hdr['AZSTART'] - hdr['AZEND']
    delel = hdr['ELSTART'] - hdr['ELEND']
    #print(seq,exptime,imgtype,domeaz,rotpa,delaz,delel)
    #break
    if int(seq) in bad:
        print("\033[91m"+"%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f"%(seq,exptime,imgtype,domeaz,rotpa,delaz,delel))
    else:
        print("\033[90m"+"%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f"%(seq,exptime,imgtype,domeaz,rotpa,delaz,delel))
```

```python
bad = list(range(106,157))+[305,306]+list(range(321,329))
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    hdr=hdulist[0].header
    seq = hdr['OBSID'].split('_')[3]
    exptime = hdr['EXPTIME']
    imgtype = hdr['IMGTYPE']
    domeaz = hdr['DOMEAZ']
    if domeaz is None:
        domeaz= 0.0
    rotpa = hdr['ROTPA']
    delaz = hdr['AZSTART'] - hdr['AZEND']
    delel = hdr['ELSTART'] - hdr['ELEND']
    #print(seq,exptime,imgtype,domeaz,rotpa,delaz,delel)
    #break
    if int(seq) in bad:
        print("\033[91m"+"%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f"%(seq,exptime,imgtype,domeaz,rotpa,delaz,delel))
    else:
        print("\033[90m"+"%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f"%(seq,exptime,imgtype,domeaz,rotpa,delaz,delel))
```

```python

```

```python
DATA_DIR = '/project/shared/auxTel/'
filedir = DATA_DIR+'_parent/raw/'
files = glob.glob(filedir+'2021-01-26/*-det000.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles, files[-1])
```

```python
bad = range(90,114)
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    hdr=hdulist[0].header
    seq = hdr['OBSID'].split('_')[3]
    exptime = hdr['EXPTIME']
    imgtype = hdr['IMGTYPE']
    domeaz = hdr['DOMEAZ']
    rotpa = hdr['ROTPA']
    delaz = hdr['AZSTART'] - hdr['AZEND']
    delel = hdr['ELSTART'] - hdr['ELEND']
    if int(seq) in bad:
        print("\033[91m"+"%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f"%(seq,exptime,imgtype,domeaz,rotpa,delaz,delel))
    else:
        print("\033[90m"+"%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6f"%(seq,exptime,imgtype,domeaz,rotpa,delaz,delel))
```
