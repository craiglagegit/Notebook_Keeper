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

# Notebook for querying BOT data.

Initially written 27 May 2020 by Craig Lage\
Allows inspecting the image type and exposure time of the \
BOT images used for characterizing BF.\
Re-tested 08 Sep 2020 with latest code.

```python
! eups list -s | grep lsst_distrib
! eups list -s cp_pipe
```

```python
import sys, os, glob, subprocess
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter as gf
from scipy.signal import correlate2d as corr2d
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
```

```python
def removeLinearPart(image):
    (nx, ny) = image.shape
    m = nx
    X1, X2 = np.mgrid[:m, :m]
    X = np.hstack(   ( np.reshape(X1, (m*m, 1)) , np.reshape(X2, (m*m, 1)) ) )
    X = np.hstack(   ( np.ones((m*m, 1)) , X ))
    YY = np.reshape(image, (m*m, 1))
    theta = np.dot(np.dot( np.linalg.pinv(np.dot(X.transpose(), X)), X.transpose()), YY)
    plane = np.reshape(np.dot(X, theta), (m, m))
    return image - plane
    
```

```python
corr_images = {}
nx = 100
ny = 100
xmin = 150
xmax = xmin + nx
ymin = 500
ymax = ymin + ny

sigma = 6.0
ampNum = 16
```

```python
# First, a pair of random Poisson images
level = 150000
images = []
for n in range(2):
    im = np.zeros([nx, ny])
    for i in range(nx):
        for j in range(ny):
            im[i,j] = np.random.poisson(lam=level)
    im_median = np.median(im)
    im = removeLinearPart(im)
    im = im / im_median * 100.0
    smoothed_im = gf(im, sigma = sigma)
    print(smoothed_im.min(), smoothed_im.max())
    images.append(smoothed_im)

corr_images["Poisson"] = [images, level]
```

```python
# Next, images from the 13 raft run
DATA_DIR = '/project/shared/BOT/' 
filedir = DATA_DIR+'_parent/raw/'
rafts = ['R12', 'R02']
names = ['BOT-12543-R12S02-E2V', 'BOT-12543-R12S02-ITL']
SENSOR = 'S02'
for n, RAFT in enumerate(rafts):
    files = glob.glob(filedir+'*/*/3020090200370-%s-%s-det???.fits'%(RAFT,SENSOR))
    files += glob.glob(filedir+'*/*/3020090200371-%s-%s-det???.fits'%(RAFT,SENSOR))
    #files = glob.glob(filedir+'*/*/3020090200361-%s-%s-det???.fits'%(RAFT,SENSOR))
    #files += glob.glob(filedir+'*/*/3020090200362-%s-%s-det???.fits'%(RAFT,SENSOR))
    files = np.sort(files)
    numFiles = len(files)
    print(numFiles)
    images = []
    for file in files:
        hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
        phdr=hdulist[0].header
        filenumber = file.split('/')[-1][0:13]
        seq = int(file.split('/')[-1][8:13])
        exptime = phdr['EXPTIME']
        imgtype = phdr['IMGTYPE'] 
        print("%s\t%s\t%f\t%d"%(filenumber, imgtype, exptime, seq))
        im = hdulist[ampNum].data[ymin:ymax, xmin:xmax]
        im_median = np.median(im)
        im = removeLinearPart(im)
        im = im / im_median * 100.0
        smoothed_im = gf(im, sigma = sigma)
        print("Median = %.1f"%im_median, smoothed_im.min(), smoothed_im.max())
        images.append(smoothed_im)
    corr_images[names[n]] = [images, im_median]

```

```python
# Next, images from the 9 raft run
DATA_DIR = '/project/shared/BOT/' 
filedir = DATA_DIR+'_parent/raw/'
rafts = ['R12', 'R02']
names = ['BOT-6790D-R12S02-E2V', 'BOT-6790D-R12S02-ITL']
SENSOR = 'S02'
for n, RAFT in enumerate(rafts):
    files = glob.glob(filedir+'*/*/3019101300054-%s-%s-det???.fits'%(RAFT,SENSOR))
    files += glob.glob(filedir+'*/*/3019101300055-%s-%s-det???.fits'%(RAFT,SENSOR))
    files = np.sort(files)
    numFiles = len(files)
    print(numFiles)
    images = []
    for file in files:
        hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
        phdr=hdulist[0].header
        filenumber = file.split('/')[-1][0:13]
        seq = int(file.split('/')[-1][8:13])
        exptime = phdr['EXPTIME']
        imgtype = phdr['IMGTYPE'] 
        print("%s\t%s\t%f\t%d"%(filenumber, imgtype, exptime, seq))
        im = hdulist[ampNum].data[ymin:ymax, xmin:xmax]
        im_median = np.median(im)
        im = removeLinearPart(im)
        im = im / im_median * 100.0
        smoothed_im = gf(im, sigma = sigma)
        print("Median = %.1f"%im_median, smoothed_im.min(), smoothed_im.max())
        images.append(smoothed_im)
    corr_images[names[n]] = [images, im_median]

```

```python
# Next, images from ComCam
DATA_DIR = '/project/shared/comCam/'
filedir = DATA_DIR+'_parent/raw/'
rafts = ['R22']
names = ['ComCam-2020-08-13-ITL']
SENSOR = 'S02'
for n, RAFT in enumerate(rafts):
    
    filedir = DATA_DIR+'_parent/raw/'
    files = glob.glob(filedir+'*/*/2020081300054-%s-%s-det???.fits'%(RAFT,SENSOR))
    files += glob.glob(filedir+'*/*/2020081300055-%s-%s-det???.fits'%(RAFT,SENSOR))
    files = np.sort(files)
    numFiles = len(files)
    print(numFiles)
    images = []
    for file in files:
        hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
        phdr=hdulist[0].header
        filenumber = file.split('/')[-1][0:13]
        seq = int(file.split('/')[-1][8:13])
        exptime = phdr['EXPTIME']
        imgtype = phdr['IMGTYPE'] 
        print("%s\t%s\t%f\t%d"%(filenumber, imgtype, exptime, seq))
        im = hdulist[ampNum].data[ymin:ymax, xmin:xmax]
        im_median = np.median(im)
        im = removeLinearPart(im)
        im = im / im_median * 100.0
        smoothed_im = gf(im, sigma = sigma)
        print("Median = %.1f"%im_median, smoothed_im.min(), smoothed_im.max())
        images.append(smoothed_im)
    corr_images[names[n]] = [images, im_median]

```

```python
# Next, images from UCDavis E2V
DATA_DIR = '/project/bootcamp/cslage/e2v_fits_files/flats/20190709_e2v_flats/'
name = 'UCD-2019-07-09-E2V'


files = glob.glob(DATA_DIR +'E2V-CCD250-112-09_flat_flat_348_*.fits')
files += glob.glob(DATA_DIR +'E2V-CCD250-112-09_flat_flat_351_*.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)
images = []
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    phdr=hdulist[0].header
    exptime = phdr['EXPTIME']
    imgtype = phdr['IMGTYPE'] 
    print("%s\t%s"%(imgtype, exptime))
    im = pf.getdata(file,ampNum)[ymin:ymax, xmin:xmax]
    im_median = np.median(im)
    im = removeLinearPart(im)
    im = im / im_median * 100.0
    smoothed_im = gf(im, sigma = sigma)
    print("Median = %.1f"%im_median, smoothed_im.min(), smoothed_im.max())
    images.append(smoothed_im)
corr_images[name] = [images, im_median]

```

```python
# Next, images from UCDavis ITL
DATA_DIR = '/project/bootcamp/cslage/itl_fits_files/flats/20180525_002_flats/'
name = 'UCD-2018-05-25-ITL'


files = glob.glob(DATA_DIR +'ITL-3800C-002_flat_flat_1300_*.fits')
files += glob.glob(DATA_DIR +'ITL-3800C-002_flat_flat_1301_*.fits')
files = np.sort(files)
numFiles = len(files)
print(numFiles)
images = []
for file in files:
    hdulist = pf.open(file, mode='readonly', do_not_scale_image_data=True)
    phdr=hdulist[0].header
    exptime = phdr['EXPTIME']
    imgtype = phdr['IMGTYPE'] 
    print("%s\t%s"%(imgtype, exptime))
    im = pf.getdata(file,ampNum)[ymin:ymax, xmin:xmax]
    im_median = np.median(im)
    im = removeLinearPart(im)
    im = im / im_median * 100.0
    smoothed_im = gf(im, sigma = sigma)
    print("Median = %.1f"%im_median, smoothed_im.min(), smoothed_im.max())
    images.append(smoothed_im)
corr_images[name] = [images, im_median]

```

```python

print(corr_images.keys())

```

```python
for key in corr_images.keys():
    [images, flux] = corr_images[key]
    fig = plt.figure(figsize = (16,8))

    plotCounter = 1
    corr = corr2d(images[0], images[1])
    scalar_corr = corr.max()
    for image in images:
        plt.subplot(1,2,plotCounter)
        plt.imshow(image, vmin = -0.10, vmax = 0.10)
        plotCounter += 1
    plt.suptitle(key + " Flux = %.1f, Scalar_correlation = %.3f"%(flux,scalar_corr), fontsize = 18)
    print(scalar_corr)

```

```python

```

```python
# First, a pair of random Poisson images
for m in range(1,40):
    level = 10000 * m
    images = []
    for n in range(2):
        im = np.zeros([nx, ny])
        for i in range(nx):
            for j in range(ny):
                im[i,j] = np.random.poisson(lam=level)
        im_median = np.median(im)
        im = (im - im_median) / im_median * 100.0
        smoothed_im = gf(im, sigma = sigma)
        #print(im_median, smoothed_im.min(), smoothed_im.max())
        images.append(smoothed_im)
    corr = corr2d(images[0], images[1])
    scalar_corr = corr.max()
    print(level, scalar_corr)
```

```python

```
