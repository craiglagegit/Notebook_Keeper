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

## Bad ComCam files

Investigating file corruption reported by Lupton\
Craig Lage - 07-Sep-21

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
```

```python
# These are the ones Lupton reported an issue
# The discrepancy always begins at 14,128
file1 = '/lsstdata/offline/instrument/LSSTComCam/storage/2021-08-19/CC_O_20210819_000003-R22S01.fits'
file2 = '/lsstdata/offline/instrument/LSSTComCam-ccs/storage/20210819/CC_O_20210819_000003/CC_O_20210819_000003_R22_S01.fits'
```

```python
amp = 3
hdulist1 = pf.open(file1, mode='readonly', do_not_scale_image_data=True)
dat1 = hdulist1[amp].data
hdulist2 = pf.open(file2, mode='readonly', do_not_scale_image_data=True)
dat2 = hdulist2[amp].data

```

```python
# These are discrepant starting at (14.128)
match = True
for i in range(2048):
    for j in range(576):
        if dat1[i,j] != dat2[i,j]:
            print(i,j, dat1[i,j], dat2[i,j])
            match = False
            break
        else:
            continue
    if match:
        continue
    else:
        break
```

```python
from matplotlib.colors import LogNorm
# Now let's look at ithem
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

plt.figure(figsize=(16,8))
plt.subplot(1,2,1)
plt.title(f"{file1}\n Mean = {dat1.mean()}; Std = {dat1.std()}",fontsize=8)
img = plt.imshow(dat1, norm=LogNorm(vmin=23000, vmax=23500), interpolation='Nearest', cmap='gray')
colorbar(img)
plt.subplot(1,2,2)
plt.title(f"{file2}\n Mean = {dat2.mean()}; Std = {dat2.std()}",fontsize=8)
img = plt.imshow(dat2, norm=LogNorm(vmin=23000, vmax=23500), interpolation='Nearest', cmap='gray')
colorbar(img)

plt.tight_layout(h_pad=1)
plt.savefig("/project/cslage/ComCam/bad_files/20210819_000003-R22S01.png")
```

```python
# These are new files
# They look OK
file1 = '/lsstdata/offline/instrument/LSSTComCam/storage/2021-09-16/CC_O_20210916_000060-R22S02.fits'
file2 = '/lsstdata/offline/instrument/LSSTComCam-ccs/storage/20210916/CC_O_20210916_000060/CC_O_20210916_000060_R22_S02.fits'
```

```python
amp = 3
hdulist1 = pf.open(file1, mode='readonly', do_not_scale_image_data=True)
dat1 = hdulist1[amp].data
hdulist2 = pf.open(file2, mode='readonly', do_not_scale_image_data=True)
dat2 = hdulist2[amp].data

```

```python
# These look OK
match = True
for i in range(2048):
    for j in range(576):
        if dat1[i,j] != dat2[i,j]:
            print(i,j, dat1[i,j], dat2[i,j])
            match = False
            break
        else:
            continue
    if match:
        continue
    else:
        break
```

```python
from matplotlib.colors import LogNorm
# Now let's look at ithem
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

plt.figure(figsize=(16,8))
plt.subplot(1,2,1)
plt.title(f"{file1}\n Mean = {dat1.mean()}; Std = {dat1.std()}",fontsize=8)
img = plt.imshow(dat1, norm=LogNorm(vmin=21000, vmax=23000), interpolation='Nearest', cmap='gray')
colorbar(img)
plt.subplot(1,2,2)
plt.title(f"{file2}\n Mean = {dat2.mean()}; Std = {dat2.std()}",fontsize=8)
img = plt.imshow(dat2, norm=LogNorm(vmin=21000, vmax=23000), interpolation='Nearest', cmap='gray')
colorbar(img)

plt.tight_layout(h_pad=1)
plt.savefig("/project/cslage/ComCam/bad_files/20210916_000060-R22S02.png")
```

```python

```
