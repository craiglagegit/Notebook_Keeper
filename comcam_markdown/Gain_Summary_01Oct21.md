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

# Notebook for plotting ComCam gains.

Initially written 28 Aug 2020 by Craig Lage.

```python
import sys, os, glob, time
import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from lsst.daf.persistence import Butler
```

```python
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']

plt.figure(figsize=(16,16))
plt.subplots_adjust(hspace=0.5)
plt.suptitle("ComCam Gains - 20210930 - R band", fontsize=24)
plotcounter = 0
for detector, sensor in enumerate(sensors):
    plotcounter += 1
    filename = "/repo/main/u/cslage/ptc_20210930/20211001T194239Z/ptc/ptc_LSSTComCam_R22_%s_u_cslage_ptc_20210930_20211001T194239Z.fits"%sensor
    #filename = "/repo/main/u/cslage/ptc_20210402A/20210413T170828Z/ptc/ptc_LSSTComCam_R22_%s_u_cslage_ptc_20210402A_20210413T170828Z.fits"%sensor
    hdulist = pf.open(filename, mode='readonly', do_not_scale_image_data=True)
    plt.subplot(3,3,plotcounter)
    plt.title("Detector%d"%detector, fontsize = 12)
    hdulist = pf.open(filename, mode='readonly', do_not_scale_image_data=True)
    data=hdulist[1].data
    gain_data = data['GAIN']
    gain_err_data = data['GAIN_ERR']
    gains = []
    gain_err = []
    names = []
    amp_nums = list(range(16))
    for ii in amp_nums:
        gains.append(gain_data[ii])
        gain_err.append(gain_err_data[ii])
        names.append(data['AMPLIFIER_NAME'][ii])
    plt.errorbar(amp_nums, gains, yerr=gain_err)
    plt.ylim(1.0, 1.8)
    plt.ylabel("Gain", fontsize = 12)
    plt.xticks(amp_nums,names, fontsize=8)
    #plt.legend(loc = 'upper right', fontsize = 12)
plt.savefig('/repo/main/u/cslage/ptc_20210930/plots/Gain_Summary_09Sep21.pdf')

```

```python

```
