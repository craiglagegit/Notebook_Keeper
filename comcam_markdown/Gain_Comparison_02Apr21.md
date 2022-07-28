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
! eups list -s | grep lsst_distrib
! eups list -s | grep cp_pipe
```

```python
import sys, os, glob, time
import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
from lsst.daf.persistence import Butler
```

```python
DIR_1 = '/project/shared/comCam/rerun/cslage/PTC_2021-04-01/'
DIR_2 = '/project/shared/comCam/rerun/cslage/PTC_2021-04-02/'
DIR_3 = '/project/shared/comCam-CCS/rerun/cslage/PTC_2020-12-29/'
RAFT = 'R22'
#dirs = [DIR_1, DIR_2, DIR_3]
dirs = [DIR_2, DIR_3]
#dirs = [DIR_3]
#run_names = ['Summit-2021-04-01-Gband', 'Summit-2021-04-02-Iband', 'Summit-2020-12-29-Rband']
run_names = ['Summit-2021-04-02-Iband', 'Summit-2020-12-29-Rband']
#run_names = ['Summit-2020-12-29-Iband?']
markers = ['o', 'x', '+', '*', '^', 'v']
```

```python
plt.figure(figsize=(16,16))
plt.subplots_adjust(hspace=0.5)
plt.suptitle("ComCam Gains - Summit", fontsize=24)
for i,dir in enumerate(dirs):
    butler = Butler(dir)
    run_name = run_names[i]
    plotcounter = 0
    for detector in range(9):
        plotcounter += 1
        plt.subplot(3,3,plotcounter)
        plt.title("Detector%d"%detector, fontsize = 12)
        #ptcDataset = butler.get('photonTransferCurveDataset', raftName=RAFT, detector=detector)
        #gain_data = ptcDataset.gain
        #gain_err_data = ptcDataset.gainErr
        #amps = gain_data.keys()
        datasetFile = dir+'/calibrations/ptc/ptcDataset-det%03d.fits'%detector
        hdulist = pf.open(datasetFile, mode='readonly', do_not_scale_image_data=True)
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
        plt.errorbar(amp_nums, gains, yerr=gain_err, marker = markers[i], label = run_name)
        plt.ylim(1.0, 2.0)
        plt.ylabel("Gain", fontsize = 12)
        plt.xticks(amp_nums,names, fontsize=8)
        plt.legend(loc = 'upper right', fontsize = 12)
plt.savefig(DIR_2+'plots/Gain_Summary_06Apr20.pdf')

```

```python

```
