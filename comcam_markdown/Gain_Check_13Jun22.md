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

# Checking gain and noise values on all three cameras

Initially written 13 Jul 2022 by Craig Lage

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.io.fits as pf
from lsst.daf.butler import Butler
import lsst.afw.math as afwMath
```

```python
pdf = PdfPages("/project/cslage/gain_noise/Gain_Noise_Summary_13Jul22.pdf")
```

```python jupyter={"outputs_hidden": true} tags=[]
butler = Butler("/repo/main", collections=["LSSTCam/raw/all", "LSSTCam/calib/unbounded", \
                                          "u/cslage/bps_13144M"])
LSSTCam_camera = butler.get('camera', instrument='LSSTCam')
```

```python tags=[]
expId = 3021120600575
plotcounter = 0
fig1 = plt.figure(figsize=(16,16))
plt.subplots_adjust(hspace=0.3)

for plotData in ["Gain", "Noise"]:
    plt.suptitle("LSSTCam Gain and Noise - 13Jul22", fontsize=24)
    for detector in [54, 55, 56, 73, 74, 77]:
        plotcounter += 1
        plt.subplot(4,3,plotcounter)
        plt.title("%s"%LSSTCam_camera[detector].getName(), fontsize = 12)
        ptcDataset = butler.get('ptc', exposure=expId, detector=detector)
        if plotData == "Gain":
            data = ptcDataset.gain
            err_data = ptcDataset.gainErr
        else:
            data = ptcDataset.noise
            err_data = ptcDataset.noiseErr
        amps = data.keys()

        yvals = []
        stored_yvals = []
        err = []
        amp_nums = []
        for ii, amp in enumerate(amps):
            yvals.append(data[amp])
            err.append(err_data[amp])
            amp_nums.append(ii)
            amps = data.keys()
            if plotData == "Gain":
                stored_yvals.append(LSSTCam_camera[detector].getAmplifiers()[ii].getGain())
            else:
                stored_yvals.append(LSSTCam_camera[detector].getAmplifiers()[ii].getReadNoise())
        plt.errorbar(amp_nums, yvals, yerr=err, label="Meas_13144")
        plt.plot(amp_nums, stored_yvals, marker='x', label="DM Camera Values")
        if plotData == "Gain":
            plt.ylim(1.0, 2.4)
            plt.ylabel("Gain (e-/ADU)", fontsize = 12)
        else:
            plt.ylim(0.0,20.0)
            plt.ylabel("Read noise (electrons)", fontsize = 12)
        plt.xticks(amp_nums,amps, fontsize=8)
        plt.legend(loc = 'upper right', fontsize = 12)
pdf.savefig(fig1)
```

```python tags=[]
butler = Butler("/repo/main", collections=["LSSTComCam/raw/all", "LSSTComCam/calib/unbounded",  \
                                           "u/cslage/comcam/ptc_20220218"])
LSSTComCam_camera = butler.get('camera', instrument='LSSTComCam')
```

```python tags=[]
expId = 2022021800078
plotcounter = 0
fig2 = plt.figure(figsize=(16,16))
plt.subplots_adjust(hspace=0.3)

for plotData in ["Gain", "Noise"]:
    plt.suptitle("LSSTComCam Gain and Noise - 13Jul22", fontsize=24)
    for detector in range(9):
        plotcounter += 1
        plt.subplot(6,3,plotcounter)
        plt.title("%s"%LSSTComCam_camera[detector].getName(), fontsize = 12)
        ptcDataset = butler.get('ptc', exposure=expId, detector=detector)
        if plotData == "Gain":
            data = ptcDataset.gain
            err_data = ptcDataset.gainErr
        else:
            data = ptcDataset.noise
            err_data = ptcDataset.noiseErr
        amps = data.keys()

        yvals = []
        stored_yvals = []
        err = []
        amp_nums = []
        for ii, amp in enumerate(amps):
            yvals.append(data[amp])
            err.append(err_data[amp])
            amp_nums.append(ii)
            amps = data.keys()
            if plotData == "Gain":
                stored_yvals.append(LSSTComCam_camera[detector].getAmplifiers()[ii].getGain())
            else:
                stored_yvals.append(LSSTComCam_camera[detector].getAmplifiers()[ii].getReadNoise())
        plt.errorbar(amp_nums, yvals, yerr=err, label="Meas_20220218")
        plt.plot(amp_nums, stored_yvals, marker='x', label="DM Camera Values")
        if plotData == "Gain":
            plt.ylim(1.0, 2.4)
            plt.ylabel("Gain (e-/ADU)", fontsize = 12)
        else:
            plt.ylim(0.0,50.0)
            plt.ylabel("Read noise (electrons)", fontsize = 12)
        plt.xticks(amp_nums,amps, fontsize=8)
        plt.legend(loc = 'upper right', fontsize = 12)
pdf.savefig(fig2)
```

```python tags=[]
butler = Butler("/repo/main", collections=["LATISS/raw/all", "LATISS/calib/unbounded", \
                                          "u/cslage/latiss/ptc_20210217"])
LATISS_camera = butler.get('camera', instrument='LATISS')
```

```python
expId = 2021021700096
plotcounter = 0
fig3 = plt.figure(figsize=(16,8))
plt.subplots_adjust(hspace=0.3)

for plotData in ["Gain", "Noise"]:
    plt.suptitle("LATISS Gain and Noise - 13Jul22", fontsize=24)
    for detector in range(1):
        plotcounter += 1
        plt.subplot(1,2,plotcounter)
        plt.title("%s"%LATISS_camera[detector].getName(), fontsize = 12)
        ptcDataset = butler.get('ptc', exposure=expId, detector=detector)
        if plotData == "Gain":
            data = ptcDataset.gain
            err_data = ptcDataset.gainErr
        else:
            data = ptcDataset.noise
            err_data = ptcDataset.noiseErr
        amps = data.keys()

        yvals = []
        stored_yvals = []
        err = []
        amp_nums = []
        for ii, amp in enumerate(amps):
            yvals.append(data[amp])
            err.append(err_data[amp])
            amp_nums.append(ii)
            amps = data.keys()
            if plotData == "Gain":
                stored_yvals.append(LATISS_camera[detector].getAmplifiers()[ii].getGain())
            else:
                stored_yvals.append(LATISS_camera[detector].getAmplifiers()[ii].getReadNoise())
        plt.errorbar(amp_nums, yvals, yerr=err, label="Meas_20210217")
        plt.plot(amp_nums, stored_yvals, marker='x', label="DM Camera Values")
        if plotData == "Gain":
            plt.ylim(0.5, 4.0)
            plt.ylabel("Gain (e-/ADU)", fontsize = 12)
        else:
            plt.ylim(0.0,50.0)
            plt.ylabel("Read noise (electrons)", fontsize = 12)
        plt.xticks(amp_nums,amps, fontsize=8)
        plt.legend(loc = 'upper right', fontsize = 12)
pdf.savefig(fig3)
pdf.close()
```

```python
LATISS_camera[0][7].getGain()
```

```python
LATISS_camera[0].getName()
```

```python
LSSTComCam_fits_Filename = '/repo/main/LSSTComCam/calib/DM-28636/unbounded/camera/camera_LSSTComCam_LSSTComCam_calib_DM-28636_unbounded.fits'
```

```python
hdulist = pf.open(LSSTComCam_fits_Filename)
```

```python
hdulist[6].header
```

```python

```
