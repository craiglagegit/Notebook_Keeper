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

Initially written 18 Feb 2022 by Craig Lage.

```python
import sys, os, glob, time
import numpy as np
import astropy.io.fits as pf
import matplotlib.pyplot as plt
import lsst.daf.butler as dafButler
```

```python
butler1 = dafButler.Butler("/repo/main", collections=["LSSTComCam/raw/all", "LSSTComCam/calib/unbounded",
                                                      "u/cslage/comcam/ptc_20220216"])
exposure1 = 2022021600192
butler2 = dafButler.Butler("/repo/main", collections=["LSSTComCam/raw/all", "LSSTComCam/calib/unbounded",
                                                      "u/cslage/comcam/ptc_20220217"])
exposure2 = 2022021700200
butler3 = dafButler.Butler("/repo/main", collections=["LSSTComCam/raw/all", "LSSTComCam/calib/unbounded",
                                                      "u/cslage/comcam/ptc_20220218"])
exposure3 = 2022021800078
```

```python
RAFT = 'R22'
butlers = [butler1, butler2, butler3]
exposures = [exposure1, exposure2, exposure3]
run_names = ['2022-02-16-Rband', '2022-02-17-Zband', '2022-02-18-Iband']
markers = ['o', 'x', '+', '*', '^', 'v']
```

```python
plt.figure(figsize=(16,16))
plt.subplots_adjust(hspace=0.5)
plt.suptitle("ComCam Gains", fontsize=24)
for i,butler in enumerate(butlers):
    run_name = run_names[i]
    plotcounter = 0
    for detector in range(9):
        plotcounter += 1
        plt.subplot(3,3,plotcounter)
        plt.title("Detector%d"%detector, fontsize = 12)
        ptcDataset = butler.get('ptc', exposure=exposures[i], detector=detector)
        gain_data = ptcDataset.gain
        gain_err_data = ptcDataset.gainErr
        amps = gain_data.keys()
        gains = []
        gain_err = []
        amp_nums = []
        for ii, amp in enumerate(amps):
            gains.append(gain_data[amp])
            gain_err.append(gain_err_data[amp])
            amp_nums.append(ii)
        plt.errorbar(amp_nums, gains, yerr=gain_err, marker = markers[i], label = run_name)
        plt.ylim(1.0, 2.4)
        plt.ylabel("Gain", fontsize = 12)
        plt.xticks(amp_nums,amps, fontsize=8)
        plt.legend(loc = 'upper right', fontsize = 12)
plt.savefig('/repo/main/u/cslage/comcam/ptc_20220217/plots/Gain_Summary_21Feb22.pdf')

```

```python

```
