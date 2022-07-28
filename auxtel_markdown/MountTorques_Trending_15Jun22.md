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

# Using the RubinTV Mount Plotting to identify large errors

Craig Lage - 15-Jun-22

```python
import nest_asyncio
nest_asyncio.apply()
import sys, time, os, asyncio
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from lsst.rubintv.production import mountTorques

from lsst_efd_client import EfdClient
from lsst.daf.butler import Butler
import lsst.log as log
```

```python
client = EfdClient('ldf_stable_efd')
butler = Butler('/repo/main', collections="LATISS/raw/all")
logger = log.getLogger('myLogger')
```

```python
figure = plt.figure(figsize=(16,22))
saveFilename = None

errs = []
for dayObs in [20220609, 20220503]:
    exposureList = []

    for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
        if record.observation_type not in ['bias', 'flat', 'dark']:
            exposureList.append(record.id)
    exposureList = sorted(exposureList, key=lambda x: x[0])
    for expId in exposureList:
        saveFilename = None
        dataId = {'detector':0, 'exposure':expId}
        try:
            err = mountTorques.plotMountTracking(dataId, butler, client, figure, saveFilename, logger)
            print(expId, err[0], err[1], err[2])
            errs.append([expId, err[0], err[1], err[2]])
            if err[0] > 0.25 or err[1] > 0.25:
                print(f"Plotting expId {expId}")
                saveFilename = f'/project/cslage/AuxTel/mount_graphs/large_errors_07jul22/Mount_Torques_{expId}.pdf'
                err = mountTorques.plotMountTracking(dataId, butler, client, figure, saveFilename, logger)
        except:
            continue
```

I converted this to a Python script, and ran it against 20220503, 20220504, 20220505, 20220628, 20220629, and 20220630.  I then manually sorted through the plots (which are in AuxTel/mount_graphs/large_errors_07jul22) and categorized them in a spreadsheet.  The results are as follows:

```python
import numpy as np
import matplotlib.pyplot as plt
```

```python
categories = ["Not Stable", "Jitter", "Slew during exp", "Crazy mount", "Oscillation", "Other"]
nights = ["20220503", "20220504", "20220505", "20220628", "20220629", "20220630"]
windSpeed = [10.8, 7.9, 3.0, 6.4, 8.2, 3.9]
xplot1 = np.array(list(range(len(nights))))
colors = ['red', 'blue', 'green', 'orange', 'cyan', 'magenta']
failures = {}
failures["20220503"] = [0.011396,0.142450,0.002849,0.000000,0.000000,0.000000]
failures["20220504"] = [0.045113,0.005013,0.005013,0.000000,0.000000,0.010025]
failures["20220505"] = [0.014474,0.000000,0.000000,0.000000,0.007895,0.000000]
failures["20220628"] = [0.001464,0.046852,0.004392,0.001464,0.001464,0.000000]
failures["20220629"] = [0.003802,0.005703,0.001901,0.001901,0.011407,0.000000]
failures["20220630"] = [0.006443,0.002577,0.001289,0.003866,0.001289,0.001289]
```

```python
fig = plt.figure(figsize=(11,8))
ax = plt.subplot(1,1,1)
shift = 0.1
for j in range(len(categories)):
    xplot = []
    barplot = []
    for i in range(len(nights)):
        xplot.append(i + j * shift)
        barplot.append(failures[nights[i]][j])
    ax.bar(xplot, barplot, width = shift, color = colors[j], label = categories[j])
ax.set_xticks(xplot1)
ax.set_xticklabels(nights)
ax.legend()
ax.set_ylabel("Fraction of exposures with RMS>0.25 arcseconds")

ax2 = ax.twinx()
ax2.plot(xplot1, windSpeed, marker = 'x', label = "Wind Speed")
ax2.legend(loc = 'upper center')
ax2.set_ylabel("Wind speed (m/s)")
plt.savefig("/project/cslage/AuxTel/mount_graphs/large_errors_07jul22/Fail_Categories_07Jul22.pdf")
```

```python

```
