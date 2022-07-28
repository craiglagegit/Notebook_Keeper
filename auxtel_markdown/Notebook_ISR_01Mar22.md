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

## ISR - Gen 3

In this notebook, we test running ISR with the Gen3 butler\
Craig Lage - 01-Mar-22

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import lsst.afw.cameraGeom.utils as camGeomUtils
import lsst.afw.display as afwDisplay
from lsst.ip.isr.isrTask import IsrTask, IsrTaskConfig
from lsst.daf.butler import Butler
```

```python
butler = Butler('/repo/main', collections="LATISS/raw/all")
```

```python jupyter={"outputs_hidden": true} tags=[]
# Gen3 butler
dayObs = '2022-02-16'
dayObs = int(dayObs.replace('-', ''))

exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
    exposureList.append([record.id, record])
exposureList.sort(key=lambda x: x[0])
for [id,record] in exposureList:
    print(record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name)

```

```python

```

```python
# Config for processing raw bias image
isrConfig = IsrTaskConfig()
isrConfig.doLinearize = False
isrConfig.doBias = False
isrConfig.doFlat = False
isrConfig.doDark = False
isrConfig.doFringe = False
isrConfig.doDefect = False
isrConfig.doWrite = False
isrConfig.doBrighterFatter = False
isrTask = IsrTask(config=isrConfig)
```

```python
biasExpId = 2022021600015
rawBiasImage = butler.get('raw', detector=0, exposure=biasExpId)
```

```python
biasPostISR = isrTask.run(rawBiasImage).exposure
```

```python
flatExpId=2022021600060
rawFlatImage = butler.get('raw', detector=0, exposure=flatExpId)

```

```python
isrConfig.doBias = True
#isrConfig.doApplyGains = True
isrTask = IsrTask(config=isrConfig)
```

```python tags=[]
flatPostISR = isrTask.run(rawFlatImage, bias=biasPostISR).exposure
```

```python
arrayData = flatPostISR.image.array
mean = arrayData.mean()
std = arrayData.std()
print(mean, std)
```

```python
# Look at the data with matplotlib
# The raw data doesn't look very good, because of the large pedestal of about 15,000 ADU
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

plt.figure(figsize=(8,8))
plt.suptitle(f"Image",fontsize=18)
arr = flatPostISR.image.array
img = plt.imshow(arr, norm=LogNorm(vmin=mean-3.0*std, vmax=mean+3.0*std), interpolation='Nearest', cmap='gray')
colorbar(img)
plt.tight_layout(h_pad=1)
#plt.savefig(REPO_DIR+"/plots/NGC4755_17Feb21.pdf")
```

```python

```
