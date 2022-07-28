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

## Running ISR from a notebook

Craig Lage - 16-Jun-22


    
        Here's an example of running ISR with a pipetask using the yaml file:

          isr:
            class: lsst.ip.isr.IsrTask
            config:
              connections.ccdExposure: raw
              connections.outputExposure: parameters.exposureName
              doWrite: true
              doOverscan: true
              doAssembleCcd: true
              doBias: true
              doVariance: true
              doLinearize: false
              doCrosstalk: false
              doBrighterFatter: false
              doDark: true
              doStrayLight: false
              doFlat: false
              doFringe: false
              doApplyGains: false
              doDefect: true
              doNanMasking: true
              doInterpolate: false
              doSaturation: false
              doSaturationInterpolation: false
              growSaturationFootprintSize: 0
     

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.ip.isr import IsrTask, IsrTaskConfig
from lsst.daf.butler import Butler
```

```python
butler = Butler('/repo/LATISS', collections=["LATISS/raw/all", "LATISS/calib"])
```

```python
isrConfig = IsrTaskConfig()
isrConfig.doLinearize=False
isrConfig.doOverscan=True
isrConfig.doAssembleCcd=True
isrConfig.doBias=True
isrConfig.doVariance=True
isrConfig.doLinearize=False
isrConfig.doCrosstalk=False
isrConfig.doBrighterFatter=False
isrConfig.doDark=False
isrConfig.doStrayLight=False
isrConfig.doFlat=False
isrConfig.doFringe=False
isrConfig.doApplyGains=False
isrConfig.doDefect=False
isrConfig.doNanMasking=True
isrConfig.doInterpolate=False
isrConfig.doSaturation=False
isrConfig.doSaturationInterpolation=False


# Adjust these as needed and add as many more as you want
```

```python
isrTask = IsrTask(config=isrConfig)
```

```python
expId = 2022060900070
exp = butler.get('raw', detector=0, exposure=expId)
biasExp = butler.get('bias', detector=0, exposure=expId)
```

```python
isrResult = isrTask.run(exp, bias=biasExp)
```

```python
# Now look at the data with matplotlib
from matplotlib.colors import LogNorm

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
arr = isrResult.exposure.image.array
img = plt.imshow(arr, norm=LogNorm(vmin=10, vmax=1000), interpolation='Nearest', cmap='gray')
colorbar(img)
plt.tight_layout(h_pad=1)

```

```python

```
