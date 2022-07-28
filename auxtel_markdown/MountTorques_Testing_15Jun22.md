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

# Testing the RubinTV Mount Plotting

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
expId = 2022052600012
dataId = {'detector':0, 'exposure':expId}
figure = plt.figure(figsize=(16,22))
saveFilename = f'/project/cslage/AuxTel/mount_graphs/Mount_Torques_{expId}.pdf'
mountTorques.plotMountTracking(dataId, butler, client, figure, saveFilename, logger)
```

```python
expId = 2022050300244
dataId = {'detector':0, 'exposure':expId}
figure = plt.figure(figsize=(16,22))
saveFilename = f'/project/cslage/AuxTel/mount_graphs/Mount_Torques_{expId}.pdf'
mountTorques.plotMountTracking(dataId, butler, client, figure, saveFilename, logger)
```

```python

```
