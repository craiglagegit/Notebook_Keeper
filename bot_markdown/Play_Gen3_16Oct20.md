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

# Notebook for studying Gen3 butler

Initially written 20 Oct 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
! eups list -s | grep cp_pipe
```

```python
import sys, os, glob, time
import numpy as np
import matplotlib.pyplot as plt
from lsst.daf.butler import Butler
```

```python
REPO_DIR = '/lsstdata/offline/teststand/BOT/gen3repo'
butler = Butler(REPO_DIR)
```

```python
list(butler.registry.queryCollections())
```

```python
print(dir(butler))
print(dir(butler.registry))
```

```python
list(butler.registry.queryDatasetTypes())
```

```python jupyter={"outputs_hidden": true}
list(butler.registry.queryDatasets(datasetType='raw', collections=['LSSTCam/raw/all'], dataId={'detector':47, 'instrument':'LSSTCam'}))
```

```python
test = butler.get('raw', collections=['LSSTCam/raw/all'], dataId={'detector':47,'instrument':'LSSTCam','exposure':3020090100162})
```

```python
print(type(test))
print(dir(test))
print(test.getMetadata())
```

```python
arr = test.image.array
plt.imshow(arr)
```

```python

```
