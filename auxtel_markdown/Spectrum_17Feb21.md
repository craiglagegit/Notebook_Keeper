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

```python
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler
from lsst.rapid.analysis.summarizeImage import SummarizeImage
%matplotlib inline
```

```python
butler = Butler('/project/shared/auxTel/rerun/quickLook/')
```

```python
exp = butler.get('quickLookExp', dataId={'expId': 2021021700332, 'detector': 0})
```

```python
summarise = SummarizeImage(exp)
summarise.run()
```

```python

```

```python
REPO_DIR = '/readonly/lsstdata/auxtel/base/auxtel/oods/butler/repo'
butler = Butler(REPO_DIR)
dayObs = '2021-02-17'
```

```python
visits = butler.queryMetadata('raw', ['expId', 'EXPTIME', 'imagetype', 'filter', 'DATE'],\
                              detector=0, dayObs=dayObs)
```

```python
for visit in visits:
    print(visit)
```
