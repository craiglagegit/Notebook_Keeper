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

# This notebook demonstrates a functionality issue with the select_packed_time_series of the efd-client

```python
import pandas as pd
import numpy as np
from astropy.time import Time, TimeDelta
from lsst_efd_client import EfdClient, resample
```

```python
client = EfdClient('summit_efd')
```

<!-- #region tags=[] -->
### Declare timestamp(s) used for EFD queries
<!-- #endregion -->

```python
# From an observeration last night where a fault occurred
t1 = Time('2022-02-16T05:15:00', scale='utc')
window = TimeDelta(40, format='sec')
t2=t1+window
nas2_single = await client.select_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", \
                                              ["nasmyth2CalculatedAngle0","nasmyth2CalculatedAngle99", \
                                               "private_sndStamp" ,"private_kafkaStamp", "cRIO_timestamp"], \
                                              t1, t2)
nas2_single_index_time=Time(nas2_single.index)

nas2_packed = await client.select_packed_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", \
                                                     ["nasmyth2CalculatedAngle"], t1, t2)
nas2_packed_index_time=Time(nas2_packed.index)
```

```python
(nas2_packed_index_time.utc[0] - nas2_single_index_time.utc[0]).sec
```

```python

```
