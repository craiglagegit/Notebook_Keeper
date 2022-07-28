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

## Timebase testing - NCSA
Craig Lage - Mar 8, 2022

```python
import pandas as pd
from astropy.time import Time, TimeDelta
from lsst_efd_client import EfdClient
```

```python
client = EfdClient('summit_efd')
```

```python
tstart = Time("2022-03-08T14:12:58.283Z", scale='utc')
tend = Time("2022-03-08T14:13:03.523Z", scale='utc')

print(f"tstart={tstart}, tend={tend}")

az = await client.select_packed_time_series('lsst.sal.ATMCS.mount_AzEl_Encoders', \
                                            'azimuthCalculatedAngle', tstart, tend)
print(f"az_start={az.index[0]}, az_end={az.index[-1]}")
offset = (tstart.jd - az.index[0].to_julian_date()) * 86400.0
print(f"Time Offset = {offset:.2f} seconds") 
```

```python
tstart = Time("2022-03-08T14:12:58.283Z", scale='utc')
tend = Time("2022-03-08T14:13:03.523Z", scale='utc')

print(f"tstart={tstart}, tend={tend}")

az = await client.select_packed_time_series('lsst.sal.ATMCS.mount_AzEl_Encoders', \
                                            'azimuthCalculatedAngle', tstart, tend)
az0 = await client.select_time_series('lsst.sal.ATMCS.mount_AzEl_Encoders', \
                                            'azimuthCalculatedAngle0', tstart, tend)
crio = await client.select_time_series('lsst.sal.ATMCS.mount_AzEl_Encoders', \
                                            'cRIO_timestamp', tstart, tend)
print(f"az_start={az.index[0]}, az_end={az.index[-1]}")
print(f"az0_start={az0.index[0]}, az0_end={az0.index[-1]}")
print(f"crio_start={crio.index[0]}, crio_end={crio.index[-1]}")
offset1 = (tstart.jd - az.index[0].to_julian_date()) * 86400.0
offset2 = (tstart.jd - az0.index[0].to_julian_date()) * 86400.0
offset3 = (tstart.jd - crio.index[0].to_julian_date()) * 86400.0
offset4 = (az0.index[0].to_julian_date() - az.index[0].to_julian_date()) * 86400.0
offset5 = (az0.index[0].to_julian_date() - crio.index[0].to_julian_date()) * 86400.0

print(f"tstart-az[0] = {offset1:.2f} seconds") 
print(f"tstart-az0[0] = {offset2:.2f} seconds") 
print(f"tstart-crio[0] = {offset3:.2f} seconds") 
print(f"az0[0]-az[0] = {offset4:.2f} seconds") 
print(f"az0[0]-crio[0] = {offset5:.2f} seconds") 
```

```python

```
