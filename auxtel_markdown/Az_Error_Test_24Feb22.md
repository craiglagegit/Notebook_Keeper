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

# AuxTel Plot tracking - 22-Feb-22

In this notebook, investigate impact of fan on mount tracking.

```python
import sys, time, os, asyncio

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
%matplotlib widget
import pandas as pd
from astropy.time import Time, TimeDelta
from lsst_efd_client import EfdClient
```

```python
# Get EFD client
client = EfdClient('summit_efd')
```

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-02-24T17:40:00", scale='utc')
tend = Time("2022-02-24T17:42:00", scale='utc')
az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "azimuthCalculatedAngle",  tstart, tend)
el = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "elevationCalculatedAngle",  tstart, tend)
```

```python
# Calculate the tracking errors
az_vals = np.array(az.values.tolist())[:,0]
el_vals = np.array(el.values.tolist())[:,0]
times = np.array(az.values.tolist())[:,1]
times = times - times [0]

# Fit with a quartic
az_fit = np.polyfit(times, az_vals, 4)
el_fit = np.polyfit(times, el_vals, 4)

az_model = az_fit[0] * times * times * times * times + az_fit[1] * times * times * times \
    + az_fit[2] * times *times + az_fit[3] * times + az_fit[4]
el_model = el_fit[0] * times * times * times * times + el_fit[1] * times * times * times \
    + el_fit[2] * times * times + el_fit[3] * times + el_fit[4]

# Errors in arcseconds
az_error = (az_vals - az_model) * 3600
el_error = (el_vals - el_model) * 3600

# Calculate RMS
az_rms = np.sqrt(np.mean(az_error * az_error))
el_rms = np.sqrt(np.mean(el_error * el_error))

```

```python
# Plot it
fig = plt.figure(figsize = (8,8))
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.suptitle(f"Mount Tracking Az Error 20220224\nStarting Az=90, Starting El=45", fontsize = 18)
plt.subplot(2,2,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.set_ylabel("Degrees")
plt.subplot(2,2,2)
ax3 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
ax3.set_ylabel("Degrees")

plt.subplot(2,2,3)
plt.plot(times, az_error, color='red')
plt.title(f"Azimuth RMS error = {az_rms:.2f} arcseconds")
plt.ylim(-10.0,10.0)
plt.xticks([])
plt.ylabel("ArcSeconds")
plt.subplot(2,2,4)
plt.plot(times, el_error, color='green')
plt.title(f"Elevation RMS error = {el_rms:.2f} arcseconds")
plt.ylim(-10.0,10.0)
plt.xticks([])
plt.ylabel("ArcSeconds")
plt.savefig(f"/project/cslage/AuxTel/mount_graphs/Az_Error_Debug_1_24Feb22.pdf")

```

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-02-24T17:48:00", scale='utc')
tend = Time("2022-02-24T17:50:00", scale='utc')
az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "azimuthCalculatedAngle",  tstart, tend)
el = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "elevationCalculatedAngle",  tstart, tend)
```

```python
# Calculate the tracking errors
az_vals = np.array(az.values.tolist())[:,0]
el_vals = np.array(el.values.tolist())[:,0]
times = np.array(az.values.tolist())[:,1]
times = times - times [0]

# Fit with a quartic
az_fit = np.polyfit(times, az_vals, 4)
el_fit = np.polyfit(times, el_vals, 4)

az_model = az_fit[0] * times * times * times * times + az_fit[1] * times * times * times \
    + az_fit[2] * times *times + az_fit[3] * times + az_fit[4]
el_model = el_fit[0] * times * times * times * times + el_fit[1] * times * times * times \
    + el_fit[2] * times * times + el_fit[3] * times + el_fit[4]

# Errors in arcseconds
az_error = (az_vals - az_model) * 3600
el_error = (el_vals - el_model) * 3600

# Calculate RMS
az_rms = np.sqrt(np.mean(az_error * az_error))
el_rms = np.sqrt(np.mean(el_error * el_error))

```

```python
# Plot it
fig = plt.figure(figsize = (8,8))
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.suptitle(f"Mount Tracking Az Error 20220224\nStarting Az=0, Starting El=80", fontsize = 18)
plt.subplot(2,2,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.set_ylabel("Degrees")
plt.subplot(2,2,2)
ax3 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
ax3.set_ylabel("Degrees")

plt.subplot(2,2,3)
plt.plot(times, az_error, color='red')
plt.title(f"Azimuth RMS error = {az_rms:.2f} arcseconds")
plt.ylim(-10.0,10.0)
plt.xticks([])
plt.ylabel("ArcSeconds")
plt.subplot(2,2,4)
plt.plot(times, el_error, color='green')
plt.title(f"Elevation RMS error = {el_rms:.2f} arcseconds")
plt.ylim(-10.0,10.0)
plt.xticks([])
plt.ylabel("ArcSeconds")
plt.savefig(f"/project/cslage/AuxTel/mount_graphs/Az_Error_Debug_2_24Feb22.pdf")

```

```python

```
