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
client = EfdClient('ldf_stable_efd')
```

```python tags=[]
# Now get the mount tracking info for the time of the test
tstart = Time("2022-02-24T17:06:00", scale='utc')
tend = Time("2022-02-24T17:20:00", scale='utc')
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
# Fan powers
p0Power = Time("2022-02-24T17:08:00", scale='utc')
p25Power = Time("2022-02-24T17:10:00", scale='utc')
p50Power = Time("2022-02-24T17:12:00", scale='utc')
p75Power = Time("2022-02-24T17:14:00", scale='utc')
fullPower = Time("2022-02-24T17:16:00", scale='utc')
zeroPower = Time("2022-02-24T17:18:00", scale='utc')
```

```python
# Plot it
fig = plt.figure(figsize = (8,8))
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.suptitle(f"Mount Tracking vs Fan speed 20220224", fontsize = 18)
plt.subplot(2,2,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.axvline(p0Power.to_datetime(), color="blue", ls = '--')
ax1.text(p0Power.to_datetime(), 0.5, 'On-0%', rotation=90)
ax1.axvline(p25Power.to_datetime(), color="blue", ls = '--')
ax1.text(p25Power.to_datetime(), 0.5, '25%', rotation=90)
ax1.axvline(p50Power.to_datetime(), color="cyan", ls = '--')
ax1.text(p50Power.to_datetime(), 0.5, '50%', rotation=90)
ax1.axvline(p75Power.to_datetime(), color="magenta", ls = '--')
ax1.text(p75Power.to_datetime(), 0.5, '75%', rotation=90)
ax1.axvline(fullPower.to_datetime(), color="black", ls = '--')
ax1.text(fullPower.to_datetime(), 0.5, '100%', rotation=90)
ax1.axvline(zeroPower.to_datetime(), color="green", ls = '--')
ax1.text(zeroPower.to_datetime(), 0.5, 'Off', rotation=90)
ax1.set_ylabel("Degrees")
plt.subplot(2,2,2)
ax3 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
ax3.axvline(p0Power.to_datetime(), color="blue", ls = '--')
ax3.text(p0Power.to_datetime(), 0.5, 'On-0%', rotation=90)
ax3.axvline(p25Power.to_datetime(), color="blue", ls = '--')
ax3.text(p25Power.to_datetime(), 80.11, '25%', rotation=90)
ax3.axvline(p50Power.to_datetime(), color="cyan", ls = '--')
ax3.text(p50Power.to_datetime(), 80.11, '50%', rotation=90)
ax3.axvline(p75Power.to_datetime(), color="magenta", ls = '--')
ax3.text(p75Power.to_datetime(), 80.11, '75%', rotation=90)
ax3.axvline(fullPower.to_datetime(), color="black", ls = '--')
ax3.text(fullPower.to_datetime(), 80.11, '100%', rotation=90)
ax3.axvline(zeroPower.to_datetime(), color="green", ls = '--')
ax3.text(zeroPower.to_datetime(), 80.11, '0%', rotation=90)
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
plt.savefig(f"/project/cslage/AuxTel/mount_graphs/Mount_Errors_Fan_0_80_24Feb22.pdf")

```

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-02-24T17:06:00", scale='utc')
tend = Time("2022-02-24T17:08:00", scale='utc')
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
plt.suptitle(f"Mount Tracking vs Fan speed 20220224", fontsize = 18)
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
plt.savefig(f"/project/cslage/AuxTel/mount_graphs/Mount_Errors_Fan_BlowUp_0_80_24Feb22.pdf")

```

```python tags=[]
# Running at 90/45
# Now get the mount tracking info for the time of the test
tstart = Time("2022-02-24T17:57:00", scale='utc')
tend = Time("2022-02-24T18:11:00", scale='utc')
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
# Fan powers
p0Power = Time("2022-02-24T17:59:00", scale='utc')
p25Power = Time("2022-02-24T18:01:00", scale='utc')
p50Power = Time("2022-02-24T18:03:00", scale='utc')
p75Power = Time("2022-02-24T18:05:00", scale='utc')
fullPower = Time("2022-02-24T18:07:00", scale='utc')
zeroPower = Time("2022-02-24T18:09:00", scale='utc')
```

```python
# Plot it
fig = plt.figure(figsize = (8,8))
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.suptitle(f"Mount Tracking vs Fan speed 20220224", fontsize = 18)
plt.subplot(2,2,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.axvline(p0Power.to_datetime(), color="blue", ls = '--')
ax1.text(p0Power.to_datetime(), 90.0, 'On-0%', rotation=90)
ax1.axvline(p25Power.to_datetime(), color="blue", ls = '--')
ax1.text(p25Power.to_datetime(), 90.0, '25%', rotation=90)
ax1.axvline(p50Power.to_datetime(), color="cyan", ls = '--')
ax1.text(p50Power.to_datetime(), 90.0, '50%', rotation=90)
ax1.axvline(p75Power.to_datetime(), color="magenta", ls = '--')
ax1.text(p75Power.to_datetime(), 90.0, '75%', rotation=90)
ax1.axvline(fullPower.to_datetime(), color="black", ls = '--')
ax1.text(fullPower.to_datetime(), 90.0, '100%', rotation=90)
ax1.axvline(zeroPower.to_datetime(), color="green", ls = '--')
ax1.text(zeroPower.to_datetime(), 90.0, 'Off', rotation=90)
ax1.set_ylabel("Degrees")
plt.subplot(2,2,2)
ax3 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
ax3.axvline(p0Power.to_datetime(), color="blue", ls = '--')
ax3.text(p0Power.to_datetime(), 48.7, 'On-0%', rotation=90)
ax3.axvline(p25Power.to_datetime(), color="blue", ls = '--')
ax3.text(p25Power.to_datetime(), 48.7, '25%', rotation=90)
ax3.axvline(p50Power.to_datetime(), color="cyan", ls = '--')
ax3.text(p50Power.to_datetime(), 48.7, '50%', rotation=90)
ax3.axvline(p75Power.to_datetime(), color="magenta", ls = '--')
ax3.text(p75Power.to_datetime(), 80.11, '75%', rotation=90)
ax3.axvline(fullPower.to_datetime(), color="black", ls = '--')
ax3.text(fullPower.to_datetime(), 80.11, '100%', rotation=90)
ax3.axvline(zeroPower.to_datetime(), color="green", ls = '--')
ax3.text(zeroPower.to_datetime(), 80.11, '0%', rotation=90)
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
plt.savefig(f"/project/cslage/AuxTel/mount_graphs/Mount_Errors_Fan_90_45_24Feb22.pdf")

```

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-02-24T17:57:00", scale='utc')
tend = Time("2022-02-24T17:59:00", scale='utc')
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
plt.suptitle(f"Mount Tracking vs Fan speed 20220224", fontsize = 18)
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
plt.savefig(f"/project/cslage/AuxTel/mount_graphs/Mount_Errors_Fan_BlowUp_90_45_24Feb22.pdf")

```

```python tags=[]
# Running at 90/45
# Slit open ~50%
# Now get the mount tracking info for the time of the test
tstart = Time("2022-02-24T18:22:00", scale='utc')
tend = Time("2022-02-24T18:36:00", scale='utc')
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
# Fan powers
p0Power = Time("2022-02-24T18:24:00", scale='utc')
p25Power = Time("2022-02-24T18:26:00", scale='utc')
p50Power = Time("2022-02-24T18:28:00", scale='utc')
p75Power = Time("2022-02-24T18:30:00", scale='utc')
fullPower = Time("2022-02-24T18:32:00", scale='utc')
zeroPower = Time("2022-02-24T18:34:00", scale='utc')
```

```python
# Plot it
fig = plt.figure(figsize = (8,8))
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.suptitle(f"Mount Tracking vs Fan speed 20220224\n Dome open", fontsize = 18)
plt.subplot(2,2,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.axvline(p0Power.to_datetime(), color="blue", ls = '--')
ax1.text(p0Power.to_datetime(), 89.9, 'On-0%', rotation=90)
ax1.axvline(p25Power.to_datetime(), color="blue", ls = '--')
ax1.text(p25Power.to_datetime(), 89.9, '25%', rotation=90)
ax1.axvline(p50Power.to_datetime(), color="cyan", ls = '--')
ax1.text(p50Power.to_datetime(), 89.9, '50%', rotation=90)
ax1.axvline(p75Power.to_datetime(), color="magenta", ls = '--')
ax1.text(p75Power.to_datetime(), 89.9, '75%', rotation=90)
ax1.axvline(fullPower.to_datetime(), color="black", ls = '--')
ax1.text(fullPower.to_datetime(), 89.9, '100%', rotation=90)
ax1.axvline(zeroPower.to_datetime(), color="green", ls = '--')
ax1.text(zeroPower.to_datetime(), 89.9, 'Off', rotation=90)
ax1.set_ylabel("Degrees")
plt.subplot(2,2,2)
ax3 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
ax3.axvline(p0Power.to_datetime(), color="blue", ls = '--')
ax3.text(p0Power.to_datetime(), 48.7, 'On-0%', rotation=90)
ax3.axvline(p25Power.to_datetime(), color="blue", ls = '--')
ax3.text(p25Power.to_datetime(), 48.7, '25%', rotation=90)
ax3.axvline(p50Power.to_datetime(), color="cyan", ls = '--')
ax3.text(p50Power.to_datetime(), 48.7, '50%', rotation=90)
ax3.axvline(p75Power.to_datetime(), color="magenta", ls = '--')
ax3.text(p75Power.to_datetime(), 48.7, '75%', rotation=90)
ax3.axvline(fullPower.to_datetime(), color="black", ls = '--')
ax3.text(fullPower.to_datetime(), 48.7, '100%', rotation=90)
ax3.axvline(zeroPower.to_datetime(), color="green", ls = '--')
ax3.text(zeroPower.to_datetime(), 48.7, '0%', rotation=90)
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
plt.savefig(f"/project/cslage/AuxTel/mount_graphs/Mount_Errors_Fan_Dome_90_45_24Feb22.pdf")

```

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-02-24T18:22:00", scale='utc')
tend = Time("2022-02-24T18:36:00", scale='utc')
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
plt.suptitle(f"Mount Tracking vs Fan speed 20220224\n Dome open", fontsize = 18)
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
plt.savefig(f"/project/cslage/AuxTel/mount_graphs/Mount_Errors_Fan_Dome_BlowUp_90_45_24Feb22.pdf")

```

```python

```
