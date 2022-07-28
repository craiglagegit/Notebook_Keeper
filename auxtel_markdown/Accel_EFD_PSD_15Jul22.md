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

# AuxTel accelerometer PSD - 15-Jul-22
Testing the unpacking code for lsst-efd-client \
Craig Lage

```python
import sys, time, os, asyncio
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
from astropy.time import Time, TimeDelta
from lsst.daf.butler import Butler
from lsst_efd_client.efd_helper import EfdClient
from lsst_efd_client.efd_utils import merge_packed_PSD
```

```python
# Get EFD client
client = EfdClient('ldf_stable_efd')
```

```python
# Times to start looking at PSD data
start = Time("2022-07-14 12:00:00Z", scale='utc')
end = Time("2022-07-14 12:00:10Z", scale='utc')
```

# First, test the merge_packed_PSD function.

```python
accel_data = await client.select_time_series("lsst.sal.ESS.accelerometerPSD", ["*"], start, end)

indexCounter = 0
plt.figure(figsize=(16,16))
plt.subplots_adjust(wspace=0.5, hspace=0.3)
axes = ['X', 'Y', 'Z']
sensors = ["AuxTel-M1", "AuxTel-M2", "AuxTel-Truss"]
plotCounter = 1
for sensor in sensors:
    for axis in axes:
        base_field = f"accelerationPSD{axis}"
        plt.subplot(3,3,plotCounter)
        plt.title(f"{sensor} - {axis}", fontsize=12)
        df = merge_packed_PSD(accel_data, base_field, sensor)
        row = df.iloc[indexCounter][2:]
        row.plot()
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('PSD [m^2/(Hz s^4)]')
        plotCounter += 1
timestamp = df.index[0].strftime("%Y%m%dT%H%M%SZ")
plt.suptitle(f"Accelerometer Power Spectral Density - {timestamp}", fontsize=16)
plt.savefig(f"/project/cslage/AuxTel/accel_data/Accel_PSD_{timestamp}.pdf")
```

# Next, test the select_packed_PSD function with a single axis and sensor.

```python
indexCounter = 0
plt.figure(figsize=(16,16))
plt.subplots_adjust(wspace=0.5, hspace=0.3)
xes = ['X', 'Y', 'Z']
sensors = ["AuxTel-M1", "AuxTel-M2", "AuxTel-Truss"]
plotCounter = 1
for sensor in sensors:
    for axis in axes:
        base_field = f"accelerationPSD{axis}"
        plt.subplot(3,3,plotCounter)
        plt.title(f"{sensor} - {axis}", fontsize=12)
        df = await client.select_packed_PSD("lsst.sal.ESS.accelerometerPSD", base_field, sensor, start, end)
        row = df.iloc[indexCounter][2:]
        row.plot()
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('PSD [m^2/(Hz s^4)]')
        plotCounter += 1
timestamp = df.index[0].strftime("%Y%m%dT%H%M%SZ")
plt.suptitle(f"Accelerometer Power Spectral Density - {timestamp}", fontsize=16)
plt.savefig(f"/project/cslage/AuxTel/accel_data/Accel_PSD_{timestamp}.pdf")
```

# Next, test the select_packed_PSD function with multiple axes and sensors.

```python
indexCounter = 0
plt.figure(figsize=(16,16))
plt.subplots_adjust(wspace=0.5, hspace=0.3)
axes = ['X', 'Y', 'Z']
sensors = ["AuxTel-M1", "AuxTel-M2", "AuxTel-Truss"]
base_fields = []
for axis in axes:
    base_fields.append(f"accelerationPSD{axis}")
df = await client.select_packed_PSD("lsst.sal.ESS.accelerometerPSD", base_fields, sensors, start, end)
plotCounter = 1
for sensor in sensors:
    for axis in axes:
        SensorName = f"{sensor}-{axis}"
        plt.subplot(3,3,plotCounter)
        plt.title(f"{sensor} - {axis}", fontsize=12)
        plot_df = df[df["SensorName"] == SensorName]
        row = plot_df.iloc[indexCounter][2:]
        row.plot()
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('PSD [m^2/(Hz s^4)]')
        plotCounter += 1
timestamp = df.index[0].strftime("%Y%m%dT%H%M%SZ")
plt.suptitle(f"Accelerometer Power Spectral Density - {timestamp}", fontsize=16)
plt.savefig(f"/project/cslage/AuxTel/accel_data/Accel_PSD_{timestamp}.pdf")
```

# This shows the structure of the dataframe after unpacking.

```python
df.head(9)
```

```python

```
