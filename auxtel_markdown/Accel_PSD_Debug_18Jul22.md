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

# AuxTel accelerometer PSD - 18-Jul-22
Comparing PSD plots \
Craig Lage

```python
import sys, time, os, asyncio
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import pandas as pd
import pickle as pkl
from astropy.time import Time, TimeDelta
from lsst.daf.butler import Butler
from lsst_efd_client.efd_helper import EfdClient
from lsst_efd_client.efd_utils import merge_packed_PSD
```

```python
# First plot data with the notebook
# This is stored data from before the CSC implementation
filename = '/project/cslage/AuxTel/mount_graphs/20220504T033329Z.pkl'
# Unpickle the accel dataframe
file = open(filename, 'rb')
df = pkl.load(file)
file.close()

startingIndex = 5000
endingIndex = startingIndex + 400
# Grab a 2 second subset of the data
subdf = df.iloc[startingIndex:endingIndex]
timestamp = subdf.index[0].strftime("%Y%m%dT%H%M%SZ")
num_samples = len(subdf)
sampling_frequency = 200.0
sampling_interval = 1 / sampling_frequency
psd_frequencies = np.fft.rfftfreq(num_samples, sampling_interval)

plt.figure(figsize=(16,16))
plt.subplots_adjust(wspace=0.5, hspace=0.7)
plt.suptitle(f"Accelerometer Power Spectral Density \n \n Notebook - {timestamp}", fontsize=24)

axes = ['AZ', 'EL', 'Z']
sensors = ["M1", "M2", "T"]
plotCounter = 1
for sensor in sensors:
    for axis in axes:
        SensorName = f"{axis}{sensor}"
        scaled_data = subdf[SensorName].to_list()
        scaled_data = np.array(scaled_data) * 9.8 # Convert to m/s^2
        psd = np.abs(np.fft.rfft(scaled_data)) ** 2
        plt.subplot(7,3,plotCounter)
        plt.plot(psd_frequencies[1:-2], psd[1:-2], color='blue')
        plt.title(f"{SensorName}", fontsize=12)
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('PSD [m^2/(Hz s^4)]')
        plotCounter += 1
        
# Now plot CSC data        
# Get EFD client
client = EfdClient('ldf_stable_efd') 

# Times to start looking at PSD data
start = Time("2022-07-14 12:00:00Z", scale='utc')
end = Time("2022-07-14 12:00:10Z", scale='utc') 

indexCounter = 0
#plt.figure(figsize=(16,16))
#plt.subplots_adjust(wspace=0.5, hspace=0.3)
axes = ['X', 'Y', 'Z']
sensors = ["AuxTel-M1", "AuxTel-M2", "AuxTel-Truss"]
base_fields = []
for axis in axes:
    base_fields.append(f"accelerationPSD{axis}")
df = await client.select_packed_PSD("lsst.sal.ESS.accelerometerPSD", base_fields, sensors, start, end)
plotCounter += 3
for sensor in sensors:
    for axis in axes:
        SensorName = f"{sensor}-{axis}"
        plt.subplot(7,3,plotCounter)
        plt.title(f"{sensor} - {axis}", fontsize=12)
        plot_df = df[df["SensorName"] == SensorName]
        row = plot_df.iloc[indexCounter][2:]
        row.plot(color='red')
        plt.xlabel('Frequency [Hz]')
        plt.ylabel('PSD [m^2/(Hz s^4)]')
        plotCounter += 1

plt.subplot(7,3,11, frame_on = False)
plt.axis('off')
timestamp = df.index[0].strftime("%Y%m%dT%H%M%SZ")
plt.text(-.2, 0, f"CSC - {timestamp}", fontsize=24)
plt.savefig(f"/project/cslage/AuxTel/accel_data/Accel_PSD_{timestamp}.pdf")
plt.savefig(f"/project/cslage/AuxTel/accel_data/Accel_PSD_Comparison_18Jul22.pdf")
```

```python

```
