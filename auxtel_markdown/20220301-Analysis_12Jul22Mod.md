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

## Motion Tests - 8s azimuth oscillation investigation



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
# for tab completion to work in current notebook instance
%config IPCompleter.use_jedi = False
```

```python
import logging
stream_handler = logging.StreamHandler(sys.stdout)
logger = logging.getLogger()
logger.addHandler(stream_handler)
logger.level = logging.DEBUG
```

```python jupyter={"outputs_hidden": true} tags=[]
# Get EFD client and bring in Lupton's unpacking code
client = EfdClient('summit_efd')
```

# EL 80

```python jupyter={"outputs_hidden": true} tags=[]
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T19:27:00", scale='utc')
tend = Time("2022-03-01T19:29:00", scale='utc')
az = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "azimuthCalculatedAngle",  tstart, tend)
el = await client.select_packed_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", "elevationCalculatedAngle",  tstart, tend)
rot = await client.select_packed_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", 'nasmyth2CalculatedAngle', tstart, tend)
az_torque_1 = await client.select_packed_time_series("lsst.sal.ATMCS.measuredTorque", 'azimuthMotor1Torque', tstart, tend)
az_torque_2 = await client.select_packed_time_series("lsst.sal.ATMCS.measuredTorque", 'azimuthMotor2Torque', tstart, tend)
el_torque = await client.select_packed_time_series("lsst.sal.ATMCS.measuredTorque", 'elevationMotorTorque', tstart, tend)
rot_torque = await client.select_packed_time_series("lsst.sal.ATMCS.measuredTorque", 'nasmyth2MotorTorque', tstart, tend)

```

```python
# Calculate the tracking errors
az_vals = np.array(az.values.tolist())[:,0]
el_vals = np.array(el.values.tolist())[:,0]
rot_vals = np.array(rot.values.tolist())[:,0]
times = np.array(az.values.tolist())[:,1]
times = times - times [0]
# The fits are much better if the time variable                                                                                                    
# is centered in the interval                                                                                                                      
fit_times = times - times[int(len(az.values[:, 1]) / 2)]

# Fit with a quartic
az_fit = np.polyfit(fit_times, az_vals, 4)
el_fit = np.polyfit(fit_times, el_vals, 4)
rot_fit = np.polyfit(fit_times, rot_vals, 2)

az_model = np.polyval(az_fit, fit_times)
el_model = np.polyval(el_fit, fit_times)
rot_model = np.polyval(rot_fit, fit_times)

# Errors in arcseconds
az_error = (az_vals - az_model) * 3600
el_error = (el_vals - el_model) * 3600
rot_error = (rot_vals - rot_model) * 3600

# Calculate RMS                                                                                                                      
az_rms = np.sqrt(np.mean(az_error * az_error))
el_rms = np.sqrt(np.mean(el_error * el_error))
rot_rms = np.sqrt(np.mean(rot_error * rot_error))

# Calculate Image impact RMS                                                                                                         
image_az_rms = az_rms * np.cos(el_vals[0] * np.pi / 180.0)
image_el_rms = el_rms
image_rot_rms = rot_rms * 280.0 * np.pi / 180.0 / 3600.0

```

```python
saveFilename = "/project/cslage/AuxTel/mount_graphs/Mount_Tracking_20220301_12Jul22.pdf"
figure = plt.figure(figsize=(11,16))
plt.subplots_adjust(wspace=0.5)
if saveFilename is not None:
    # Plotting                                                                                                                                     
    figure.clear()
    title = f"Mount Tracking 20220301, Azimuth = {az_vals[0]:.1f}, Elevation = {el_vals[0]:.1f}"
    plt.suptitle(title, fontsize=18)
    # Azimuth axis                                                                                                                                 
    plt.subplot(3, 3, 1)
    ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
    ax1.set_title("Azimuth axis", fontsize=16)
    ax1.axvline(az.index[0], color="red", linestyle="--")
    ax1.set_xticks([])
    ax1.set_ylabel("Degrees")
    plt.subplot(3, 3, 4)
    plt.plot(fit_times, az_error, color='red')

    plt.title(f"Azimuth RMS error = {az_rms:.2f} arcseconds\n"
              f"  Image RMS error = {image_az_rms:.2f} arcseconds", fontsize=10)
    plt.ylim(-10.0, 10.0)
    plt.xticks([])
    plt.ylabel("Arcseconds")
    plt.subplot(3, 3, 7)
    ax7 = az_torque_1['azimuthMotor1Torque'].plot(legend=True, color='blue')
    ax7 = az_torque_2['azimuthMotor2Torque'].plot(legend=True, color='green')
    ax7.axvline(az.index[0], color="red", linestyle="--")
    ax7.set_ylabel("Torque (motor current in amps)")

    # Elevation axis                                                                                                                               
    plt.subplot(3, 3, 2)
    ax2 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
    ax2.set_title("Elevation axis", fontsize=16)
    ax2.axvline(az.index[0], color="red", linestyle="--")
    ax2.set_xticks([])
    plt.subplot(3, 3, 5)
    plt.plot(fit_times, el_error, color='green')
    plt.title(f"Elevation RMS error = {el_rms:.2f} arcseconds\n"
              f"    Image RMS error = {image_el_rms:.2f} arcseconds", fontsize=10)
    plt.ylim(-10.0, 10.0)
    plt.xticks([])
    plt.subplot(3, 3, 8)
    ax8 = el_torque['elevationMotorTorque'].plot(legend=True, color='blue')
    ax8.axvline(az.index[0], color="red", linestyle="--")
    ax8.set_ylabel("Torque (motor current in amps)")

    # Nasmyth2 rotator axis                                                                                                                        
    plt.subplot(3, 3, 3)
    ax3 = rot['nasmyth2CalculatedAngle'].plot(legend=True, color='blue')
    ax3.set_title("Nasmyth2 axis", fontsize=16)
    ax3.axvline(az.index[0], color="red", linestyle="--")
    ax3.set_xticks([])
    plt.subplot(3, 3, 6)
    plt.plot(fit_times, rot_error, color='blue')
    plt.title(f"Nasmyth2 RMS error = {rot_rms:.2f} arcseconds\n"
              f"  Image RMS error <= {image_rot_rms:.2f} arcseconds",fontsize=10)
    plt.ylim(-10.0, 10.0)
    plt.xticks([])
    plt.subplot(3, 3, 9)
    ax9 = rot_torque['nasmyth2MotorTorque'].plot(legend=True, color='blue')
    ax9.axvline(az.index[0], color="red", linestyle="--")
    ax9.set_ylabel("Torque (motor current in amps)")
    plt.savefig(saveFilename)

```

```python tags=[]
# Plot it
fig = plt.figure(figsize = (8,8))
plt.subplots_adjust(wspace=0.5, hspace=0.5)
plt.suptitle(f"Mount Tracking Az Error 20220301\nStarting Az=0, Starting El=80", fontsize = 18)
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
#plt.savefig(f"./Az_Error_Debug_1_1Mar22_az0el80.pdf")
```

# El 70 

```python jupyter={"outputs_hidden": true} tags=[]
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T19:31:30", scale='utc')
tend = Time("2022-03-01T19:33:30", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220301\nStarting Az=0, Starting El=70", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az0el70.pdf")
```

# EL 60 

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T19:36:30", scale='utc')
tend = Time("2022-03-01T19:38:30", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220103\nStarting Az=0, Starting El=60", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az0el60.pdf")
```

# EL 50 

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T19:43:00", scale='utc')
tend = Time("2022-03-01T19:45:00", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220103\nStarting Az=0, Starting El=50", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az0el50.pdf")
```

# El 40 

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T19:48:00", scale='utc')
tend = Time("2022-03-01T19:50:00", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220103\nStarting Az=0, Starting El=40", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az0el40.pdf")
```

# El 25

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T19:53:00", scale='utc')
tend = Time("2022-03-01T19:55:00", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220301\nStarting Az=0, Starting El=25", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az0el25.pdf")
```

# EL 80 Back to 80

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T19:59:30", scale='utc')
tend = Time("2022-03-01T20:01:30", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220301\Back to Az=0, Starting El=80", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az0el80_Backto.pdf")
```

# Modifying azimuths Az = 10, EL = 80 

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T20:04:00", scale='utc')
tend = Time("2022-03-01T20:06:00", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220301\nStarting Az=10, Starting El=80", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az10el80.pdf")
```

# Az = 20 , EL = 80

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T20:09:00", scale='utc')
tend = Time("2022-03-01T20:11:00", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220224\nStarting Az=20, Starting El=80", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az20el80.pdf")
```

# Az = 30 , EL = 80

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T20:14:30", scale='utc')
tend = Time("2022-03-01T20:16:30", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220103\nStarting Az=30, Starting El=80", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az30el80.pdf")
```

```python

```

# Az = 45 , EL = 80

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T20:19:30", scale='utc')
tend = Time("2022-03-01T20:21:30", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220103\nStarting Az=45, Starting El=80", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az45el80.pdf")
```

```python

```

```python

```

# Az = 60 , EL = 80

```python
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T20:24:30", scale='utc')
tend = Time("2022-03-01T20:25:10", scale='utc')
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
plt.suptitle(f"Mount Tracking Az Error 20220103\nStarting Az=60, Starting El=80", fontsize = 18)
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
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az60el80.pdf")
```

```python

```
