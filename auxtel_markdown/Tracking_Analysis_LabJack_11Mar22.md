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

## Anayzing mount tracking with accelerometer data

Craig Lage - Mar 11, 2022

```python
import sys, time, datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pickle as pkl
from astropy.time import Time, TimeDelta
import astropy.units as u
from lsst_efd_client import EfdClient
from lsst.daf.butler import Butler
```

```python tags=[]
# Get EFD client and butler
client = EfdClient('summit_efd')
butler = Butler('/repo/LATISS', collections="LATISS/raw/all")
```

```python
# Choose an expId
expId = 2022030800027

makeGraph = True
doOffset = True #Enter a time offset tocompensate for astropy TAI issue 
```

```python tags=[]
# Find the time of exposure   

start = time.time()
mData = butler.get('raw.metadata', detector=0, exposure=expId)
imgType = mData['IMGTYPE']
tStart = Time(mData['DATE-BEG'], format='isot', scale='tai')
tEnd = Time(mData['DATE-END'], format='isot', scale='tai')
elevation = mData['ELSTART']
azimuth = mData['AZSTART']
print(f"expId = {expId}, imgType = {imgType}, Times = {tStart}, {tEnd}")
if (imgType not in ['OBJECT', 'SKYEXP', 'ENGTEST', 'DARK']):
    sys.exit()
end = time.time()
elapsed = end-start
print(f"Elapsed time for butler query = {elapsed}")
start = time.time()
print(f"tStart={tStart}, tEnd={tEnd}")

# Get the data                                                                                             
t_start = tStart.utc
t_end = tEnd.utc
print(f"t_start={t_start}, t_end={t_end}")

az = await client.select_packed_time_series('lsst.sal.ATMCS.mount_AzEl_Encoders', \
                                            'azimuthCalculatedAngle', t_start, t_end)
el = await client.select_packed_time_series('lsst.sal.ATMCS.mount_AzEl_Encoders', \
                                            'elevationCalculatedAngle', t_start, t_end)
if doOffset:
    offset = (t_start.jd - az.index[0].to_julian_date()) * 86400.0
    initial_offset = offset
    az.index += pd.DateOffset(seconds=offset)
    el.index += pd.DateOffset(seconds=offset)
else:
    initial_offset = 0.0
    
rot = await client.select_packed_time_series('lsst.sal.ATMCS.mount_Nasmyth_Encoders', \
                                            'nasmyth2CalculatedAngle', t_start, t_end)
if doOffset:
    offset = (t_start.jd - rot.index[0].to_julian_date()) * 86400.0
    rot.index += pd.DateOffset(seconds=offset)

az_torque_1 = await client.select_packed_time_series('lsst.sal.ATMCS.measuredTorque', \
                                            'azimuthMotor1Torque', t_start, t_end)
az_torque_2 = await client.select_packed_time_series('lsst.sal.ATMCS.measuredTorque', \
                                            'azimuthMotor2Torque', t_start, t_end)
el_torque = await client.select_packed_time_series('lsst.sal.ATMCS.measuredTorque', \
                                            'elevationMotorTorque', t_start, t_end)
rot_torque = await client.select_packed_time_series('lsst.sal.ATMCS.measuredTorque', \
                                            'nasmyth2MotorTorque', t_start, t_end)
if doOffset:
    offset = (t_start.jd - az_torque_1.index[0].to_julian_date()) * 86400.0
    az_torque_1.index += pd.DateOffset(seconds=offset)
    az_torque_2.index += pd.DateOffset(seconds=offset)
    el_torque.index += pd.DateOffset(seconds=offset)
    rot_torque.index += pd.DateOffset(seconds=offset)

end = time.time()
elapsed = end-start
print(f"Elapsed time to get the data = {elapsed}")
start = time.time()

# Calculate the tracking errors                                                                            
az_vals = np.array(az.values.tolist())[:,0]
el_vals = np.array(el.values.tolist())[:,0]
rot_vals = np.array(rot.values.tolist())[:,0]
times = np.array(az.values.tolist())[:,1]
rot_times = np.array(rot.values.tolist())[:,1]
times = times - times[0]
rot_times = rot_times - rot_times[0]
print("LengthAz", len(az_vals))
# Fit with a quartic
az_fit = np.polyfit(times, az_vals, 4)
el_fit = np.polyfit(times, el_vals, 4)
rot_fit = np.polyfit(rot_times, rot_vals, 2)
az_model = az_fit[0] * times * times * times * times + az_fit[1] * times * times * times \
+ az_fit[2] * times *times + az_fit[3] * times + az_fit[4]
el_model = el_fit[0] * times * times * times * times + el_fit[1] * times * times * times \
+ el_fit[2] * times * times + el_fit[3] * times + el_fit[4]
rot_model = rot_fit[0] * rot_times * rot_times + rot_fit[1] * rot_times + rot_fit[2]

# Errors in arcseconds                                                                                     
az_error = (az_vals - az_model) * 3600
el_error = (el_vals - el_model) * 3600
rot_error = (rot_vals - rot_model) * 3600

# Calculate RMS                                                                                            
az_rms = np.sqrt(np.mean(az_error * az_error))
el_rms = np.sqrt(np.mean(el_error * el_error))
rot_rms = np.sqrt(np.mean(rot_error * rot_error))

end = time.time()
elapsed = end-start
print(f"Elapsed time for error calculations = {elapsed}")
start = time.time()

# Unpickle the accel dataframe
file = open(f'/scratch/labJackData/{expId}.pkl', 'rb')
df = pkl.load(file)
file.close()

# Add a buffer to the accel times
t_accel_start = t_start - TimeDelta(1.0, format='sec')
t_accel_end = t_end + TimeDelta(1.0, format='sec')
if makeGraph:

    fig = plt.figure(figsize = (16,16))
    plt.subplots_adjust(hspace=0.5)
    offset = (t_start.jd - az.index[0].to_julian_date()) * 86400.0
    plt.suptitle(f"Mount Tracking {expId}, Azimuth = {azimuth:.1f}, Elevation = {elevation:.1f}\nTime offset = {initial_offset:.2f} seconds\n Fan = 0%", fontsize = 18)
    # Azimuth axis                                                                                         
    plt.subplot(4,3,1)
    ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
    ax1.set_title("Azimuth axis", fontsize=16)
    ax1.axvline(t_start.to_datetime(), color="green", linestyle="--")
    ax1.axvline(t_end.to_datetime(), color="red", linestyle="--")
    #ax1.set_xticks([])
    ax1.set_ylabel("Degrees")
    plt.subplot(4,3,4)
    plt.plot(times, az_error, color='red')
    plt.title(f"Azimuth RMS error = {az_rms:.2f} arcseconds")
    plt.ylim(-10.0,10.0)
    #plt.xticks([])
    plt.ylabel("ArcSeconds")
    plt.subplot(4,3,7)
    ax7 = az_torque_1['azimuthMotor1Torque'].plot(color='blue')
    ax7 = az_torque_2['azimuthMotor2Torque'].plot(color='green')
    ax7.set_title("Azimuth Motor Torque", fontsize=12)
    ax7.axvline(t_start.to_datetime(), color="green", linestyle="--")
    ax7.axvline(t_end.to_datetime(), color="red", linestyle="--")
    plt.subplot(4,3,10)
    ax10 = df['AZM2'].plot(color='red')
    ax10.set_title("M2 Az accel Std = %.6f g"%np.std(df['AZM2']), fontsize=12)
    ax10.axvline(t_start.to_datetime(), color="green", linestyle="--")
    ax10.axvline(t_end.to_datetime(), color="red", linestyle="--")
    ax10.set_xlim(t_accel_start.to_datetime(), t_accel_end.to_datetime())

    # Elevation axis                                                                                       
    plt.subplot(4,3,2)
    ax2 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
    ax2.set_title("Elevation axis", fontsize=16)
    ax2.axvline(t_start.to_datetime(), color="green", linestyle="--")
    ax2.axvline(t_end.to_datetime(), color="red", linestyle="--")
    #ax2.set_xticks([])
    plt.subplot(4,3,5)
    plt.plot(times, el_error, color='green')
    plt.title(f"Elevation RMS error = {el_rms:.2f} arcseconds")
    plt.ylim(-10.0,10.0)
    #plt.xticks([])
    plt.subplot(4,3,8)
    ax8 = el_torque['elevationMotorTorque'].plot(color='blue')
    ax8.set_title("Elevation Motor Torque", fontsize=12)
    ax8.axvline(t_start.to_datetime(), color="green", linestyle="--")
    ax8.axvline(t_end.to_datetime(), color="red", linestyle="--")
    plt.subplot(4,3,11)
    ax11 = df['ELM2'].plot(color='red')
    ax11.set_title("M2 El accel Std = %.6f g"%np.std(df['ELM2']), fontsize=12)
    ax11.axvline(t_start.to_datetime(), color="green", linestyle="--")
    ax11.axvline(t_end.to_datetime(), color="red", linestyle="--")
    ax11.set_xlim(t_accel_start.to_datetime(), t_accel_end.to_datetime())

    # Nasmyth2 rotator axis                                                                                
    plt.subplot(4,3,3)
    ax3 = rot['nasmyth2CalculatedAngle'].plot(legend=True, color='blue')
    ax3.set_title("Nasmyth2 axis", fontsize=16)
    ax3.axvline(t_start.to_datetime(), color="green", linestyle="--")
    ax3.axvline(t_end.to_datetime(), color="red", linestyle="--")
    #ax3.set_xticks([])
    plt.subplot(4,3,6)
    plt.plot(rot_times, rot_error, color='blue')
    plt.title(f"Nasmyth RMS error = {rot_rms:.2f} arcseconds")
    plt.ylim(-10000.0,10000.0)
    plt.subplot(4,3,9)
    ax9 = rot_torque['nasmyth2MotorTorque'].plot(legend=True, color='blue')
    ax9.set_title("Nasmyth2 Motor Torque", fontsize=12)
    ax9.axvline(t_start.to_datetime(), color="green", linestyle="--")
    ax9.axvline(t_end.to_datetime(), color="red", linestyle="--")
    plt.subplot(4,3,12)
    ax12 = df['ZM2'].plot(color='red')
    ax12.set_title("M2 Optical Axis accel Std = %.6f g"%np.std(df['ZM2']), fontsize=12)
    ax12.axvline(t_start.to_datetime(), color="green", linestyle="--")
    ax12.axvline(t_end.to_datetime(), color="red", linestyle="--")
    ax12.set_xlim(t_accel_start.to_datetime(), t_accel_end.to_datetime())


    plt.savefig(f"/scratch/labJackData/Mount_Accel_{expId}_11Mar22.pdf")                 

end = time.time()
elapsed = end-start
print(f"Elapsed time for plots = {elapsed}")
start = time.time()

        
```

```python

```
