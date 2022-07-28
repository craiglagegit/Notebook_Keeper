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

## AuxTel Mount fails - 18-Feb-21

In this notebook, investigate observed mount fails during the observing night 18-Feb-21

```python
import sys, time, os, asyncio

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from astropy.time import Time, TimeDelta
from astropy.coordinates import AltAz, ICRS, EarthLocation, Angle, FK5
import astropy.units as u
from lsst_efd_client import EfdClient
from lsst.pipe.tasks.quickFrameMeasurement import QuickFrameMeasurementTask
from lsst.geom import PointD
import lsst.daf.persistence as dafPersist

#from lsst.daf.persistence.dafPersist import Butler
```

```python
import lsst_efd_client
lsst_efd_client.__file__
```

```python
print(dir(EfdClient))
EfdClient.__module__
```

```python
# Set Cerro Pachon location
location = EarthLocation.from_geodetic(lon=-70.747698*u.deg,
                                       lat=-30.244728*u.deg,
                                       height=2663.0*u.m)
```

```python
# Get EFD client and bring in Lupton's unpacking code
client = EfdClient('summit_efd')

def merge_packed_time_series(packed_dataframe, base_field, stride=1, 
                             ref_timestamp_col="cRIO_timestamp", internal_time_scale="tai"):
    """Select fields that are time samples and unpack them into a dataframe.
            Parameters
            ----------
            packedDF : `pandas.DataFrame`
                packed data frame containing the desired data
            base_field :  `str`
                Base field name that will be expanded to query all
                vector entries.
            stride : `int`, optional
                Only use every stride value when unpacking.  Must be a factor
                of the number of packed values.
                (1 by default)
            ref_timestamp_col : `str`, optional
                Name of the field name to use to assign timestamps to unpacked
                vector fields (default is 'cRIO_timestamp').
            internal_time_scale : `str`, optional
                Time scale to use when converting times to internal formats
                ('tai' by default). Equivalent to EfdClient.internal_scale
        Returns
            -------
            result : `pandas.DataFrame`
                A `pandas.DataFrame` containing the results of the query.
            """
    
    packed_fields = [k for k in packed_dataframe.keys() if k.startswith(base_field)]
    packed_fields = sorted(packed_fields, key=lambda k: int(k[len(base_field):]))  # sort by pack ID
    npack = len(packed_fields)
    if npack%stride != 0:
        raise RuntimeError(f"Stride must be a factor of the number of packed fields: {stride} v. {npack}")
    packed_len = len(packed_dataframe)
    n_used = npack//stride   # number of raw fields being used
    output = np.empty(n_used*packed_len)
    times = np.empty_like(output, dtype=packed_dataframe[ref_timestamp_col][0])
    
    if packed_len == 1:
        dt = 0
    else:
        dt = (packed_dataframe[ref_timestamp_col][1] - packed_dataframe[ref_timestamp_col][0])/npack
    for i in range(0, npack, stride):
        i0 = i//stride
        output[i0::n_used] = packed_dataframe[f"{base_field}{i}"]
        times[i0::n_used] = packed_dataframe[ref_timestamp_col] + i*dt
     
    timestamps = Time(times, format='unix', scale=internal_time_scale).datetime64
    return pd.DataFrame({base_field:output, "times":times}, index=timestamps)
```

```python
async def MountTracking(expId, fail_limit = 0.5, makeGraph=True):
    # Find the time of exposure
    REPO_DIR = '/readonly/lsstdata/auxtel/base/auxtel/oods/butler/repo'
    butler = Butler(REPO_DIR)
    mData = butler.get('raw', dataId={'detector':0, 'expId':expId}).getMetadata()
    imgType = mData['IMGTYPE']
    expTime = mData['EXPTIME']
    if (imgType not in ['OBJECT', 'SKYEXP', 'ENGTEST']) or (expTime < 1.0):
        return True
    tStart = mData['DATE-BEG']
    tEnd = mData['DATE-END']
    # Get the data
    t_start = Time(tStart, scale='tai')
    t_end = Time(tEnd, scale='tai')
    mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                              t_start, t_end)
    nasmyth_position = await client.select_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", ['*'],
                                              t_start, t_end)
    time.sleep(5.0)
    az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1)
    el = merge_packed_time_series(mount_position, 'elevationCalculatedAngle', stride=1)
    rot = merge_packed_time_series(nasmyth_position, 'nasmyth2CalculatedAngle', stride=1)  
    # Calculate the tracking errors
    az_vals = np.array(az.values.tolist())[:,0]
    el_vals = np.array(el.values.tolist())[:,0]
    rot_vals = np.array(rot.values.tolist())[:,0]
    times = np.array(az.values.tolist())[:,1]
    times = times - times [0]

    # Fit with a quadratic
    az_fit = np.polyfit(times, az_vals, 2)
    el_fit = np.polyfit(times, el_vals, 2)
    rot_fit = np.polyfit(times, rot_vals, 2)

    az_model = az_fit[0] * times * times + az_fit[1] * times + az_fit[2]
    el_model = el_fit[0] * times * times + el_fit[1] * times + el_fit[2]
    rot_model = rot_fit[0] * times * times + rot_fit[1] * times + rot_fit[2]

    # Errors in arcseconds
    az_error = (az_vals - az_model) * 3600
    el_error = (el_vals - el_model) * 3600
    rot_error = (rot_vals - rot_model) * 3600

    # Calculate RMS
    az_rms = np.sqrt(np.mean(az_error * az_error))
    el_rms = np.sqrt(np.mean(el_error * el_error))
    rot_rms = np.sqrt(np.mean(rot_error * rot_error))
    # Plot it
    if makeGraph:
        fig = plt.figure(figsize = (16,8))
        plt.suptitle(f"Mount Tracking {expId}", fontsize = 18)
        plt.subplot(3,2,1)
        ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
        ax1.axvline(tStart, color="red", linestyle="--")
        ax1.set_xticks([])
        ax1.set_ylabel("Degrees")
        plt.subplot(3,2,3)
        ax3 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
        ax3.axvline(tStart, color="red", linestyle="--")
        ax3.set_xticks([])
        ax3.set_ylabel("Degrees")
        plt.subplot(3,2,5)
        ax5 = rot['nasmyth2CalculatedAngle'].plot(legend=True, color='blue')
        ax5.axvline(tStart, color="red", linestyle="--")
        ax5.set_ylabel("Degrees")

        plt.subplot(3,2,2)
        plt.plot(times, az_error, color='red')
        plt.title(f"Azimuth RMS error = {az_rms:.2f} arcseconds")
        plt.ylim(-10.0,10.0)
        plt.xticks([])
        plt.ylabel("ArcSeconds")
        plt.subplot(3,2,4)
        plt.plot(times, el_error, color='green')
        plt.title(f"Elevation RMS error = {el_rms:.2f} arcseconds")
        plt.ylim(-10.0,10.0)
        plt.xticks([])
        plt.ylabel("ArcSeconds")
        plt.subplot(3,2,6)
        plt.plot(times, rot_error, color='blue')
        plt.title(f"Nasmyth RMS error = {rot_rms:.2f} arcseconds")
        plt.ylim(-10.0,10.0)
        plt.ylabel("ArcSeconds")
        plt.savefig(f"/home/craiglagegit/DATA/mount_graphs/Mount_Errors_{expId}_20Feb21.pdf")
    if (az_rms > fail_limit) or (el_rms > fail_limit) or (rot_rms > fail_limit):
        return False
    else:
        return True

```

```python
await MountTracking(2021021800634)
```

```python
# Find the time of exposure
REPO_DIR = '/project/shared/auxTel/rerun/quickLook'
butler = dafPersist.Butler(REPO_DIR)
expId=2021021800621
exp = butler.get('quickLookExp', detector=0, expId=expId)
mData = exp.getMetadata()
tStart = mData['DATE-BEG']
tEnd = mData['DATE-END']


# Get the data
t_start = Time(tStart, scale='tai')
t_end = Time(tEnd, scale='tai')

mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                          t_start, t_end)
nasmyth_position = await client.select_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", ['*'],
                                          t_start, t_end)
torques = await client.select_time_series("lsst.sal.ATMCS.measuredTorque", ['*'],
                                          t_start, t_end)

az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1)
el = merge_packed_time_series(mount_position, 'elevationCalculatedAngle', stride=1)
rot = merge_packed_time_series(nasmyth_position, 'nasmyth2CalculatedAngle', stride=1) 
az_torque_1 =  merge_packed_time_series(torques, 'azimuthMotor1Torque', stride=1)
az_torque_2 =  merge_packed_time_series(torques, 'azimuthMotor2Torque', stride=1)
el_torque =  merge_packed_time_series(torques, 'elevationMotorTorque', stride=1)
rot_torque =  merge_packed_time_series(torques, 'nasmyth2MotorTorque', stride=1)

# Calculate the tracking errors
az_vals = np.array(az.values.tolist())[:,0]
el_vals = np.array(el.values.tolist())[:,0]
rot_vals = np.array(rot.values.tolist())[:,0]
times = np.array(az.values.tolist())[:,1]
times = times - times [0]

# Fit with a quadratic
az_fit = np.polyfit(times, az_vals, 2)
el_fit = np.polyfit(times, el_vals, 2)
rot_fit = np.polyfit(times, rot_vals, 2)

az_model = az_fit[0] * times * times + az_fit[1] * times + az_fit[2]
el_model = el_fit[0] * times * times + el_fit[1] * times + el_fit[2]
rot_model = rot_fit[0] * times * times + rot_fit[1] * times + rot_fit[2]

# Errors in arcseconds
az_error = (az_vals - az_model) * 3600
el_error = (el_vals - el_model) * 3600
rot_error = (rot_vals - rot_model) * 3600

# Calculate RMS
az_rms = np.sqrt(np.mean(az_error * az_error))
el_rms = np.sqrt(np.mean(el_error * el_error))
rot_rms = np.sqrt(np.mean(rot_error * rot_error))

# Plot it
fig = plt.figure(figsize = (16,16))
plt.suptitle(f"Mount Tracking {expId}", fontsize = 18)
# Azimuth axis
plt.subplot(3,3,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.set_title("Azimuth axis", fontsize=16)
ax1.axvline(tStart, color="red", linestyle="--")
ax1.set_xticks([])
ax1.set_ylabel("Degrees")
plt.subplot(3,3,4)
plt.plot(times, az_error, color='red')
plt.title(f"Azimuth RMS error = {az_rms:.2f} arcseconds")
plt.ylim(-10.0,10.0)
plt.xticks([])
plt.ylabel("ArcSeconds")
plt.subplot(3,3,7)
ax7 = az_torque_1['azimuthMotor1Torque'].plot(legend=True, color='blue')
ax7 = az_torque_2['azimuthMotor2Torque'].plot(legend=True, color='green')
ax7.axvline(tStart, color="red", linestyle="--")

# Elevation axis
plt.subplot(3,3,2)
ax2 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
ax2.set_title("Elevation axis", fontsize=16)
ax2.axvline(tStart, color="red", linestyle="--")
ax2.set_xticks([])
plt.subplot(3,3,5)
plt.plot(times, el_error, color='green')
plt.title(f"Elevation RMS error = {el_rms:.2f} arcseconds")
plt.ylim(-10.0,10.0)
plt.xticks([])
plt.subplot(3,3,8)
ax8 = el_torque['elevationMotorTorque'].plot(legend=True, color='blue')
ax8.axvline(tStart, color="red", linestyle="--")

# Nasmyth2 rotator axis
plt.subplot(3,3,3)
ax3 = rot['nasmyth2CalculatedAngle'].plot(legend=True, color='blue')
ax3.set_title("Nasmyth2 axis", fontsize=16)
ax3.axvline(tStart, color="red", linestyle="--")
ax3.set_xticks([])
plt.subplot(3,3,6)
plt.plot(times, rot_error, color='blue')
plt.title(f"Nasmyth RMS error = {rot_rms:.2f} arcseconds")
plt.ylim(-10.0,10.0)
plt.subplot(3,3,9)
ax9 = rot_torque['nasmyth2MotorTorque'].plot(legend=True, color='blue')
ax9.axvline(tStart, color="red", linestyle="--")

plt.savefig(f"/project/shared/auxTel/rerun/cslage/mount_plots/Mount_Errors_{expId}.pdf")

```

```python

expId=2021021800644
REPO_DIR = '/project/shared/auxTel/rerun/quickLook'
butler = dafPersist.Butler(REPO_DIR)
exp = butler.get('quickLookExp', detector=0, expId=expId)

arr = exp.getMaskedImage().getArrays()[0]
```

```python
qm_config = QuickFrameMeasurementTask.ConfigClass()
#qm_config.imageIsDispersed=True
qm = QuickFrameMeasurementTask(config=qm_config)
```

```python
result = qm.run(exp)
current_position = PointD(result.brightestObjCentroid[0], result.brightestObjCentroid[1])
print(result.brightestObjCentroid)
cpx = current_position[0]
cpy = current_position[1]
```

```python
result.imageMotion
```

```python
delta = 50
xp = 3783; yp = 1330
arr_plot = arr[yp-delta:yp+delta, xp-delta:xp+delta]
print(np.max(arr_plot), np.min(arr_plot))
```

```python
delta = 50
xp = 3183; yp = 1106
arr_plot = arr[yp-delta:yp+delta, xp-delta:xp+delta]
print(np.max(arr_plot), np.min(arr_plot))
```

```python
print(result.medianPsf)
print(result.brightestObj_xXyY)
```

```python
delta=100
plt.plot(arr[int(current_position[1]-delta):int(current_position[1]+delta), int(current_position[0])], color='red')
plt.plot(arr[int(current_position[1]), int(current_position[0]-delta):int(current_position[0]+delta)], color='green')
```

```python
from matplotlib.colors import LogNorm
# Now let's look at ithem
def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

fig, (ax1) = plt.subplots(ncols=1, figsize=(8,8))
img1 = ax1.imshow(arr_plot, norm=LogNorm(vmin=10, vmax=120000), cmap='gray')
#plt.plot([cpx-20, cpx+20], [cpy, cpy], color='red')
#plt.plot([cpx, cpx], [cpy-20, cpy+20], color='red')
colorbar(img1)
ax1.set_title("%d Object"%expId)
plt.tight_layout(h_pad=1)
plt.savefig(f"/project/shared/auxTel/rerun/cslage/mount_plots/Object_{expId}.pdf")


```

```python
from matplotlib.colors import LogNorm
# Now let's look at ithem
def colorbar(mappable):
    from mpl_toolkits.axes_grid1 import make_axes_locatable
    last_axes = plt.gca()
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    cbar = fig.colorbar(mappable, cax=cax)
    plt.sca(last_axes)
    return cbar

fig, (ax1) = plt.subplots(ncols=1, figsize=(8,8))
img1 = ax1.imshow(arr_plot, norm=LogNorm(vmin=10, vmax=120000), cmap='gray')
#plt.plot([cpx-20, cpx+20], [cpy, cpy], color='red')
#plt.plot([cpx, cpx], [cpy-20, cpy+20], color='red')
colorbar(img1)
ax1.set_title("%d brightestObjCentroid"%expId)
plt.tight_layout(h_pad=1)
plt.savefig(f"/home/craiglagegit/DATA/mount_graphs/BrightestObjCentroid_{expId}_25Feb21.pdf")

```

```python
plt.imshow(arr, cmap='gray')
```

```python
dx = 20
plt.imshow(np.transpose(arr[int(current_position[1]-dx):int(current_position[1]+dx), \
                            int(current_position[0]-dx):int(current_position[0]+dx)]), cmap='gray')
```

```python
expId=2021021800637
exp = butler.get('quickLookExp', detector=0, expId=expId)
psf = exp.getPsf()
psfShape = psf.computeShape()

ixx = psf.computeShape().getIxx()
iyy = psf.computeShape().getIyy()
ixx = np.sqrt(ixx)*2.355*.1
iyy = np.sqrt(iyy)*2.355*.1
print(f"Psf shape from imChar task (x,y) = ({ixx:.3f}, {iyy:.3f}) FWHM arcsec")


```

```python
TimeDelta?
```

```python
tStart = mData['MJD-BEG']
tEnd = mData['MJD-END']

print(ts, te, (te-ts)*86400)

TimeDelta(ts, te)
print((t_start - t_end)*86400)
```

```python
# Find the time of exposure
REPO_DIR = '/project/shared/auxTel/rerun/quickLook'
butler = dafPersist.Butler(REPO_DIR)
expId=2021021800201
exp = butler.get('quickLookExp', detector=0, expId=expId)
mData = exp.getMetadata()
for key in mData.keys():
    print(key, mData[key])
```

```python
# Find the time of exposure
REPO_DIR = '/project/shared/auxTel/rerun/quickLook'
butler = dafPersist.Butler(REPO_DIR)
expId=2021021800327
exp = butler.get('quickLookExp', detector=0, expId=expId)
mData = exp.getMetadata()
tStart = mData['DATE-BEG']
tEnd = mData['DATE-END']
imgType = mData['IMGTYPE']
expTime = (mData['MJD-END'] - mData['MJD-BEG']) * 86400
focusZ = mData['FOCUSZ']
if (imgType not in ['OBJECT', 'SKYEXP', 'ENGTEST']) or (expTime < 5.0) or (abs(focusZ) > 0.10):
    print("Not a good image")

qm_config = QuickFrameMeasurementTask.ConfigClass()
qm = QuickFrameMeasurementTask(config=qm_config)
result = qm.run(exp)
(xMotion, yMotion) = result.imageMotion

# Get the data
t_start = Time(tStart, scale='tai')
t_end = Time(tEnd, scale='tai')

mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                          t_start, t_end)
nasmyth_position = await client.select_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", ['*'],
                                          t_start, t_end)
torques = await client.select_time_series("lsst.sal.ATMCS.measuredTorque", ['*'],
                                          t_start, t_end)

az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1)
el = merge_packed_time_series(mount_position, 'elevationCalculatedAngle', stride=1)
rot = merge_packed_time_series(nasmyth_position, 'nasmyth2CalculatedAngle', stride=1) 
az_torque_1 =  merge_packed_time_series(torques, 'azimuthMotor1Torque', stride=1)
az_torque_2 =  merge_packed_time_series(torques, 'azimuthMotor2Torque', stride=1)
el_torque =  merge_packed_time_series(torques, 'elevationMotorTorque', stride=1)
rot_torque =  merge_packed_time_series(torques, 'nasmyth2MotorTorque', stride=1)

# Calculate the tracking errors
az_vals = np.array(az.values.tolist())[:,0]
el_vals = np.array(el.values.tolist())[:,0]
rot_vals = np.array(rot.values.tolist())[:,0]
times = np.array(az.values.tolist())[:,1]
times = times - times [0]

# Fit with a quadratic
az_fit = np.polyfit(times, az_vals, 2)
el_fit = np.polyfit(times, el_vals, 2)
rot_fit = np.polyfit(times, rot_vals, 2)

az_model = az_fit[0] * times * times + az_fit[1] * times + az_fit[2]
el_model = el_fit[0] * times * times + el_fit[1] * times + el_fit[2]
rot_model = rot_fit[0] * times * times + rot_fit[1] * times + rot_fit[2]

# Errors in arcseconds
az_error = (az_vals - az_model) * 3600
el_error = (el_vals - el_model) * 3600
rot_error = (rot_vals - rot_model) * 3600

# Calculate RMS
az_rms = np.sqrt(np.mean(az_error * az_error))
el_rms = np.sqrt(np.mean(el_error * el_error))
rot_rms = np.sqrt(np.mean(rot_error * rot_error))

# Plot it
fig = plt.figure(figsize = (16,16))
plt.suptitle(f"Mount Tracking {expId}", fontsize = 18)
# Azimuth axis
plt.subplot(3,3,1)
ax1 = az['azimuthCalculatedAngle'].plot(legend=True, color='red')
ax1.set_title("Azimuth axis", fontsize=16)
ax1.axvline(tStart, color="red", linestyle="--")
ax1.set_xticks([])
ax1.set_ylabel("Degrees")
plt.subplot(3,3,4)
plt.plot(times, az_error, color='red')
plt.title(f"Azimuth RMS error = {az_rms:.1f} arcseconds \n \
            Estimated image motion= {yMotion:.1f} arcseconds")
plt.ylim(-10.0,10.0)
plt.xticks([])
plt.ylabel("ArcSeconds")
plt.subplot(3,3,7)
ax7 = az_torque_1['azimuthMotor1Torque'].plot(legend=True, color='blue')
ax7 = az_torque_2['azimuthMotor2Torque'].plot(legend=True, color='green')
ax7.axvline(tStart, color="red", linestyle="--")

# Elevation axis
plt.subplot(3,3,2)
ax2 = el['elevationCalculatedAngle'].plot(legend=True, color='green')
ax2.set_title("Elevation axis", fontsize=16)
ax2.axvline(tStart, color="red", linestyle="--")
ax2.set_xticks([])
plt.subplot(3,3,5)
plt.plot(times, el_error, color='green')
plt.title(f"Elevation RMS error = {el_rms:.1f} arcseconds \n \
            Estimated image motion= {xMotion:.1f} arcseconds")
plt.ylim(-10.0,10.0)
plt.xticks([])
plt.subplot(3,3,8)
ax8 = el_torque['elevationMotorTorque'].plot(legend=True, color='blue')
ax8.axvline(tStart, color="red", linestyle="--")

# Nasmyth2 rotator axis
plt.subplot(3,3,3)
ax3 = rot['nasmyth2CalculatedAngle'].plot(legend=True, color='blue')
ax3.set_title("Nasmyth2 axis", fontsize=16)
ax3.axvline(tStart, color="red", linestyle="--")
ax3.set_xticks([])
plt.subplot(3,3,6)
plt.plot(times, rot_error, color='blue')
plt.title(f"Nasmyth RMS error = {rot_rms:.2f} arcseconds")
plt.ylim(-10.0,10.0)
plt.subplot(3,3,9)
ax9 = rot_torque['nasmyth2MotorTorque'].plot(legend=True, color='blue')
ax9.axvline(tStart, color="red", linestyle="--")

plt.savefig(f"/project/shared/auxTel/rerun/cslage/mount_plots/Mount_Errors_{expId}.pdf")

```

```python

```
