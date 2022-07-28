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

## Hexapod Position Offset 11-Jun-21

Investigate how much changing the XY hexapod offset impacts the image position\
Images from the night of 20210608

```python
import sys, time, os, asyncio

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from lsst.pipe.tasks.quickFrameMeasurement import QuickFrameMeasurementTask
from lsst.geom import PointD
import lsst.daf.persistence as dafPersist
```

```python
expId=2021060800432
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
