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
%matplotlib inline
import pandas as pd
from astropy.time import Time, TimeDelta
from lsst_efd_client import EfdClient

import numpy, scipy.optimize
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

```python tags=[]
# Get EFD client and bring in Lupton's unpacking code
client = EfdClient('summit_efd')
```

```python
# Define functions 
```

```python
def get_period(az_error):
    fs = 100 
    acf = np.correlate(az_error, az_error, 'full')[-len(az_error):]
    inflection = np.diff(np.sign(np.diff(acf))) # 
    peaks = (inflection < 0).nonzero()[0] + 1 # Find where they are negative
    delay = peaks[acf[peaks].argmax()] # Of those, find the index with the maximum value

    signal_freq = fs/delay 
    print(f'Frequency is {signal_freq} Hz')
    period = 1/signal_freq
    print(f'Period is {period} s')
```

```python
def get_frequency_fft(az_error, times):
    ft = np.fft.rfft(az_error)
    freqs = np.fft.rfftfreq(len(az_error), times[1]-times[0]) # Get frequency axis from the time axis
    mags = abs(ft) 
    plt.plot(freqs, mags)
    plt.xlim(0,0.25)
    plt.show()
    print(f'Frequency is {freqs[mags.argmax()]} \n'
     f'Period is {1/freqs[mags.argmax()]} s')
```

```python
def fit_sin(times, az_error):
    '''Fit sin to the input az_error time sequence, and return fitting parameters "amp", "omega", "phase", "offset", "freq", "period" and "fitfunc"'''
    tt = (times)
    yy = (az_error)
    ff = numpy.fft.fftfreq(len(tt), (tt[1]-tt[0]))  
    Fyy = abs(numpy.fft.fft(yy))
    # Guess parameters for the curve fit. 
    guess_freq = abs(ff[numpy.argmax(Fyy[1:])+1])   # excluding the zero frequency "peak",
    guess_amp = numpy.std(yy) * 2.**0.5
    guess_offset = numpy.mean(yy)
    guess = numpy.array([guess_amp, 2.*numpy.pi*guess_freq, 0., guess_offset])

    def sinfunc(t, A, w, p, c):  return A * numpy.sin(w*t + p) + c
    popt, pcov = scipy.optimize.curve_fit(sinfunc, tt, yy, p0=guess)
    A, w, p, c = popt
    f = w/(2.*numpy.pi)
    fitfunc = lambda t: A * numpy.sin(w*t + p) + c
    return {"amp": A, "omega": w, "phase": p, "offset": c, "freq": f, 
            "period": 1./f, "fitfunc": fitfunc, "maxcov": numpy.max(pcov), 
            "rawres": (guess,popt,pcov)}

```

```python
periods = []
speeds = []
```

# EL 80

```python tags=[]
# Blow up the first two minutes to look at the periodic error
tstart = Time("2022-03-01T19:27:00", scale='utc')
tend = Time("2022-03-01T19:29:00", scale='utc')
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

P80_speed = abs((az_vals[0] - az_vals[-1]) / (times[0] - times[-1]))
print(P80_speed)
```

```python
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
plt.ylim(-2.0,2.0)
plt.xticks([])
plt.ylabel("ArcSeconds")
plt.subplot(2,2,4)
plt.plot(times, el_error, color='green')
plt.title(f"Elevation RMS error = {el_rms:.2f} arcseconds")
plt.ylim(-10.0,10.0)
plt.xticks([])
plt.ylabel("ArcSeconds")
plt.savefig(f"./Az_Error_Debug_1_1Mar22_az0el80.pdf")
```

```python
results_80El = fit_sin(times,az_error)
P80EL = results_80El['period']
A80EL = results_80El['amp']
periods.append(P80EL)
speeds.append(P80_speed)
print(P80EL, A80EL, P80_speed)
```

# El 70 

```python tags=[]
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

P70_speed = abs((az_vals[0] - az_vals[-1]) / (times[0] - times[-1]))
print(P70_speed)
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

```python
results_70El = fit_sin(times,az_error)
P70EL = results_70El['period']
A70EL = results_70El['amp']
periods.append(P70EL)
speeds.append(P70_speed)
print(P70EL, A70EL, P70_speed)
```

# EL 60 

```python tags=[]
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

P60_speed = abs((az_vals[0] - az_vals[-1]) / (times[0] - times[-1]))
print(P60_speed)
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

```python
results_60El = fit_sin(times,az_error)
P60EL = results_60El['period']
A60EL = results_60El['amp']
periods.append(P60EL)
speeds.append(P60_speed)
print(P60EL, A60EL, P60_speed)

```

# EL 50 

```python tags=[]
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

P50_speed = abs((az_vals[0] - az_vals[-1]) / (times[0] - times[-1]))
print(P50_speed)
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

```python
results_50El = fit_sin(times,az_error)
P50EL = results_50El['period']
A50EL = results_50El['amp']
periods.append(P50EL)
speeds.append(P50_speed)
print(P50EL, A50EL, P50_speed)
```

# El 40 

```python tags=[]
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

P40_speed = abs((az_vals[0] - az_vals[-1]) / (times[0] - times[-1]))
print(P40_speed)
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

```python
results_40El = fit_sin(times,az_error)
P40EL = results_40El['period']
A40EL = results_40El['amp']
periods.append(P40EL)
speeds.append(P40_speed)
print(P40EL, A40EL, P40_speed)
```

# El 25

```python tags=[]
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

P25_speed = abs((az_vals[0] - az_vals[-1]) / (times[0] - times[-1]))
print(P25_speed)
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

```python
results_25El = fit_sin(times,az_error)
P25EL = results_25El['period']
A25EL = results_25El['amp']
periods.append(P25EL)
speeds.append(P25_speed)
print(P25EL, A25EL, P25_speed)
```

```python
periods # Seconds
```

```python
speeds # deg/sec
```

```python
period_in_arcseconds = np.array(periods) * np.array(speeds) * 3600
```

```python
angles = [80, 70, 60, 50, 40, 25]
for i, angle in enumerate(angles):
    print(f"El = {angle}, period = {period_in_arcseconds[i]:.1f} arseconds")

```

# EL 80 - Back to 80

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
plt.ylim(-2.0,2.0)
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

```python
results_80El2 = fit_sin(times,az_error)
P80EL2 = results_80El2['period']
A80EL2 = results_80El2['amp']
print(P80EL2, A80EL2)
```

# Plot amplitude and period vs elevation (AZ = 0)

```python
Periods = [P80EL, P70EL, P60EL, P50EL, P40EL, P25EL, P80EL2]
print(Periods)
Amplituds = [A80EL, A70EL, abs(A60EL), A50EL, A40EL, A25EL, abs(A80EL2)]
print(Amplituds)
elevation = [80, 70, 60, 50, 40, 25, 80]
```

```python
fig, axs = plt.subplots(2,1, figsize = (10,10))
plt.suptitle("Oscillations Period and Amplitude \n vs Elevation @ Az = 0", fontsize = 18)
axs[0].plot(elevation, Periods, "xk", color='green', markersize=10)
axs[0].set_ylabel("Period (s)")
#axs[0].set_xlabel("Elevation (degrees)")

axs[1].plot(elevation, Amplituds, "*r", color = 'blue', markersize=10)
axs[1].set_ylabel("Amplitud (arcsec)")
axs[1].set_xlabel("Elevation (degrees)")
plt.savefig(f"/project/cslage/AuxTel/mount_graphs/Oscillations_PeriodandAmplitude_vsElevation_01Mar22.pdf")
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

```python
results_az10 = fit_sin(times,az_error)
P10AZ = results_az10['period']
A10AZ = results_az10['amp']
print(P10AZ, A10AZ)
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

```python tags=[]
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

```python
results_az20 = fit_sin(times,az_error)
P20AZ = results_az20['period']
A20AZ = results_az20['amp']
print(P20AZ, A20AZ)
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
results_az30 = fit_sin(times,az_error)
P30AZ = results_az30['period']
A30AZ = results_az30['amp']
print(P30AZ, A30AZ)
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
results_az45 = fit_sin(times,az_error)
P45AZ = results_az45['period']
A45AZ = results_az45['amp']
print(P45AZ, A45AZ)
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
results_az60 = fit_sin(times,az_error)
P60AZ = results_az60['period']
A60AZ = results_az60['amp']
print(P60AZ, A60AZ)
```

# Plot amplitude and period vs azimuth (EL = 80) 

```python
Periods2 = [P10AZ, P20AZ, P30AZ, P45AZ, P80EL2]
print(Periods2)
Amplituds2 = [A10AZ, abs(A20AZ), A30AZ, A45AZ, abs(A80EL2)]
print(Amplituds2)
azimuth = [10, 20, 30, 45, 0]
```

```python
fig, axs = plt.subplots(2,1, figsize = (10,10))
plt.suptitle("Oscillations Period and Amplitude \n vs azimuth @ EL = 80", fontsize = 18)
axs[0].plot(azimuth, Periods2, "xk", color='green', markersize=10)
axs[0].set_ylabel("Period (s)")
#axs[0].set_xlabel("Elevation (degrees)")

axs[1].plot(azimuth, Amplituds2, "*r", color = 'blue', markersize=10)
axs[1].set_ylabel("Amplitud (arcsec)")
axs[1].set_xlabel("Azimuth (degrees)")
plt.savefig(f"./1Mar22Oscillations_PeriodandAmplitude_vsAzimuth.pdf")
```

```python

```
