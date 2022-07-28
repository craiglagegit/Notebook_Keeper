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

## Elevation Motor Torques
In this notebook, I look at the impact of removing the AuxTel main\
mirror cover on the torques needed to move the elevation axis.

```python
import sys, time, os

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.fft import fft, fftfreq

from astropy.time import Time, TimeDelta
from lsst_efd_client import EfdClient, merge_packed_time_series
```

We'll access the EFD instance deployed at NCSA.

```python
#client = EfdClient('summit_efd')
client = EfdClient('ldf_stable_efd')
```

```python
t_ends = [Time("2021-06-10T05:00:00", scale='tai') , Time("2021-02-18T10:00:00", scale='tai')]
notes = ["Mirror cover removed", "Mirror cover in place"]
```

```python
# Get the data first, which takes a while
nsec = 6.0*3600 # how many seconds of data to retrieve
torqueLists = []
angleLists = []
for i, t_end in enumerate(t_ends):
    t_start = t_end - TimeDelta(nsec, format='sec')
    elevation_torque = await client.select_time_series("lsst.sal.ATMCS.measuredTorque", ['*'],
                                                  t_start, t_end)
    torque = merge_packed_time_series(elevation_torque, 'elevationMotorTorque', stride=1)
    elevation_angle = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                                  t_start, t_end)
    angle = merge_packed_time_series(elevation_angle, 'elevationCalculatedAngle', stride=1)
    torqueList = torque.values.tolist()
    angleList = angle.values.tolist()
    torqueLists.append(torqueList)
    angleLists.append(angleList)
```

```python
# Now plot it
plt.figure(figsize = (16,8))
for i, t_end in enumerate(t_ends):
    date = t_end.isot.split('T')[0]
    torqueList = torqueLists[i]
    angleList = angleLists[i]
    plt.subplot(1,2,i + 1)
    plt.title(f"Elevation angle vs Torque - {date}\n {notes[i]}", fontsize = 18)
    plt.plot(np.array(angleList)[:,0],np.array(torqueList)[:,0])
    plt.arrow(20, 2.5, 20,0, width=0.1,head_length = 5.0, color='green')
    plt.arrow(85, -2.5, -20,0, width=0.1,head_length = 5.0, color='green')
    plt.xlabel("Elevation angle(degrees)", fontsize = 18)
    plt.ylabel("Torque (amps)", fontsize = 18)
plt.savefig("/project/cslage/AuxTel/torques/Elevation_Torque_vs_Angle_18Jun21.pdf")
```

```python

```
