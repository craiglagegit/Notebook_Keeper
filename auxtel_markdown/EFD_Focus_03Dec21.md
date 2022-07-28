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

# AuxTel Focus Study - 03-Dec-21

In this notebook, investigate focus settings and temp on 03-Dec-21

```python
import sys, time, os, asyncio

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
%matplotlib widget
import pandas as pd
from astropy.time import Time, TimeDelta
from lsst.daf.butler import Butler
```

```python
# Get EFD client
from lsst_efd_client import EfdClient
from lsst_efd_client import  __version__ as efdVersion
print(efdVersion)
client = EfdClient('ldf_stable_efd')
```

```python
# Get the butler
butler = Butler('/repo/main', collections="LATISS/raw/all")
```

```python
before = 15.0
tai_offset = 37.0

expId = 2021090800134
mData = butler.get('raw.metadata', detector=0, exposure=expId)
print(f"{expId} \t {mData['TIMESYS']} \t {mData['DATE']} \t {mData['DATE-BEG']} \t {mData['DATE-END']}")
date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
elevation = mData['ELSTART']
print(elevation)
start = date_beg.utc - TimeDelta(before, format='sec') + TimeDelta(tai_offset, format='sec')
end = date_beg.utc + TimeDelta(tai_offset, format='sec')
```

```python
total_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "total", start, end)
print(total_off.values[-1][0])
```

```python
z_position = await client.select_time_series("lsst.sal.ATHexapod.command_moveToPosition", "z", start, end)
print(z_position.values[-1][0])
```

```python
temp_air = await client.select_time_series("lsst.sal.ESS.temperature4Ch", "temperatureC02", start, end)
temp_truss = await client.select_time_series("lsst.sal.ESS.temperature4Ch", "temperatureC03", start, end)
temp_m2 = await client.select_time_series("lsst.sal.ESS.temperature4Ch", "temperatureC04", start, end)
temp_before = 60.0
temp_start = date_beg.utc - TimeDelta(temp_before, format='sec') + TimeDelta(tai_offset, format='sec')
temp_ext = await client.select_time_series("lsst.sal.WeatherStation.airTemperature", "avg1M", temp_start, end)
print(temp_air.values[-1][0])
print(temp_truss.values[-1][0])
print(temp_m2.values[-1][0])
print(temp_ext.values[-1][0])
```

```python
before = 15.0
temp_before = 60.0
tai_offset = 0.0#37.0

els = []
offs = []
poss = []

filename = '/project/cslage/AuxTel/efd_temp/EFD_Temp_20210908.txt'
#outfile = open(filename, 'w')
#outfile.write(f"expId\t\tElevation\tOffset\t\tHex_z\t\tT_air\tT_truss\tT_M2\n")

dayObs = 20210908
seqNos = [128, 134, 145, 149, 153, 161, 165, 489, 614, 641, 793]
dayObs = 20210909
seqNos = [152, 243, 348, 470, 542, 674, 773, 800]
dayObs = 20211005
seqNos = [ 297, 302, 307, 310, 316, 398, 415, 422, 662]
dayObs = 20211006
seqNos = [146, 545, 552]
dayObs = 20211102
seqNos = [93, 333, 346, 351, 374, 377, 383, 399, 498, 567]
dayObs = 20211103
seqNos = [70, 161, 176, 288, 446, 551, 621]

for seqNo in seqNos:
    expId = dayObs * 100000 + seqNo
    mData = butler.get('raw.metadata', detector=0, exposure=expId)
    date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
    elevation = mData['ELSTART']
    els.append(elevation)
    start = date_beg.utc - TimeDelta(before, format='sec') + TimeDelta(tai_offset, format='sec')
    temp_start = date_beg.utc - TimeDelta(temp_before, format='sec') + TimeDelta(tai_offset, format='sec')
    end = date_beg.utc + TimeDelta(tai_offset, format='sec')
    
    disp_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "disperser", start, end)
    disp_off = disp_off.values[-1][0]
    filter_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "filter", start, end)
    filter_off = filter_off.values[-1][0]
    user_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "userApplied", start, end)
    user_off = user_off.values[-1][0]
    total_off = await client.select_time_series("lsst.sal.ATAOS.logevent_focusOffsetSummary", "total", start, end)
    total_off = total_off.values[-1][0]
    offs.append(total_off)
    z_position = await client.select_time_series("lsst.sal.ATHexapod.command_moveToPosition", "z", start, end)
    z_position = z_position.values[-1][0]
    poss.append(z_position)
    temp_air = await client.select_time_series("lsst.sal.ESS.temperature4Ch", "temperatureC02", start, end)
    temp_truss = await client.select_time_series("lsst.sal.ESS.temperature4Ch", "temperatureC03", start, end)
    temp_m2 = await client.select_time_series("lsst.sal.ESS.temperature4Ch", "temperatureC04", start, end)
    temp_air = temp_air.values[-1][0]
    temp_truss = temp_truss.values[-1][0]
    temp_m2 = temp_m2.values[-1][0]
    temp_ext = await client.select_time_series("lsst.sal.WeatherStation.airTemperature", "avg1M", temp_start, end)
    temp_ext = temp_ext.values[-1][0]
    wind_spd = await client.select_time_series("lsst.sal.WeatherStation.windSpeed", "avg10M", temp_start, end)
    wind_spd = wind_spd.values[-1][0]
    wind_dir = await client.select_time_series("lsst.sal.WeatherStation.windDirection", "avg10M", temp_start, end)
    wind_dir = wind_dir.values[-1][0]
    try:
        dimm_fwhm = await client.select_time_series("lsst.sal.DIMM.logevent_dimmMeasurement", "fwhm", temp_start, end)
        dimm_fwhm = dimm_fwhm.values[-1][0]
    except:
        dimm_fwhm = None
    #outfile.write(f"{expId}\t{elevation:.4f}\t\t{total_off:.6f}\t{z_position:.6f}\t{temp_air:.2f}\t{temp_truss:.2f}\t{temp_m2:.2f}\n")

    print(expId, elevation, total_off, z_position, temp_air, temp_truss, temp_m2)
    print(disp_off, filter_off, user_off, temp_ext, wind_spd, wind_dir, dimm_fwhm)
#outfile.close()
```

```python
fig = plt.figure(figsize = (8,6))
plt.suptitle(f"Focus vs Elevation - dayObs {dayObs}", fontsize = 18)
plt.subplots_adjust(wspace=0.5)
plt.subplot(1,2,1)
plt.scatter(els, offs)
plt.xlabel("Elevation(Degrees)")
plt.ylabel("Focus Offset")
plt.subplot(1,2,2)
plt.scatter(els, poss)
plt.xlabel("Elevation(Degrees)")
plt.ylabel("Hexapod Z")
#plt.savefig(f"/project/cslage/AuxTel/efd_temp/Focus_vs_Elevation_{dayObs}_03Dec21.pdf")

```

```python

```
