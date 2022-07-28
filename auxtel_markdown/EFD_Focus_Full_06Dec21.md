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

# AuxTel Focus Study - 06-Dec-21

In this notebook, investigate focus settings and temp on observations between\
08-Sep21 and 04-Nov-21.

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

```python tags=[]
before = 15.0
temp_before = 60.0


els = []
offs = []
poss = []
temp_airs = []
temp_trusss = []
temp_m2s = []
temp_exts = []
filt = []
user = []

filename = '/project/cslage/AuxTel/efd_temp/EFD_Temp_Full_08Dec21.txt'
outfile = open(filename, 'w')
outfile.write(f"expId\t\tElev\t\tOff_disp\tOff_filt\tOff_user\tOff_tot\t\tHex_z\t\tT_air\tT_truss\tT_M2\tT_ext\tW_spd\tW_dir\tDIMM\n")

dayObss = [20210908, 20210909, 20211005, 20211006, 20211102, 20211103, 20211104]
seqNos = [[128, 134, 145, 149, 153, 161, 165, 489, 614, 641, 793], [152, 243, 348, 470, 542, 674, 773, 800],\
         [ 297, 302, 307, 310, 316, 398, 415, 422, 662], [146, 545, 552], [93, 333, 346, 351, 374, 377, 383, 399, 498, 567],\
         [70, 161, 176, 288, 446, 551, 621], [186, 253, 280, 298, 362, 383, 412, 581, 613, 942, 949, 960]]
for n, dayObs in enumerate(dayObss):
    if dayObs < 20211013:
        # EFD was switched to UTC on 20211013
        tai_offset = 37.0
    else:
        tai_offset = 0.0
        
    for seqNo in seqNos[n]:
        expId = dayObs * 100000 + seqNo
        mData = butler.get('raw.metadata', detector=0, exposure=expId)
        date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
        elevation = mData['ELSTART']
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
        z_position = await client.select_time_series("lsst.sal.ATHexapod.command_moveToPosition", "z", start, end)
        z_position = z_position.values[-1][0]
        
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
            dimm_fwhm = 999.99
        # Offsets are corrupted about half the time.  Discard the bad ones
        offset_error = total_off - (disp_off+filter_off+user_off)
        if abs (offset_error) > 1.0E-6:
            continue
        else:
            els.append(elevation)
            filt.append(filter_off)
            user.append(user_off)
            offs.append(total_off)
            poss.append(z_position)
            temp_airs.append(temp_air)
            temp_trusss.append(temp_truss)
            temp_m2s.append(temp_m2)
            temp_exts.append(temp_ext)
            outfile.write(f"{expId}\t{elevation:.4f}\t\t{disp_off:.6f}\t{filter_off:.6f}\t{user_off:.6f}\t{total_off:.6f}\t{z_position:.6f}\t{temp_air:.2f}\t{temp_truss:.2f}\t{temp_m2:.2f}\t{temp_ext:.2f}\t{wind_spd:.2f}\t{wind_dir:.2f}\t{dimm_fwhm:.2f}\n")
        
        #print(expId, elevation, total_off, z_position, temp_air, temp_truss, temp_m2)
        #print(disp_off, filter_off, user_off, temp_ext, wind_spd, wind_dir, dimm_fwhm)
outfile.close()
els = np.array(els)
filt = np.array(filt)
user = np.array(user)
offs = np.array(offs)
poss = np.array(poss)
temp_airs = np.array(temp_airs)
temp_trusss = np.array(temp_trusss)
temp_m2s = np.array(temp_m2s)
temp_exts = np.array(temp_exts)
```

```python
fig = plt.figure(figsize = (12,6))
plt.suptitle(f"Focus vs Elevation - 20210908 - 20211104", fontsize = 18)
plt.subplots_adjust(wspace=0.5)
plt.subplot(1,3,1)
plt.scatter(els, poss)
plt.xlabel("Elevation(Degrees)")
plt.ylabel("Hexapod Z")
plt.ylim(-1.30, -1.10)
plt.xlim(20.0, 90.0)

plt.subplot(1,3,2)
plt.scatter(els, poss-filt)
plt.xlabel("Elevation(Degrees)")
plt.ylabel("Hexapod Z - Filter offset")
plt.ylim(-1.30, -1.10)
plt.xlim(20.0, 90.0)

# Plot the focus ve elevation model
poly_z = -0.217
poly_int = -1.045
xplot = np.linspace(5.5, 89.0, 100)
yplot = poly_int + poly_z * np.cos(np.radians(90.0 - xplot))
plt.subplot(1,3,3)
plt.scatter(els, poss-filt-user)
plt.plot(xplot, yplot, ls = '--', color = 'red', label = 'Model')
plt.xlabel("Elevation(Degrees)")
plt.ylabel("Hexapod Z - Filter offset - User_applied offset")
plt.ylim(-1.30, -1.10)
plt.xlim(20.0, 90.0)
plt.legend()

#plt.savefig(f"/project/cslage/AuxTel/efd_temp/Focus_vs_Elevation_1_08Dec21.pdf")

```

```python
from scipy.optimize import curve_fit

def func(X, intercept, slopeE, slopeT):
    el, T = X
    return (intercept + slopeE * el + slopeT * T)

def func_cos(X, intercept, slopeE, slopeT):
    el, T = X
    return (intercept + slopeE * np.cos(np.radians(90.0 - el)) + slopeT * T)


def chiSquared(myFunc, els, temps, vals, sigma):
    # Chi squared per DOF for the given model
    chiSq = 0.0
    intercept, slopeE, slopeT = fit
    model = myFunc((els, temps), intercept, slopeE, slopeT)
    for i, val in enumerate(vals):
        chiSq += (model[i] - val)**2 / sigma**2
    # Degrees of freedom are number of data points - number of fit parameters
    dof = len(vals) - 3
    return chiSq / dof
        
```

```python
names = ['Air', 'Truss', 'M2', 'Ext']
for i, temps in enumerate([temp_airs, temp_trusss, temp_m2s, temp_exts]):
    p0 = -1.2, -2E-3, 5E-3
    fit, covs = curve_fit(func, (els, temps), poss - filt, p0)
    intercept, slopeE, slopeT = fit
    #print(fit)
    print(f"Linear fit, {names[i]}: ChiSquared = {chiSquared(func, els, temps, poss - filt, 0.01)}")
    fit, covs = curve_fit(func_cos, (els, temps), poss - filt, p0)
    intercept, slopeE, slopeT = fit
    #print(fit)
    print(f"Cosine fit, {names[i]}: ChiSquared = {chiSquared(func_cos, els, temps, poss - filt, 0.01)}")
```

```python
temps = temp_airs
fit, covs = curve_fit(func, (els, temps), poss - filt, p0)
intercept, slopeE, slopeT = fit

fig = plt.figure(figsize = (8,6))
plt.suptitle(f"Focus vs Elevation - 20210908 - 20211104", fontsize = 18)
plt.subplots_adjust(wspace=0.5)
plt.subplot(1,2,1)
plt.scatter(els, poss - filt)
plt.xlabel("Elevation(Degrees)")
plt.ylabel("Hexapod Z - Filter offset")
plt.ylim(-1.30, -1.10)
plt.xlim(20.0, 90.0)
plt.subplot(1,2,2)
plt.scatter(els, poss - filt - slopeT * temps)
xplot = np.linspace(0,90,100)
yplot = intercept + slopeE * xplot
plt.plot(xplot, yplot, ls = '--', color = 'red')
plt.ylim(-1.30, -1.10)
plt.xlim(20.0, 90.0)
plt.xlabel("Elevation(Degrees)")
plt.ylabel("Hexapod Z - Filter offset - slopeT * T_air")
plt.text(25.0, -1.125, f"Intercept = {intercept:0.3f} mm")
plt.text(25.0, -1.1375, f"SlopeT = {slopeT:0.6f} mm/deg C")
plt.text(25.0, -1.15, f"SlopeE = {slopeE:0.6f} mm/deg El")
#plt.savefig(f"/project/cslage/AuxTel/efd_temp/Focus_vs_Elevation_5_08Dec21.pdf")

```

```python
fig = plt.figure(figsize = (12,12))
plt.suptitle(f"Focus vs Elevation - 20210908 - 20211104", fontsize = 18)
plt.subplots_adjust(wspace=0.3, hspace=0.5)

for i, temps in enumerate([temp_airs, temp_trusss, temp_m2s, temp_exts]):
    p0 = -1.2, -2E-3, 5E-3
    fit, covs = curve_fit(func, (els, temps), poss - filt, p0)
    intercept, slopeE, slopeT = fit
    chi2 = chiSquared(func, els, temps, poss - filt, 0.01)
    
    plt.subplot(4,2,2*i+1)
    plt.title(f"{names[i]}, Linear", fontsize=12)
    plt.scatter(els, poss - filt - slopeT * temps)
    xplot = np.linspace(0,90,100)
    yplot = intercept + slopeE * xplot
    plt.plot(xplot, yplot, ls = '--', color = 'red')
    plt.ylim(-1.30, -1.10)
    plt.xlim(20.0, 90.0)
    plt.xlabel("Elevation(Degrees)")
    plt.ylabel("HexZ")
    plt.text(40.0, -1.12, f"Intercept = {intercept:0.3f} mm")
    plt.text(40.0, -1.14, f"SlopeT = {slopeT:0.6f} mm/deg C")
    plt.text(40.0, -1.16, f"SlopeE = {slopeE:0.6f} mm/deg El")
    plt.text(50.0, -1.18, f"Chi^2/DOF = {chi2:0.2f}")

    fit, covs = curve_fit(func_cos, (els, temps), poss - filt, p0)
    intercept, slopeE, slopeT = fit
    chi2 = chiSquared(func_cos, els, temps, poss - filt, 0.01)
    plt.subplot(4,2,2*i+2)
    plt.title(f"{names[i]}, Cosine", fontsize=12)
    plt.scatter(els, poss - filt - slopeT * temps)
    xplot = np.linspace(0,90,100)
    yplot = intercept + slopeE * np.cos(np.radians(90.0 - xplot))
    plt.plot(xplot, yplot, ls = '--', color = 'red')
    plt.ylim(-1.30, -1.10)
    plt.xlim(20.0, 90.0)
    plt.xlabel("Elevation(Degrees)")
    plt.ylabel("HexZ")
    plt.text(40.0, -1.12, f"Intercept = {intercept:0.3f} mm")
    plt.text(40.0, -1.14, f"SlopeT = {slopeT:0.6f} mm/deg C")
    plt.text(40.0, -1.16, f"SlopeE = {slopeE:0.6f} mm/deg El")
    plt.text(50.0, -1.18, f"Chi^2/DOF = {chi2:0.2f}")
plt.savefig(f"/project/cslage/AuxTel/efd_temp/Focus_vs_Elevation_Models_13Dec21.pdf")
```

```python

```
