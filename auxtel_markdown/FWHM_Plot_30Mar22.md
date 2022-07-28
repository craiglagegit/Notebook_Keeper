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

## AuxTel Image elongation due to azimuth oscillation

Craig Lage 17-Mar-22

```python
import sys, time, os, asyncio, glob
from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from scipy.stats import linregress
import pickle as pkl
from astropy.time import Time, TimeDelta
from astropy.coordinates import AltAz, ICRS, EarthLocation, Angle, FK5
import astropy.units as u
from lsst.obs.lsst.translators.latiss import AUXTEL_LOCATION
from lsst.daf.butler import Butler
from lsst_efd_client import EfdClient
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python
# Get EFD client and the butler
client = EfdClient('ldf_stable_efd')
butler = Butler('/repo/main', collections="LATISS/runs/quickLook")
```

```python
# Gen3 butler - get the observations from a given night
dayObs = '2022-04-07'
dayObs = int(dayObs.replace('-', ''))

exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
    exposureList.append([record.id, record])
exposureList.sort(key=lambda x: x[0])
FWHM_list = []
for [id,record] in exposureList:
    #print(record.id, record.observation_type, record.exposure_time, record.physical_filter, record.target_name)
    type = record.observation_type
    grating = record.physical_filter.split('~')[1]
    if type == 'science' and grating == 'empty':
        FWHM_list.append(record.id)
print(len(FWHM_list), "images will be characterized")
```

```python
# Set up the source catalog task
charConfig = CharacterizeImageConfig()
charConfig.doMeasurePsf = False
charConfig.doApCorr = False
charConfig.doDeblend = False
charConfig.repair.doCosmicRay = True
charConfig.repair.doInterpolate = True   
charConfig.detection.minPixels = 500
charTask = CharacterizeImageTask(config=charConfig)
```

```python
def RotatedMoments(Ixx, Iyy, Ixy, theta):
    # Rotates the moments about an angle theta.
    # Formulae are from the Sextractor documentation
    # https://sextractor.readthedocs.io/en/latest/Position.html\
    # ?highlight=shape#basic-shape-parameters-a-b-theta
    c = np.cos(theta)
    s = np.sin(theta)
    IxxRot = c * c * Ixx + s * s * Iyy - 2.0 * c * s * Ixy
    IyyRot = s * s * Ixx + c * c * Iyy + 2.0 * c * s * Ixy
    IxyRot = c * s * (Ixx - Iyy) + (c * c - s * s) * Ixy
    return [IxxRot, IyyRot, IxyRot] 
```

```python tags=[]
# Now get the image data and calculate the median FWHM
data = {}
for i, expId in enumerate(FWHM_list):
    #if i%10 != 0:
    #    continue
    try:
        expData = {}
        exp = butler.get('quickLookExp', detector=0, exposure=expId)
        mData = exp.getMetadata()
        charResult = charTask.run(exp)
        sourceCatalog = charResult.sourceCat
        rotpa = Angle(mData['ROTPA'] * u.deg)
        el = Angle(mData['ELSTART'] * u.deg)
        az = Angle(mData['AZSTART'] * u.deg)
        dec = Angle(mData['DECSTART'] * u.deg)
        lat = AUXTEL_LOCATION.lat
        # This calculates the angle theta between (X,Y) and (Az,El)
        sinTheta =  np.cos(lat) / np.cos(dec) * np.sin(az)
        cosTheta = (np.sin(el) * np.sin(dec) - np.sin(lat)) / (np.cos(el) * np.cos(dec))
        theta = Angle(np.arcsin(sinTheta))
        # The following removes the ambiguity in arcsin(theta)
        if cosTheta > 0:
            rotAzEl = rotpa - theta
        else:    
            rotAzEl = rotpa - (Angle(180.0 * u.deg) - theta)
        date_beg = Time(mData['DATE-BEG'], format='isot', scale='tai')
        time = date_beg.utc
        Ixx = sourceCatalog.getIxx()
        Iyy = sourceCatalog.getIyy()
        Ixy = sourceCatalog.getIxy()
        Ip = (Ixx + Iyy) / 2.0
        Im = (Ixx - Iyy) / 2.0
        A2 = Ip + np.sqrt(Im**2 + Ixy**2)
        B2 = Ip - np.sqrt(Im**2 + Ixy**2)
        [Iaa, Iee, Iae] = RotatedMoments(Ixx, Iyy, Ixy, rotAzEl)
        FWHM_x = 2.35 * np.sqrt(Ixx) 
        FWHM_y = 2.35 * np.sqrt(Iyy) 
        FWHM_az = 2.35 * np.sqrt(Iaa) 
        FWHM_el = 2.35 * np.sqrt(Iee) 
        FWHM_a = 2.35 * np.sqrt(A2) 
        FWHM_b = 2.35 * np.sqrt(B2) 
        expData['Time'] = time.isot
        expData['FWHM_x'] = np.median(FWHM_x)
        expData['FWHM_y'] = np.median(FWHM_y)
        expData['FWHM_az'] = np.median(FWHM_az)
        expData['FWHM_el'] = np.median(FWHM_el)
        expData['FWHM_a'] = np.median(FWHM_a)
        expData['FWHM_b'] = np.median(FWHM_b)
        data[expId] = expData
    except:
        continue
outfile = open(f'/project/cslage/AuxTel/fwhm/FWHM_New_{dayObs}.pkl', 'wb')
pkl.dump(data,outfile)
outfile.close()
```

```python
# Now unpickle the data and plot against the DIMM data
```

```python
dayObs = 20220316
infile = open(f'/project/cslage/AuxTel/fwhm/FWHM_New_{dayObs}.pkl','rb')
data = pkl.load(infile)
infile.close()
print(len(data), "observations")
```

```python
times = []
fwhmx = []
fwhmy = []
fwhmaz = []
fwhmel = []
fwhma = []
fwhmb = []
for expId in data.keys():
    times.append(Time(data[expId]['Time']).to_datetime())
    fwhmx.append(data[expId]['FWHM_x'])
    fwhmy.append(data[expId]['FWHM_y'])
    fwhmaz.append(data[expId]['FWHM_az'])
    fwhmel.append(data[expId]['FWHM_el'])
    fwhma.append(data[expId]['FWHM_a'])
    fwhmb.append(data[expId]['FWHM_b'])

plt.scatter(times, fwhmb)     
```

```python
# Now get the DIMM data
tstart = Time(times[0])
tend = Time(times[-1])
dimm_fwhm = await client.select_time_series("lsst.sal.DIMM.logevent_dimmMeasurement", \
                                            "fwhm", tstart, tend)
```

```python
# Now plot it all
fig = plt.figure(figsize=(8,10))
plt.subplots_adjust(hspace=0.5)
plt.suptitle(f"AuxTel FWHM -{dayObs}", fontsize = 16)
# X and Y
plt.subplot(3,1,1)
ax1 = dimm_fwhm['fwhm'].plot(color = 'lime', label='DIMM', lw=3)
ax1.scatter(times, np.array(fwhmx) / 10.0, color='red', marker = 'x', s=10, label = 'FWHM_x')
ax1.scatter(times, np.array(fwhmy) / 10.0, color='blue', marker = 'x', s=10, label = 'FWHM_y')
ax1.set_ylabel('FWHM(arcseconds)')
ax1.set_ylim(0.0, 2.0)
ax1.legend(loc='lower left')
# Az and El
plt.subplot(3,1,2)
ax2 = dimm_fwhm['fwhm'].plot(color = 'lime', label='DIMM', lw=3)
ax2.scatter(times, np.array(fwhmaz) / 10.0, color='red', marker = 'x', s=10, label = 'FWHM_az')
ax2.scatter(times, np.array(fwhmel) / 10.0, color='blue', marker = 'x', s=10, label = 'FWHM_el')
ax2.set_ylabel('FWHM(arcseconds)')
ax2.set_ylim(0.0, 2.0)
#ax2.text(times[int(len(times) / 2)], 0.2, "NB - I might have Az and El reversed!")
ax2.legend(loc='lower left')
# a and b
plt.subplot(3,1,3)
ax3 = dimm_fwhm['fwhm'].plot(color = 'lime', label='DIMM', lw=3)
ax3.scatter(times, np.array(fwhma) / 10.0, color='red', marker = 'x', s=10, label = 'FWHM_a')
ax3.scatter(times, np.array(fwhmb) / 10.0, color='blue', marker = 'x', s=10, label = 'FWHM_b')
ax3.set_ylabel('FWHM(arcseconds)')
ax3.set_ylim(0.0, 2.0)
ax3.legend(loc='lower left')
plt.savefig(f'/project/cslage/AuxTel/fwhm/FWHM_New_{dayObs}.pdf')
```

```python
# Now plot a scatter plot of three nights of observations
#filenames = ['FWHM_New_20220405.pkl', 'FWHM_New_20220406.pkl', 'FWHM_New_20220407.pkl']
filenames = ['FWHM_New_20220405.pkl', 'FWHM_New_20220406.pkl']
#filenames = ['FWHM_New_20220316.pkl']
auxtel_fwhm = []
dimm = []
for filename in filenames:
    infile = open('/project/cslage/AuxTel/fwhm/%s'%filename,'rb')
    data = pkl.load(infile)
    infile.close()
    print(filename, len(data), "observations")
    times = []
    fwhm = []
    for expId in data.keys():
        times.append(Time(data[expId]['Time']).to_datetime())
        fwhm.append((data[expId]['FWHM_a'] + data[expId]['FWHM_b']) / 20.0)

    # Get the DIMM data
    tstart = Time(times[0])
    tend = Time(times[-1])
    dimm_fwhm = await client.select_time_series("lsst.sal.DIMM.logevent_dimmMeasurement", \
                                            "fwhm", tstart, tend)
    
    # Skip measurements early in the night
    for i, time in enumerate(times):
        if time.hour < 3 or time.hour > 14:
            continue
        # Now get the DIMM measurement closest in time
        nearest_dimm_index = dimm_fwhm.index.get_loc(time, method='nearest')
        dimm_measure = dimm_fwhm['fwhm'][nearest_dimm_index]
        if np.isnan(dimm_measure) or np.isnan(fwhm[i]):
            continue
        dimm.append(dimm_fwhm['fwhm'][nearest_dimm_index])
        auxtel_fwhm.append(fwhm[i])
                    
print(len(dimm), "observations")    
```

```python
# Now plot the scatter plot
slope, intercept, r_value, p_value, std_err = linregress(dimm, auxtel_fwhm)
xplot = np.linspace(0.5, 2.0, 100)
yplot = intercept + slope * xplot
#plt.title("FWHM Measurements, April 5,6,7, 2022")
plt.title("FWHM Measurements, March 16, 2022")
plt.plot(xplot,yplot, color='red',ls = '--', lw=2)
plt.scatter(dimm, auxtel_fwhm)
plt.xlim(0.5, 2.0)
plt.ylim(0.5, 2.0)
plt.xlabel("Rubin DIMM (arcsec)")
plt.ylabel("AuxTel FWHM (a+b)/2 (arsec)")
plt.text(0.6,1.7,f"R2={r_value:.2f}")
plt.savefig(f'/project/cslage/AuxTel/fwhm/FWHM_April_5-6_2022.pdf')
```

```python

```
