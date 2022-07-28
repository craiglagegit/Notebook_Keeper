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

## AuxTel AzEl offsets - 15-Apr-21

In this notebook, investigate az-el offsets from 11-Mar-21

```python
import sys, time, os, asyncio, glob

from datetime import datetime
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pickle as pkl
import pandas as pd
import astropy.io.fits as pf
from astropy.time import Time, TimeDelta
from astropy.time import Time, TimeDelta
from astropy.coordinates import SkyCoord, AltAz, ICRS, EarthLocation, Angle, FK5
import astropy.units as u

from lsst.daf.butler import Butler as gen3Butler
from lsst.daf.persistence import Butler as gen2Butler
from lsst_efd_client import EfdClient
from lsst.pipe.tasks.characterizeImage import CharacterizeImageTask, CharacterizeImageConfig
```

```python jupyter={"outputs_hidden": true}
# Gen3 butler
dayObs = 20210311
myExpIds = [2021031100259, 2021031100260, 2021031100261, 2021031100281, 2021031100282, 2021031100283]
REPO_DIR = '/repo/main'
butler = gen3Butler(REPO_DIR, collections="LATISS/raw/all")

exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.day_obs=%d"%dayObs):
    exposureList.append(record.id)
exposureList.sort()
myVisits = []
for exposure in exposureList:
    mData = butler.get('raw.metadata', detector=0, exposure=exposure)
    expTime = mData['EXPTIME']
    imgType = mData['IMGTYPE']
    obj = mData['OBJECT']
    filter = mData['FILTER']
    rotpa = mData['ROTPA']
    date_beg = mData['DATE-BEG']
    elstart = mData['ELSTART']
    azstart = mData['AZSTART']
    rastart = mData['RASTART']
    decstart = mData['DECSTART']
    dummy=0.0
    # Need to use DATE-BEG to get the right timestamp
    visit = (exposure, expTime, imgType, obj, filter, date_beg, rotpa, dummy, dummy, azstart, elstart, rastart, decstart)
    if (int(visit[0]) in myExpIds):
        myVisits.append(visit)
        print(visit)
```

```python jupyter={"outputs_hidden": true}
mData = butler.get('raw.metadata', detector=0, exposure=2021031100261)
for key in mData.keys():
    print(key, mData[key])
```

```python
for visit in myVisits:
    print(visit)
```

```python
# Get EFD client
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
# These are for finding the timestamps of the offset events
backUp = 120 # seconds before first image to get initial offset
start = Time(myVisits[0][5],scale='tai') - TimeDelta(backUp, format='sec')
end = Time(myVisits[-1][5],scale='tai') - TimeDelta(0, format='sec')
timestamp = f"time >= {start} AND time <= {end}"
```

```python
# Now get the offsets applied and the nasmyth angle
offsets = await client.select_time_series("lsst.sal.ATPtg.command_offsetAzEl", ['*'],
                                          start, end)


```

```python
print(len(offsets))
```

```python
for i, offset in enumerate(offsets.values):
    print(i, Time(offsets.index[i]).tai.isot,offset[0], offset[1])
    if i > 8:
        break           
```

```python jupyter={"outputs_hidden": true}
expId = 2021031100261
mData = butler.get('raw.metadata', detector=0, exposure=expId)
rotpa = mData['ROTPA']
date_beg = mData['DATE-BEG']
date_end = mData['DATE-END']
el = mData['ELSTART']
backup = 0.0
start = Time(date_beg,scale='tai') - TimeDelta(backUp, format='sec')
end = Time(date_end,scale='tai')
nasmyth_position = await client.select_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", ['*'],
                                              start, end)

rot = merge_packed_time_series(nasmyth_position, 'nasmyth2CalculatedAngle', stride=1)

```

```python
print(rot.values[0][0], rot.values[-1][0])
```

```python
bore_sight_angle = el - rot.values[0][0] + 90.0
print(bore_sight_angle)
#off = np.array([126.0, 55.0])
off = np.array([140.0, 0.0])
theta = Angle(bore_sight_angle*u.deg).rad 
c, s = np.cos(theta), np.sin(theta)
R = np.array(((c, -s), (s, c))) 
rotated_off = R.dot(off)
print(rotated_off)
print(np.arctan2(12.5958, -139.432) * 180.0 / np.pi)
```

```python
# Plot the first few to check the interleaving of the offsets with the exposures
# Blue are the times of setting the offsets, and red are the start of the exposure
startPlot = Time('2021-03-12T02:28:10',scale='tai') 
endPlot = Time('2021-03-12T02:28:55',scale='tai')
mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                          startPlot, endPlot)
time.sleep(2.0)
az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1)
el = merge_packed_time_series(mount_position, 'elevationCalculatedAngle', stride=1)
```

```python
np.array(el.values.tolist()).shape
```

```python
# Calculate the tracking errors
az_vals = np.array(az.values.tolist())[:,0]
el_vals = np.array(el.values.tolist())[:,0]
times = np.array(az.values.tolist())[:,1]
times = times - times [0]

az_vals_fit = np.array(az.values.tolist())[0:2000,0]
el_vals_fit = np.array(el.values.tolist())[0:2000,0]
times_fit = np.array(az.values.tolist())[0:2000,1]
times_fit = times_fit - times_fit [0]

# Fit with a quadratic
az_fit = np.polyfit(times_fit, az_vals_fit, 2)
el_fit = np.polyfit(times_fit, el_vals_fit, 2)

az_model = az_fit[0] * times * times + az_fit[1] * times + az_fit[2]
el_model = el_fit[0] * times * times + el_fit[1] * times + el_fit[2]

# Errors in arcseconds
az_error = (az_vals - az_model) * 3600
el_error = (el_vals - el_model) * 3600

az_shift = np.mean(az_error[3200:4000])
el_shift = np.mean(el_error[3200:4000])

```

```python
fig = plt.figure(figsize=(8,16))
plt.subplots_adjust(hspace=1.0)
plt.subplot(3,1,1)
plt.title(f"Commands - 11-Mar-21", fontsize = 18)
# Azimuth axis
ax1 = offsets['num'].plot(color='red')
ax1.set_ylim(0,1.0)
for i in range(6):
    try:
        t1 = Time(offsets.index[i]).tai.isot
        ax1.axvline(t1, ymin=0.5, ymax=0.9, color="blue")
    except:
        pass
    t2 = Time(myVisits[i][5]).tai.isot
    ax1.axvline(t2, ymin=0.1, ymax=0.5, color="red")
ax1.set_xlim(startPlot.tai.isot,endPlot.tai.isot)
plt.subplot(3,1,2)
plt.title(f"Azimuth change", fontsize = 18)
plt.plot(times, az_error, color='green')
plt.text(30.0,5.0, f"AZ_shift = {az_shift:.4f} arcsec")
plt.subplot(3,1,3)
plt.title(f"Elevation change", fontsize = 18)
plt.plot(times, el_error, color='green')
plt.text(30.0,5.0, f"EL_shift = {el_shift:.4f} arcsec")


```

```python
# Now append the applied offsets to the list of visits
# A few drop out because the offsets are not clear
backUp = 240
fullVisits = []
for i, visit in enumerate(myVisits):
    newList = list(visit)
    if i == 0:
        startTime = Time(myVisits[i][5],scale='tai',precision=0) - TimeDelta(backUp, format='sec')
        startTime = startTime.tai.isot
    elif i==3 or i==4:
        startTime = Time(myVisits[1][5],scale='tai',precision=0).tai.isot
    else:
        startTime = Time(myVisits[i - 1][5],scale='tai',precision=0).tai.isot
    endTime = Time(myVisits[i][5],scale='tai',precision=0).tai.isot
    #print(startTime, endTime)
    try:
        offset = offsets.loc[startTime:endTime].values
        if len(offset) == 1:

            newList[7] = offset[0][0]
            newList[8] = offset[0][1]
            
        else:
            print("Not = 1", len(offset))
            pass
    except:
        print("Failed the try")
        pass
    print(newList[0], newList[7], newList[8])
    fullVisits.append(newList)
```

```python
for fullVisit in fullVisits:
    print(fullVisit[0],fullVisit[7],fullVisit[8])

```

```python jupyter={"outputs_hidden": true}
# Get the raw quickLook data.  Only Gen2 works
REPO_DIR = '/project/shared/auxTel/rerun/quickLook'
gen2_butler = gen2Butler(REPO_DIR)
dayObs = '2021-03-11'
```

```python
charConfig = CharacterizeImageConfig()
charConfig.doMeasurePsf = False#True
charConfig.doApCorr = False
charConfig.doDeblend = False
charConfig.repair.doCosmicRay = True
charConfig.repair.doInterpolate = True   
charConfig.detection.minPixels = 500
charTask = CharacterizeImageTask(config=charConfig)
```

```python jupyter={"outputs_hidden": true}
# Try doing them all
charVisits = []
for fullVisit in fullVisits:
    expId = fullVisit[0]
    try:
        charVisit = {}
        charVisit['Visit'] = fullVisit
        exp = gen2_butler.get('quickLookExp', detector=0, expId=expId)
        charResult = charTask.run(exp)
        sourceCatalog = charResult.sourceCat
        #print(expId)
        #print(sourceCatalog['base_SdssShape_instFlux'])
        maxFlux = np.nanmax(sourceCatalog['base_CircularApertureFlux_3_0_instFlux'])
        selectBrightestSource = sourceCatalog['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.99
        brightestSource = sourceCatalog.subset(selectBrightestSource)
        brightestCentroid = (brightestSource['base_SdssCentroid_x'][0], \
                             brightestSource['base_SdssCentroid_y'][0])
        brightCatalog = sourceCatalog.subset(sourceCatalog['base_CircularApertureFlux_3_0_instFlux'] > maxFlux * 0.001)
        print(f"expId:{expId}. Found {len(sourceCatalog)} sources, {len(brightCatalog)} bright sources")
        print(f"Brightest centroid at {brightestCentroid}")
        charVisit['exp'] = exp
        charVisit['brightestCentroid'] = brightestCentroid
        charVisit['brightCatalog'] = sourceCatalog#brightCatalog
        charVisits.append(charVisit)
    except:
        print(f"Skipping expId {expId}.")
        continue
```

```python
for charVisit in charVisits:
    print(charVisit['Visit'][0], charVisit['Visit'][7], charVisit['Visit'][8])
    print(charVisit['brightestCentroid'])
```

```python
outfile = open('/project/cslage/AuxTel/offsets/offsets_HD75519_27apr21.pkl','wb')

pkl.dump(charVisits,outfile)
outfile.close()
```

```python
infile = open('/project/cslage/AuxTel/offsets/offsets_HD75519_27apr21.pkl','rb')
charVisits = pkl.load(infile)
infile.close()
```

```python
print([charVisits[1]['brightestCentroid'][0] - charVisits[0]['brightestCentroid'][0], \
    charVisits[1]['brightestCentroid'][1] - charVisits[0]['brightestCentroid'][1]], \
      np.array([charVisits[1]['Visit'][7] - charVisits[0]['Visit'][7], \
    charVisits[1]['Visit'][8] - charVisits[0]['Visit'][8]]) / 0.095)
```

```python
print([charVisits[2]['brightestCentroid'][0] - charVisits[1]['brightestCentroid'][0], \
    charVisits[2]['brightestCentroid'][1] - charVisits[1]['brightestCentroid'][1]], \
      np.array([charVisits[2]['Visit'][7] - charVisits[1]['Visit'][7], \
    charVisits[2]['Visit'][8] - charVisits[1]['Visit'][8]]) / 0.095)
```

```python
print([charVisits[2]['brightestCentroid'][0] - charVisits[1]['brightestCentroid'][0], \
    charVisits[2]['brightestCentroid'][1] - charVisits[1]['brightestCentroid'][1]], \
      np.array([charVisits[2]['Visit'][7], \
    charVisits[2]['Visit'][8]]) / 0.095)
```

```python
print([charVisits[5]['brightestCentroid'][0] - charVisits[4]['brightestCentroid'][0], \
    charVisits[5]['brightestCentroid'][1] - charVisits[4]['brightestCentroid'][1]], \
      np.array([charVisits[5]['Visit'][7] - charVisits[4]['Visit'][7], \
    charVisits[5]['Visit'][8] - charVisits[4]['Visit'][8]]) / 0.095)
```

```python
print([charVisits[5]['brightestCentroid'][0] - charVisits[4]['brightestCentroid'][0], \
    charVisits[5]['brightestCentroid'][1] - charVisits[4]['brightestCentroid'][1]], \
      np.array([charVisits[5]['Visit'][7], \
    charVisits[5]['Visit'][8]]) / 0.095)
```

```python
print([charVisits[4]['brightestCentroid'][0] - charVisits[3]['brightestCentroid'][0], \
    charVisits[4]['brightestCentroid'][1] - charVisits[3]['brightestCentroid'][1]], \
      np.array([charVisits[4]['Visit'][7] - charVisits[3]['Visit'][7], \
    charVisits[4]['Visit'][8] - charVisits[3]['Visit'][8]]) / 0.095)
```

```python
#%matplotlib inline
# Look at the data with matplotlib
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


#Plot these four
expIds = [2021031100260, 2021031100261, 2021031100282, 2021031100283]

plt.figure(figsize=(16,16))
plotCounter = 1
for charVisit in charVisits:
    expId = charVisit['Visit'][0]
    offset = (charVisit['Visit'][7],charVisit['Visit'][8])
    if expId not in expIds:
        continue
    plt.subplot(2,2,plotCounter)
    plotCounter += 1
    plt.title(f"Image - {expId}",fontsize=18)
    arr = charVisit['exp'].image.array
    arr = np.clip(arr, 1, 100000) # This image has some negative values, and this removes them
    img = plt.imshow(arr, norm=LogNorm(vmin=1, vmax=1000),  interpolation='Nearest', cmap='gray')
    cat = charVisit['brightCatalog']
    plt.scatter(cat['base_SdssCentroid_x'],cat['base_SdssCentroid_y']\
                ,color='red', marker='x', label=f"offset={offset}")
    plt.scatter([charVisit['brightestCentroid'][0]],[charVisit['brightestCentroid'][1]] \
                ,color='green', marker='+', s=100)
    colorbar(img)
    plt.legend()
    plt.ylim(0,4000)
plt.tight_layout(h_pad=1)
#plt.savefig(f"/project/cslage/AuxTel/offsets/Offsets_{expIds[0]}_{expIds[1]}_{expIds[2]}_{expIds[3]}_16Apr21.pdf")
```

```python
np.arctan2(12.5958,-139.43222)*180.0/np.pi
```

```python

```

```python
# Plot the first few to check the interleaving of the offsets with the exposures
# Blue are the times of setting the offsets, and red are the start of the exposure
startPlot = Time('2021-03-12T04:36:36',scale='tai') 
endPlot = Time('2021-03-12T04:37:00',scale='tai')
mount_position = await client.select_time_series("lsst.sal.ATMCS.mount_AzEl_Encoders", ['*'],
                                          startPlot, endPlot)
time.sleep(2.0)
az = merge_packed_time_series(mount_position, 'azimuthCalculatedAngle', stride=1)
el = merge_packed_time_series(mount_position, 'elevationCalculatedAngle', stride=1)
```

```python
np.array(el.values.tolist()).shape
```

```python
# Calculate the tracking errors
n1=0
n2=1000
n3=2000
n4=2200


az_vals = np.array(az.values.tolist())[:,0]
el_vals = np.array(el.values.tolist())[:,0]
times = np.array(az.values.tolist())[:,1]
times = times# - times [0]

az_vals_fit = np.array(az.values.tolist())[n1:n2,0]
el_vals_fit = np.array(el.values.tolist())[n1:n2,0]
times_fit = np.array(az.values.tolist())[n1:n2,1]
times_fit = times_fit# - times_fit [0]

# Fit with a quadratic
az_fit = np.polyfit(times_fit, az_vals_fit, 2)
el_fit = np.polyfit(times_fit, el_vals_fit, 2)

az_model = az_fit[0] * times * times + az_fit[1] * times + az_fit[2]
el_model = el_fit[0] * times * times + el_fit[1] * times + el_fit[2]

# Errors in arcseconds
az_error = (az_vals - az_model) * 3600
el_error = (el_vals - el_model) * 3600

az_shift = np.mean(az_error[n3:n4])
el_shift = np.mean(el_error[n3:n4])

```

```python
fig = plt.figure(figsize=(8,16))
plt.subplots_adjust(hspace=1.0)
plt.subplot(3,1,1)
plt.title(f"Commands - 11-Mar-21", fontsize = 18)
# Azimuth axis
ax1 = offsets['num'].plot(color='red')
ax1.set_ylim(0,1.0)
for i in range(6):
    try:
        t1 = Time(offsets.index[i]).tai.isot
        ax1.axvline(t1, ymin=0.5, ymax=0.9, color="blue")
    except:
        pass
    t2 = Time(myVisits[i][5]).tai.isot
    ax1.axvline(t2, ymin=0.1, ymax=0.5, color="red")
ax1.set_xlim(startPlot.tai.isot,endPlot.tai.isot)
plt.subplot(3,1,2)
plt.title(f"Azimuth change", fontsize = 18)
#ax2 = az['azimuthCalculatedAngle'].plot(legend=False, color='green')
plt.plot(times, az_error, color='green')
plt.text(times[n2], 20.0, f"AZ_shift = {az_shift:.4f} arcsec")
plt.subplot(3,1,3)
plt.title(f"Elevation change", fontsize = 18)
#ax3 = el['elevationCalculatedAngle'].plot(legend=False, color='green')
plt.plot(times, el_error, color='green')
plt.text(times[n2], -80.0, f"EL_shift = {el_shift:.4f} arcsec")


```

```python

```
