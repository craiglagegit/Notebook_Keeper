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

# This notebook demonstrates a functionality issue with the select_packed_time_series of the efd-client

```python
import matplotlib
%matplotlib widget
import aioinflux
import getpass
import pandas as pd
import asyncio
import numpy as np
from astropy.time import Time, TimeDelta
from matplotlib import pyplot as plt
from bokeh.plotting import figure, output_notebook, show
from bokeh.models import LinearAxis, Range1d
output_notebook()
from bokeh.models import Span, Label

from lsst_efd_client import EfdClient, resample
```

```python
client = EfdClient('summit_efd')
client.output = 'dataframe'
```

### Declare timestamp(s) used for EFD queries

```python
# From an observeration last night where a fault occurred
t1 = Time('2022-02-16T05:15:00', scale='utc')
window = TimeDelta(40, format='sec')
t2=t1+window
```

```python
# From Craig's technote (https://sitcomtn-018.lsst.io/) on timestamps - the notebook is meant for the above dataset.
# If you run this one afterwards you'll see the same issue, but the time offset is different!
# https://github.com/craiglagegit/ScratchStuff/blob/master/notebooks/Plot_Tracking_UTC_29Oct21.ipynb

# t1 = Time('2021-10-14T18:48:30.004', scale='utc')
# window = TimeDelta(13, format='sec')
# t2=t1+window
```

The following are position (angle) data reported by the nasmyth rotator. <br>
The first is a packed dataset <br>
The second is a single value of the packed dataset, but pulled out individually and therefore does not use "select_packed_time_series"

```python
nas2_packed = await client.select_packed_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", ["nasmyth2CalculatedAngle"], t1, t2)
nas2_single = await client.select_time_series("lsst.sal.ATMCS.mount_Nasmyth_Encoders", ["nasmyth2CalculatedAngle0","nasmyth2CalculatedAngle99","private_sndStamp" ,"private_kafkaStamp", "cRIO_timestamp"], t1, t2)
```

Now, the packed data reports 100 values, where the first value corresponds to the `cRIO_timestamp` -- note that the cRIO_timestamp *IS* in TAI. <br>
So, the difference between cRIO_timestamp and private_sndStamp (also in TAI) is equal to how long it takes the data to get out of the FPGA and published to DDS. <br>
This means that the private_efdStamp will always be greater than cRIO_timestamp. <br>
Also, because there are 100 values published that represent 1s, it is expected that there might be an additional 1s delay<br>
We know this value ranges between 1-3 seconds or so and probably depends on the load to the FPGA.
Essentially, if there are errors in how the single and packed timestamps match up, it should be equal to this value!

```python
nas2_single['private_sndStamp']-nas2_single['cRIO_timestamp']
plt.hist(nas2_single['private_sndStamp']-nas2_single['cRIO_timestamp'])
plt.ylabel('Occurrences')
plt.xlabel('time difference [seconds]')
plt.show()
```

So that's all good! <br>
Now, let's confirm that the returned index is equal to private_kafkaStamp.
Things are a bit messy here because the index is dtype='datetime64[ns, UTC]'

```python
nas2_single_kafkaTime=Time(nas2_single['private_kafkaStamp'].values, format='unix_tai')
nas2_single_index_time=Time(nas2_single.index)
```

```python tags=[]
nas2_single_index_time-nas2_single_kafkaTime
```

So this is consistent to ~10 microseconds... I imagine that's good enough!


### Now look at the packed data
So 40 rows, should now be 4000 rows. <br>
Also the INDEX should be in UTC.

```python
print(nas2_packed.shape)
```

```python
nas2_single.index
```

```python
nas2_packed.index
```

Note that in the packed case, dtype='datetime64[ns]' --> This is different than what is in the single value case (dtype='datetime64[ns, UTC]') <br>
This makes me think that there might be a 37s problem!


So let's look at the index times. The first and last time in the packed series must correspond to the first and last time in the single value series

```python
nas2_packed_index_time=Time(nas2_packed.index)
```

```python
# The difference in seconds between the indexes
# without the .sec it outputs it in julian days!
(nas2_packed_index_time[0]-nas2_single_index_time[0]).sec
```

```python
(nas2_packed_index_time[-1]-nas2_single_index_time[-1]).sec
```

They do not! Note that this is using the indexes, and whatever format they are.
If we have astropy convert them to a standard, then they're still off by ~7 seconds! 

```python
nas2_packed_index_time.utc[0]
```

```python
nas2_single_index_time.utc[0]
```

```python
(nas2_packed_index_time.utc[0] - nas2_single_index_time.utc[0]).sec
```

```python
nas2_packed_index_time.utc
```

I can't figure out where this comes from! and this is why I *think* it's a bug...  
If it were attributable to the FPGA, it would be equal to ~1-3 seconds, as shown in the histogram plot above.
This offset greatly exceeds that value!


Just to hammer home my point, if we plot the data, the offset is apparent. <br>
This is what we are seeing all over the place and it's really hard to know what's going on!

```python
# Do some calculations so the y-axes aren't crazy
yr_cen=np.median(nas2_packed['nasmyth2CalculatedAngle'])
dy=1.1*(np.max(nas2_packed['nasmyth2CalculatedAngle']) - np.min(nas2_packed['nasmyth2CalculatedAngle']))/2

p = figure(x_axis_type='datetime', y_range=(yr_cen-dy, yr_cen+dy), plot_width=800, plot_height=400)
p.yaxis.axis_label = "Nasmyth Position (degrees)"
p.xaxis.axis_label = "Time"
p.line(x=nas2_packed_index_time.utc.value, y=nas2_packed['nasmyth2CalculatedAngle'], color='black', line_width=2, legend_label='Packed Data')
p.cross(x=nas2_single_index_time.utc.value, y=nas2_single['nasmyth2CalculatedAngle0'], color='red', line_width=2, line_dash='dashed', legend_label='single data')#

p.legend.location = 'bottom_left'
p.legend.click_policy = 'hide'
show(p)
```

For what it's worth, the CSC faulted due to this at 2022-02-16 05:15:27 [UTC]. 
So this indicates that the packed data is the one with the timestamp issues... which probably isn't surprising.

```python

```
