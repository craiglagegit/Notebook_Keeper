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

# Notebook for eotest gain with DM gain.

Initially written 21 May 2020 by Craig Lage.

```python
! eups list -s | grep lsst_distrib
! eups list -s | grep obs_lsst
```

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.persistence import Butler

```

```python
DATA_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_12597'
[RAFT, SENSOR, DETECTOR] = ['R34', 'S11', 157]
butler = Butler(DATA_DIR)

```

```python
dataRef = butler.dataRef('raw', dataId={'detector': 157, 'expId': 3020100600337, 'dayObs': '2020-10-06'})
```

```python
exp = dataRef.get('raw')
```

```python
ccd = exp.getDetector()

for amp in ccd:
    name = amp.getName()
    gain = amp.getGain()
    print(name, gain)
```

```python
dataRef.get('bias')
```

```python
dataRef.get('dark')
```

```python
dataRef.get('flat')
```

```python
dataRef.get('defects')
```

```python
dataRef.get('dark')
```

```python

```

```python
exp1 = butler.get('postISRCCD', dataId={'detector': 94,  'visit': 3019101200270})
```

```python
md1 = exp1.getMetadata()
```

```python
md1['FILTER2']
```

```python
dataRef = butler.dataRef('postISRCCD',dataId={'detector': 94,  'visit': 3019101200270})
```

```python
exp1 = dataRef.get('postISRCCD')
```

```python
md1 = exp1.getMetadata()
```

```python
md1['MONDIODE']
```

```python
test = butler.dataRef('raw', dataId={'detector': 94,  'visit': 3019101200270})
test.get('raw')
test.get('bias')
test.get('dark')
test.get('flat')
test.get('defect')
```

<!-- #raw -->
from collections import OrderedDict
test = butler.dataRef('raw', dataId=OrderedDict([('visit', '3019101200250'), ('detector', 0)]))
<!-- #endraw -->

```python

```

```python
test = test.get('raw', dataId={'detector': 0,  'visit': 3019101200250})
```

```python
test = butler.dataRef('raw', dataId={'detector': DETECTOR,  'expId': '3019101200250'})
```

```python
test = butler.dataRef('bias', dataId={'detector': DETECTOR,  'expId': 3019101200480})
```

```python
test = butler.dataRef('flat', dataId={'detector': DETECTOR,  'expId': 3019101200480})
```

```python
test = butler.dataRef('dark', dataId={'detector': DETECTOR,  'expId': 3019101200480})
```

```python
test = butler.dataRef('defects', dataId={'detector': DETECTOR,  'expId': 3019101200480})
```

```python
test = butler.dataRef('postISRCCD', dataId={'detector': DETECTOR,  'expId': 3019101200999})
```

```python
test = butler.get('defects', dataId={'detector': DETECTOR,  'expId': 3019101200480})
print(test.getMetadata())
```

```python
test = butler.get('flat', dataId={'detector': DETECTOR,  'expId': 3019101200480})
print(test.getMetadata())
```

```python
test = butler.get('dark', dataId={'detector': DETECTOR,  'expId': 3019101200480})
print(test.getMetadata())
```

```python
test = butler.get('raw', dataId={'detector': DETECTOR,  'expId': 3019101200480})
print(dir(test))
```

```python

```

```python
REPO_DIR_1 = '/project/cslage/BOT_lspdev/E2V_6790D_Gain_Edge4_10_R22S11'
REPO_DIR_2 = '/project/shared/BOT/rerun/cslage/PTC_6790D_NewAll'
#REPO_DIR_1 = '/project/shared/BOT/rerun/cslage/Test_ISR1'
#REPO_DIR_2 = '/project/shared/BOT/rerun/cslage/ISR_Test10'
butler1 = Butler(REPO_DIR_1)
butler2 = Butler(REPO_DIR_2)
visit = 3019101300324
exp1 = butler1.get('postISRCCD', raftName=RAFT, detectorName=SENSOR, visit=visit)
exp2 = butler2.get('postISRCCD', raftName=RAFT, detectorName=SENSOR, visit=visit)
#exp1 = butler1.get('defect', raftName=RAFT, detectorName=SENSOR, visit=visit)
#exp2 = butler2.get('defect', raftName=RAFT, detectorName=SENSOR, visit=visit)
ccd1 = exp1.getDetector()
for ampObject in ccd1:
    amp = ampObject.getName()
    if amp == 'C00':
        break

img1 = exp1.maskedImage[ampObject.getBBox()]
img2 = exp2.maskedImage[ampObject.getBBox()]
names = [REPO_DIR_1, REPO_DIR_2]
arrs = [img1.getArrays()[0][100:120,100:120],img2.getArrays()[0][100:120,100:120]] 

```

```python
# Now let's look at them
from matplotlib.colors import LogNorm
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
fig=plt.figure(figsize=(16,10))
plt.suptitle("%s%s - visit:%d - Amp %s [100:120,100:120], "%(RAFT,SENSOR,visit,amp), fontsize = 18)
for i, arr in enumerate(arrs):
    plt.subplot(1,2,i+1)
    plt.title("%s \n Mean = %.3f, Std = %.3f"%(names[i],arr.mean(),arr.std()))
    img1 = plt.imshow(arr)#,norm=LogNorm(vmin=35000,vmax=38000))#, vmin=50000, vmax=75000)
    colorbar(img1)
plt.tight_layout(h_pad=1)
#plt.savefig(REPO_DIR_2+"/plots/PostISR_Images_30Jul0.pdf")

```

<!-- #raw -->
ASH_REPO = '/home/adriansh/lsst_devel/analysis/satellite/20200616-starlink1/'
butler = Butler(ASH_REPO)


<!-- #endraw -->

```python

```
