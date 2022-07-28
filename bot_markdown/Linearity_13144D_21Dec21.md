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

# Notebook for investigating linearity corrections

Initially written 20 Dec 2021 by Craig Lage\
copying from Chris Waters.

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
from lsst.daf.butler import Butler
import lsst.afw.math as afwMath
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib",\
                                                    "LSSTCam/calib/u/cslage/13144", "u/cslage/bps_13144D"])
expId=3021120600576
```

```python
camera = butler.get('camera', instrument='LSSTCam')
offset = 0
names = ["E2V", "ITL"]
for i, det in enumerate([100, 181]):
    lin = butler.get('linearity', detector=det, exposure=expId, instrument='LSSTCam')
    ptc = butler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
    plt.subplot(1,2,i+1)
    plt.title(names[i]+f" - {det}")
    for it, amp in enumerate(camera[0].getAmplifiers()):
        centers, values = np.split(lin.linearityCoeffs[amp.getName()], 2)
        plt.scatter(centers, values + it * offset, marker='+')
plt.subplots_adjust(wspace=0.5)
#plt.savefig("/repo/main/u/cslage/bps_13144D/plots/Linearity_E2V_ITL_21Dec21.pdf")
```

```python
#E2V
det = 55
lin = butler.get('linearity', detector=det, exposure=expId, instrument='LSSTCam')
ptc = butler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
camera = butler.get('camera', instrument='LSSTCam')
```

```python tags=[]
for it, amp in enumerate(camera[0].getAmplifiers()):
    print(lin.linearityType[amp.getName()])
    print(lin.linearityCoeffs[amp.getName()])
    break
```

```python
offset = 0
for it, amp in enumerate(camera[0].getAmplifiers()):
    centers, values = np.split(lin.linearityCoeffs[amp.getName()], 2)
    plt.scatter(centers, values + it * offset, marker='+', label="LinA")

```

```python
def ExpApprox(mu, g, a00, n):
    if (g < 1.0E-6) or (abs(a00) < 1.0E-9):
        return np.zeros([len(mu)])
    else:
        expFactor = 2.0 * a00 * mu * g
        if max(expFactor) > 100.0:
            return np.zeros([len(mu)])
        else:
            preFactor = 1.0 / (2.0 * g * g * a00)
            noiseTerm = n / (g * g)
            return preFactor * (np.exp(expFactor) - 1.0) + noiseTerm
```

```python
for amp in ptc.ampNames:
    centers, values = np.split(lin.linearityCoeffs[amp], 2)
    interp = afwMath.makeInterpolate(centers.tolist(), values.tolist(),
                                     afwMath.stringToInterpStyle("AKIMA_SPLINE"))
    delta = interp.interpolate(np.array(ptc.rawMeans[amp]))
    linearized = np.array(ptc.rawMeans[amp]) - np.array(delta) # ?? Adjust X-axis??
    #linearized2 = np.array(ptc.rawVars[amp]) - np.array(delta)
    gain = ptc.gain[amp]
    a00 = ptc.ptcFitPars[amp][0]
    noise = ptc.noise[amp]
    
    plt.subplots_adjust(wspace = 0.5)
    
    plt.subplot(1,3,1)
    yplot = ExpApprox(np.array(ptc.rawMeans[amp]), gain, a00, noise)
    plt.scatter(linearized, ptc.rawVars[amp], marker='o', label="Linearized")
    #plt.scatter(ptc.rawMeans[amp], linearized2, marker='o', label="Linearized")
    plt.scatter(ptc.rawMeans[amp], ptc.rawVars[amp], marker='+', label="Raw")
    plt.plot(ptc.rawMeans[amp], yplot, ls = '--', color = 'red', label = 'ExpApprox')
    plt.xlim(50000,80000)
    plt.ylim(30000, 40000)
    plt.subplot(1,3,2)
    plt.scatter(ptc.rawMeans[amp], yplot - ptc.rawVars[amp], marker='+', label="Raw")
    plt.xlim(0,100000)
    plt.ylim(-1000,1000)
    plt.legend()
    plt.subplot(1,3,3)
    plt.scatter(ptc.rawMeans[amp], delta, marker='+', label="Raw")
    plt.xlim(0,100000)
    plt.ylim(-1000,1000)
    plt.legend()

    break


```

```python
#ITL
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib",\
                                                    "LSSTCam/calib/u/cslage/13144", "u/cslage/bps_13144B"])
expId=3021120600576
det = 74
lin = butler.get('linearity', detector=det, exposure=expId, instrument='LSSTCam')
ptc = butler.get('ptc', detector=det, exposure=expId, instrument='LSSTCam')
camera = butler.get('camera', instrument='LSSTCam')
```

```python tags=[]
for it, amp in enumerate(camera[0].getAmplifiers()):
    print(lin.linearityType[amp.getName()])
    print(lin.linearityCoeffs[amp.getName()])
    break
```

```python
offset = 0
for it, amp in enumerate(camera[0].getAmplifiers()):
    centers, values = np.split(lin.linearityCoeffs[amp.getName()], 2)
    #print(centers)
    #print(values)
    #break
    plt.scatter(centers, values + it * offset, marker='+', label="Lin")

```

```python
def ExpApprox(mu, g, a00, n):
    if (g < 1.0E-6) or (abs(a00) < 1.0E-9):
        return np.zeros([len(mu)])
    else:
        expFactor = 2.0 * a00 * mu * g
        if max(expFactor) > 100.0:
            return np.zeros([len(mu)])
        else:
            preFactor = 1.0 / (2.0 * g * g * a00)
            noiseTerm = n / (g * g)
            return preFactor * (np.exp(expFactor) - 1.0) + noiseTerm
```

```python
for amp in ptc.ampNames:
    centers, values = np.split(lin.linearityCoeffs[amp], 2)
    interp = afwMath.makeInterpolate(centers.tolist(), values.tolist(),
                                     afwMath.stringToInterpStyle("AKIMA_SPLINE"))
    delta = interp.interpolate(np.array(ptc.rawMeans[amp]))
    linearized = np.array(ptc.rawMeans[amp]) - np.array(delta) # ?? Adjust X-axis??
    #linearized2 = np.array(ptc.rawVars[amp]) - np.array(delta)
    gain = ptc.gain[amp]
    a00 = ptc.ptcFitPars[amp][0]
    noise = ptc.noise[amp]
    
    plt.subplots_adjust(wspace = 0.5)
    
    plt.subplot(1,3,1)
    yplot = ExpApprox(np.array(ptc.rawMeans[amp]), gain, a00, noise)
    plt.scatter(linearized, ptc.rawVars[amp], marker='o', label="Linearized")
    #plt.scatter(ptc.rawMeans[amp], linearized2, marker='o', label="Linearized")
    plt.scatter(ptc.rawMeans[amp], ptc.rawVars[amp], marker='+', label="Raw")
    plt.plot(ptc.rawMeans[amp], yplot, ls = '--', color = 'red', label = 'ExpApprox')
    plt.xlim(50000,80000)
    plt.ylim(30000, 40000)
    plt.subplot(1,3,2)
    plt.scatter(ptc.rawMeans[amp], yplot - ptc.rawVars[amp], marker='+', label="Raw")
    plt.xlim(0,100000)
    plt.ylim(-1000,1000)
    plt.legend()
    plt.subplot(1,3,3)
    plt.scatter(ptc.rawMeans[amp], delta, marker='+', label="Raw")
    plt.xlim(0,100000)
    plt.ylim(-1000,1000)
    plt.legend()

    break


```

```python

```
