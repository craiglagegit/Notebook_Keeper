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

# Notebook for plotting extracted covariances on BOT data

Initially written 16 Oct 2019 by Craig Lage.\


```python
! eups list -s | grep lsst_distrib
! eups list -s obs_lsst 
! eups list -s cp_pipe
```

```python
import os, sys, time, datetime, glob, subprocess
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import astropy.io.fits as pf
from scipy import stats
from lsst.daf.persistence import Butler
```

```python
# To remember, an alternate way to get the ptc data:
#from lsst.ip.isr import PhotonTransferCurveDataset
#datasetFile = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_New_12606/rerun/plazas/DM-27185/calibrations/ptc/ptcDataset-det094.fits'
#datasetPtc = PhotonTransferCurveDataset.readFits(datasetFile)
#print(datasetPtc.covariances)
```

```python
run = '12606'
DATA_DIR = '/project/shared/BOT/'
RAFT = 'R12'
SENSOR = 'S02'
REPO_DIR = '/project/shared/BOT/rerun/cslage/PTC_LSSTCAM_FullCov_12606'
DETECTOR = 47#94
```

```python
butler = Butler(REPO_DIR)
ptcDataset = butler.get('photonTransferCurveDataset', dataId={'detector': DETECTOR})
gains = ptcDataset.gain
means = ptcDataset.rawMeans
xCorrs = ptcDataset.covariances
```

```python
#print(dir(ptcDataset))
print(ptcDataset.ptcFitType)
print(ptcDataset.inputExpIdPairs['C15'])
means = ptcDataset.finalMeans['C15']
print(means)
vars = ptcDataset.finalVars['C15']
print(vars)
plt.plot(means,vars, marker='x')
plt.xscale('log')
plt.yscale('log')
```

```python
# Next plot the covariance vs flux 
for PlotDelta in [5,8]: # Number of pixels to look at
    serialCov = np.zeros([len(means.keys()),PlotDelta])
    pdf = PdfPages(REPO_DIR+"/plots/Covariance_vs_Flux_%s_%s.pdf"%(PlotDelta,DETECTOR))
    for ampNum, amp in enumerate(means.keys()):
        gain = gains[amp]
        NumFluxes = int(len(means[amp]))
        fig = plt.figure(figsize = (16,8))
        plt.suptitle("Covariance vs Flux - Amp %s, %s-%s"%(amp,RAFT,SENSOR), fontsize = 24)
        plt.subplots_adjust(wspace=0.3, hspace=0.6)
        plotcounter = 0
        for jj in range(PlotDelta-1, -1, -1):
            for ii in range(PlotDelta):
                plotcounter += 1
                plt.subplot(PlotDelta, PlotDelta, plotcounter)
                cov = []
                flux = []

                for n in range(NumFluxes):
                    xcorr = xCorrs[amp][n][ii][jj]
                    mean = means[amp][n]
                    if ii == 0 and jj == 0:
                        # This isn't right yet
                        xcorr = xcorr - mean * gain
                        cov.append(-xcorr)
                    else:
                        cov.append(xcorr)
                    flux.append(mean)
                cov = np.array(cov)
                flux = np.array(flux)

                plt.scatter(flux, cov, color='blue', marker='x', s=50)
                coefs = np.polyfit(flux*flux, cov, 1)
                if jj == 1:
                    serialCov[ampNum, ii] = coefs[0]
                xplot = np.linspace(0,150000, 100)
                yplot = max(0, coefs[0])*xplot*xplot
                plt.plot(xplot,yplot, color = 'red', lw = 2)
                plt.title("Pixel: (%d, %d)"%(ii, jj), fontsize = 12)
                if jj == 0:
                    plt.xlabel("Central Pixel Charge(e-)", fontsize = 12)
                if ii == 0:
                    plt.ylabel("Correlation", fontsize = 12)
                plt.xlim(0,120000)
                plt.xticks([0,150000],fontsize = 12)

                if ii == 0 and jj == 0:
                    plt.yticks([0,10000],fontsize = 12)
                    plt.ylim(-1000,30000)
                elif ii == 0 and jj == 1:
                    plt.yticks([0,2000,4000],fontsize = 12)
                    plt.ylim(-500,4000)
                elif ii == 2 and jj == 0:
                    plt.yticks([0,500,1000],fontsize = 12)
                    plt.ylim(-100,1000)
                elif ii == 1 and jj < 2:
                    plt.yticks([0,1000,2000],fontsize = 12)
                    plt.ylim(-500,4000)
                else:
                    plt.yticks([-200,0,200],fontsize = 12)
                    plt.ylim(-200,1000)


        pdf.savefig(fig)  # saves the current figure into a pdf page
        plt.close()
    pdf.close()

```

```python
# Plot the serial covariances vs distance 
plt.figure(figsize=(16,8))
plt.subplots_adjust(hspace=0.3,wspace=0.02)
plt.suptitle("SLAC run 12606/12610, %s-%s"%(RAFT,SENSOR),fontsize = 24)
plotcounter = 0
xaxis = np.arange(PlotDelta)
for ampNum, amp in enumerate(means.keys()):
    plotcounter += 1
    plt.subplot(4,4,plotcounter)
    plt.title("Serial Covariances %s"%amp,fontsize=12)
    plt.scatter(xaxis[1:], serialCov[ampNum,1:]*1.0E7 ,marker='x',color='green')
    plt.plot([1.0,8.0],[0.0,0.0], color='black', ls = '--')
    plt.ylim(-1.0,5.0)
    #plt.xlim(0,120000)                                                                                                                
    plt.tick_params(left=False,  bottom=False, labelleft=False,  labelbottom=False)
    if plotcounter in [1,5,9,13]:
        plt.ylabel("Covariance * 1E7",fontsize=10)
        plt.tick_params(left=True, labelleft=True)
    if plotcounter in [13,14,15,16]:
        plt.xlabel("X pixel", fontsize=10)
        #plt.xticks([0,25000,50000,75000,100000])                                                                                      
        plt.tick_params(bottom=True, labelbottom=True)

plt.savefig(REPO_DIR+'/plots/Cov_vs_X_16Oct20_%s.pdf'%(DETECTOR))

```

```python
# Plot the serial covariances vs distance                                                                                              
plt.figure(figsize=(16,8))
plt.suptitle("SLAC run 12606/12610, %s-%s"%(RAFT,SENSOR),fontsize = 24)
xaxis = np.arange(PlotDelta)
for ampNum, amp in enumerate(means.keys()):
    plt.plot(xaxis[1:], np.log10(serialCov[ampNum,1:]), label=amp)
plt.ylim(-9.0,-6.0)
plt.ylabel("Log Covariance",fontsize=18)
plt.xlabel("X pixel", fontsize=18)
plt.legend()
plt.savefig(REPO_DIR+'/plots/Cov_vs_X_Log_16Oct20_%s.pdf'%(DETECTOR))

```

```python

```

```python

```
