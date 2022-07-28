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

# Notebook for plotting BF Kernels generated with Gen 3.

Initially written 22 Jun 2021 by Craig Lage.

```python
import os, sys, time, datetime, glob, subprocess
import matplotlib.pyplot as plt
import numpy as np
import astropy.io.fits as pf
from scipy import stats
from lsst.daf.butler import Butler
```

```python
butler = Butler('/repo/main', collections=
                ['LSSTComCam/raw/all', 'LSSTComCam/calib', 'LSSTComCam/calib/u/cslage/20210402A'])
```

```python
expId = 2021040200025 # I think any exposure within the set of flat pairs will work.
DETECTOR = 4
dataId={'instrument':'LSSTComCam', 'detector':4, 'exposure':expId}
```

```python
ptc_dataset = butler.get('ptc', dataId=dataId)
bf_kernel = butler.get('bfk', dataId=dataId)
gains = bf_kernel.gain
means = bf_kernel.means # Mean flux of flat pairs in electrons
rawMeans = ptc_dataset.rawMeans # Mean flux of flat pairs in electrons
rawVars = ptc_dataset.rawVars # Mean flux of flat pairs in electrons
rawXcorrs = bf_kernel.rawXcorrs # Raw extracted covariances in ADU^2. [0,0] is the variance}
meanXcorrs = bf_kernel.meanXcorrs # Extracted covariances used to extract kernel. These are per e-.
kernels = bf_kernel.ampKernels # ampwise kernel
ptcResults = ptc_dataset.ptcFitPars
```

```python
print(ptc_dataset.covariances['C01'][4][0][0], bf_kernel.rawXcorrs['C01'][4][0,0], ptc_dataset.rawVars['C01'][4], ptc_dataset.finalVars['C01'][2])
```

```python jupyter={"outputs_hidden": true}
for ii in range(8):
    for jj in range(8):
        print(ptc_dataset.covariances['C01'][4][ii][jj], bf_kernel.rawXcorrs['C01'][4][ii,jj])
```

```python
# Now plot the correlations and the kernel. 
# So much variation in the (0,0) covariance!
for amp in means.keys():
    fig = plt.figure(figsize=(16,5))
    plt.subplots_adjust(hspace=0.5)
    plt.suptitle("COVARIANCES(*1E7)       Amp %s      KERNEL(*1E7)"%amp, fontsize=24)
    plt.subplot(1,4,1)
    plt.imshow(np.log10(abs(np.array(meanXcorrs[amp]))))
    plt.subplot(1,4,2)
    plt.plot([0,16],[0,0], ls='--', color='black')
    plt.plot(-meanXcorrs[amp][:,8]*1E7, color='blue', drawstyle='steps-mid')
    plt.plot(-meanXcorrs[amp][8,:]*1E7, linestyle='--', color='red', drawstyle='steps-mid')
    plt.ylim(-20,10)
    plt.subplot(1,4,3)
    plt.imshow(kernels[amp])
    plt.subplot(1,4,4)  
    plt.plot([0,16],[0,0], ls='--', color='black')
    plt.plot(kernels[amp][:,8]*1E7, color='blue', drawstyle='steps-mid')
    plt.plot(kernels[amp][8,:]*1E7, linestyle='--', color='red', drawstyle='steps-mid')
    plt.ylim(-10,2)
```

```python

```

```python
ptcResults['C01']
```

```python
gains['C01']
```

```python
# Next, reproduce the Photon Transfer Curves

plt.figure(figsize=(16,8))
plt.subplots_adjust(hspace=0.3,wspace=0.02)
plt.suptitle("ComCam",fontsize = 24)
plotcounter = 0

for amp in means.keys():
    plotcounter += 1
    plt.subplot(4,4,plotcounter)
    plt.title("PTC %s"%amp,fontsize=12)
    data = []
    for m, flux in enumerate(rawMeans[amp]):
        var = rawXcorrs[amp][m][0,0]# / 2.0 # Divide by two because variance of a flat difference is 2X.
        data.append([flux, var])
    data.sort()
    data = np.array(data)
    quad_fit = np.polyfit(data[:,0], data[:,1], 2)
    plt.scatter(data[:,0], data[:,1],marker='x',color='green')
    gain = ptcResults[amp][1]
    rms_noise = np.sign(ptcResults[amp][2]) * np.sqrt(abs(ptcResults[amp][2])) * gain
    xplot=np.linspace(0.0, 150000.0, 100)
    yplot_quad = ptcResults[amp][0]*xplot*xplot + (1.0/ptcResults[amp][1])*xplot + ptcResults[amp][2]
    yplot_linear = (1.0/ptcResults[amp][1])*xplot + ptcResults[amp][2]
    plt.plot(xplot,yplot_linear,color='green', ls='--',label='Linear Fit;Gain=%.3f;Noise=%.1f e-'%(gain,rms_noise))
    plt.plot(xplot,yplot_quad,color='red', label='Quadratic fit')

    plt.ylim(0,100000)
    plt.xlim(0,150000)
    #plt.xticks([0,10000,20000])
    #plt.yticks([0,2000,4000])
    plt.tick_params(left=False,  bottom=False, labelleft=False,  labelbottom=False)
    if plotcounter in [1,5,9,13]:
        plt.ylabel("Variance(ADU^2)",fontsize=10)
        plt.tick_params(left=True, labelleft=True)
    if plotcounter in [13,14,15,16]:
        plt.xlabel("Flux(ADU)", fontsize=10)
        plt.xticks([0,25000,50000,75000,100000])
        plt.tick_params(bottom=True, labelbottom=True)
        
    plt.legend(loc = 'lower right', fontsize = 6)
#plt.savefig(OUTPUT_DIR+'plots/PTC_Quad_All_%s_%s.pdf'%(date,DETECTOR))
```

```python
# Next plot the covariance vs flux 
# Had to add a maxFlux here to weed out points >150K.
# Is that the best method?
#date = OUTPUT_DIR.split('/')[4]
NumPairs = 1 # Number of pairs at each flux
fluxLevel = 80000.0
maxMean = 120000.0
PlotDelta = 5 # Number of pixels to look at
#pdf = PdfPages(OUTPUT_DIR+"plots/Covariance_vs_Flux_%s_%s.pdf"%(date,DETECTOR))
for amp in ['C04']:#means.keys():
    gain = gains[amp]
    NumFluxes = int(len(rawMeans[amp]) / NumPairs)
    fig = plt.figure(figsize = (16,8))
    plt.suptitle("Covariance vs Flux - Amp %s"%amp, fontsize = 24)
    plt.subplots_adjust(wspace=0.3, hspace=0.6)
    plotcounter = 0
    for jj in range(PlotDelta-1, -1, -1):
        for ii in range(PlotDelta):
            cov_mean = []
            cov_std = []
            flux_mean = []
            plotcounter += 1
            plt.subplot(PlotDelta, PlotDelta, plotcounter)
            cov = []
            flux = []

            for n in range(NumFluxes):
                cov_per_flux = []
                flux_per_flux = []
                for m in range(NumPairs):
                    i = n * NumPairs + m
                    xcorr = rawXcorrs[amp][i][ii,jj] * gain**2
                    mean = rawMeans[amp][i] * gain
                    if mean > maxMean:
                        continue
                    if ii == 0 and jj == 0:
                        # This is right, but needs double-checking
                        xcorr = xcorr - (mean / gain / ptcResults[amp][1] + ptcResults[amp][2])*gain**2
                        cov.append(-xcorr)
                        cov_per_flux.append(-xcorr)
                    else:
                        cov.append(xcorr)
                        cov_per_flux.append(xcorr)
                    flux.append(mean)
                    flux_per_flux.append(mean)
                if len(cov_per_flux) == 0:
                    continue
                cov_per_flux = np.array(cov_per_flux)
                flux_per_flux = np.array(flux_per_flux)
                cov_mean.append(cov_per_flux.mean())
                cov_std.append(cov_per_flux.std())
                flux_mean.append(flux_per_flux.mean())
            cov = np.array(cov)
            flux = np.array(flux)
            plt.scatter(flux, cov, color='blue', marker='.', s=20)
            cov_mean = np.array(cov_mean)
            cov_std = np.array(cov_std)
            flux_mean = np.array(flux_mean)
            cov_mean = cov_mean#*flux_mean*flux_mean
            cov_std = cov_std#*flux_mean*flux_mean

            plt.errorbar(flux_mean,cov_mean,yerr=cov_std,color='green',marker='x',ls='None')  
            
            if ii == 0 and jj == 0:
                corrValue = meanXcorrs[amp][jj+8,ii+8]*fluxLevel*fluxLevel
            else:
                corrValue = -meanXcorrs[amp][jj+8,ii+8]*fluxLevel*fluxLevel
            plt.scatter(fluxLevel, corrValue, color = 'red', marker='o')
            plt.title("Pixel: (%d, %d)"%(ii, jj), fontsize = 12)
            coefs = np.polyfit(flux*flux, cov, 1)
            xplot = np.linspace(0,150000, 100)
            yplot = max(0, coefs[0])*xplot*xplot
            plt.plot(xplot,yplot, color = 'red', lw = 2)
            if jj == 0:
                plt.xlabel("Central Pixel Charge(e-)", fontsize = 12)
            if ii == 0:
                plt.ylabel("Correlation", fontsize = 12)
            plt.xlim(0,150000)
            plt.xticks([0,100000],fontsize = 12)
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
                plt.ylim(-500,2000)
            else:
                plt.yticks([-200,0,200],fontsize = 12)
                plt.ylim(-200,500)

    #pdf.savefig(fig)  # saves the current figure into a pdf page
    #plt.close()
#pdf.close()

```

```python
# Plot the correlations as a function of pixel 
# Current implementation isn't actually implementing the model, and isn't doing the
# forceZeroSum correctly.
#date = OUTPUT_DIR.split('/')[4]
#pdf = PdfPages(OUTPUT_DIR+"plots/Covariance_Matrix_%s_%s.pdf"%(date,DETECTOR))
NumPairs = 1
NumFluxes = int(len(means['C12']) / NumPairs)
for amp in ['C04']:#means.keys():
        #try:
        gain = gains[amp]
        NumFluxes = int(len(means[amp]) / NumPairs)
        posrs = []
        negrs = []
        poscs = []
        negcs = []
        fitrs = []
        fitcs = []
        poserrs = []
        negerrs = []
        for ii in range(1,16):
            for jj in range(1,16):
                r2 = (ii-8)**2 + (jj-8)**2
                value = meanXcorrs[amp][ii,jj]
                n = NumFluxes - 1
                cov = []
                for m in range(NumPairs):
                    i = n * NumPairs + m
                    xcorr = rawXcorrs[amp][i][abs(ii-8),abs(jj-8)] * gain**2
                    mean = rawMeans[amp][i] * gain
                    cov.append(xcorr)
                cov = np.array(cov)
                if ii == 8 and jj == 8:
                    negcs.append(abs(value))
                    negrs.append(0.85)
                    negerrs.append(abs(cov.std()))
                elif value < 0.0:
                    poscs.append(abs(value))
                    posrs.append(r2)
                    poserrs.append(abs(cov.std()))
                    if r2 > 1.1 and r2 < 20.0:
                        fitrs.append(np.log10(r2))
                        fitcs.append(np.log10(abs(value)))
                else:
                    negcs.append(abs(value))
                    negrs.append(r2)
                    negerrs.append(abs(cov.std()))

        slope, intercept, r_value, p_value, std_err = stats.linregress(fitrs,fitcs)

        fig = plt.figure(figsize=(16,8))
        plt.title("Covariance Matrix - %s"%amp, fontsize = 24)
        plt.errorbar(posrs, poscs, yerr=poserrs, color='blue', marker = 'o', markersize = 8, ls = '', label = 'Positive Meas')
        plt.errorbar(negrs, negcs, yerr=negerrs, color='red', marker = 'o', markersize = 8, ls = '', label = 'Negative Meas')

        plt.text(1.2, 5E-6, "C00: Meas = %.4g"%(-meanXcorrs[amp][8,8]), fontsize = 24)
        plt.text(1.2, 2.5E-6, "C01: Meas = %.4g"%(-meanXcorrs[amp][9,8]), fontsize = 24)
        plt.text(1.2, 1.25E-6, "C10: Meas = %.4g"%(-meanXcorrs[amp][8,9]), fontsize = 24)
        plt.text(1.2, 6.25E-7, "C11: Meas = %.4g"%(-meanXcorrs[amp][9,9]), fontsize = 24)
        xplot = np.linspace(1.0,200,100)
        yplot = 10.0**intercept * xplot**(slope)
        #plot(xplot,yplot,ls='--',color='green')
        plt.xscale('log')
        plt.yscale('log')
        plt.xlim(0.8,100)
        plt.ylim(1E-10,1E-5)
        plt.xticks([1.0, 10.0, 100.0], fontsize=24)
        plt.yticks([1E-10, 1E-9, 1E-8, 1E-7, 1E-6, 1E-5], fontsize=24)
        plt.xlabel("$i^2 + j^2$", fontsize = 24)
        plt.ylabel("Covariance or $\delta$ Area/Area", fontsize = 24)
        
        
        # Plot the model
        preFactor = np.sqrt(meanXcorrs[amp][9,8] * meanXcorrs[amp][8,9])
        slopeFactor = -1.35
        rplot = np.linspace(1,100, 100)
        yrplot = preFactor * rplot**(slopeFactor)
        plt.plot(rplot,yrplot, color = 'red', lw = 2, ls='--', label="Model")

        plt.legend(fontsize=24)
        
        #pdf.savefig(fig)  # saves the current figure into a pdf page
        #plt.close()
        #except:
        #print("Skipping amp %s"%amp)
        #continue
#pdf.close()
```

```python

```
