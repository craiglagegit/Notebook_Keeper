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

# Notebook for comparing eotest gain with DM gain.

Initially written 22 Nov 2021 by Craig Lage.

```python
import sys, os, glob, time
import pickle as pkl
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.io.fits as pf
from lsst.daf.butler import Butler
```

```python
butler = Butler("/repo/main", collections=["LSSTCam/raw/all","LSSTCam/calib",\
                                                    "LSSTCam/calib/u/cslage/13144", "u/cslage/bps_13144C"])
exposure=3021120600576
```

```python
# Get the eotest results
filename = "/project/cslage/BOT_LSSTCam/eotest/eotest_gain_13144_15dec21.pkl"
file = open(filename, 'rb')
#fe55_results = pkl.load(file)
ptc_results = pkl.load(file)
file.close()
print(ptc_results.keys())

rafts = [       'R01', 'R02', 'R03', \
         'R10', 'R11', 'R12', 'R13', 'R14', \
         'R20', 'R21', 'R22', 'R23', 'R24', \
         'R30', 'R31', 'R32', 'R33', 'R34', \
                'R41', 'R42', 'R43']
sensors = ['S00', 'S01', 'S02', 'S10', 'S11', 'S12', 'S20', 'S21', 'S22']
print(rafts)
print(ptc_results['ptc_gain']['R01']['S00'])
```

```python
# Get the "badAmps" results.  These are amps where the gains differed by more than 5%
filename = '/project/cslage/BOT_LSSTCam/eotest/badList_13144C.txt'
file = open(filename, 'r')
lines = file.readlines()
file.close()
badAmps = []
for line in lines:
    items = line.split()
    badAmps.append(items)
print(len(badAmps))
```

```python
def getDetector(raft, sensor):
    # Subroutine to find vendor and detector number given raft and sensor                                                                                                                                                           
    startingCol = [1,0,0,0,1] # First raft column in each row                                                                                                                                                                       
    rows = [0,3,8,13,18] # Starting raft sequence number of each row                                                                                                                                                                
    if raft in ['R11','R12','R13','R14','R21','R22','R23','R24','R30',\
                'R31','R32','R33','R34']:
        vendor = 'E2V'
    else:
        vendor = 'ITL'
    raftRow = int(list(raft)[1])
    raftCol = int(list(raft)[2])
    sensorRow = int(list(sensor)[1])
    sensorCol = int(list(sensor)[2])
    detectorNum = (rows[raftRow] + (raftCol - startingCol[raftRow])) * 9
    detectorNum += 3 * sensorRow + sensorCol
    return vendor, detectorNum

# This dictionary captures the amp naming correspondence
slacAmps = {'C10':'AMP01','C11':'AMP02','C12':'AMP03','C13':'AMP04',\
           'C14':'AMP05','C15':'AMP06','C16':'AMP07','C17':'AMP08',\
           'C07':'AMP09','C06':'AMP10','C05':'AMP11','C04':'AMP12',\
           'C03':'AMP13','C02':'AMP14','C01':'AMP15','C00':'AMP16'}
```

```python
def CalcDifferences(rawMeans, rawVars, slacMeans, slacVars):
    # weed out points where they don't match and return the percent differences
    meanMatches = []
    meanDiff = []
    varDiff = []
    for i, mean in enumerate(rawMeans):
        for j, sMean in enumerate(slacMeans):
            mDiff = abs(mean - sMean) / sMean * 100.0
            if mDiff < 5.0:
                meanMatches.append(mean)
                meanDiff.append(mDiff)
                vDiff = abs(rawVars[i] - slacVars[j]) / slacVars[j] * 100.0
                varDiff.append(vDiff)
                break
    return meanMatches, meanDiff, varDiff
    
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

```python jupyter={"outputs_hidden": true} tags=[]
pdf = PdfPages("/repo/main/u/cslage/bps_13144C/plots/Gain_Differences_GT_5Pct_13144C_21Dec21.pdf")
for ii, [RAFT,SENSOR,amp] in enumerate(badAmps):
    slacAmp = slacAmps[amp]
    slacNum = int(slacAmp.strip('AMP')) - 1
    VENDOR, DETECTOR  = getDetector(RAFT, SENSOR)
    fig = plt.figure(figsize=(16,16))
    ax1 = plt.axes([0.1,0.4,0.8,0.5])
    ptcDataset = butler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')
    dmGain = ptcDataset.gain[amp]
    if ptcDataset.ptcFitType == 'EXPAPPROXIMATION':
        dmA00 = ptcDataset.ptcFitPars[amp][0]
    if ptcDataset.ptcFitType == 'FULLCOVARIANCE':
        dmA00 = ptcDataset.aMatrix[amp][0][0]
    dmNoise = ptcDataset.noise[amp]
    rawMeans = ptcDataset.rawMeans[amp]
    rawVars = ptcDataset.rawVars[amp]
    dmMeans = np.array(ptcDataset.finalMeans[amp])
    dmMeans = dmMeans[~np.isnan(dmMeans)]
    if len(dmMeans > 0):
        maxDM = dmMeans.max()
    else:
        maxDM = 0.0
    eo_PtcTurnoff = ptc_results['ptc_turnoff'][RAFT][SENSOR][slacNum]
    filename = "/project/cslage/BOT_LSSTCam/eotest/%s_%s_13144_ptc.fits"%(RAFT,SENSOR)
    hdu = pf.open(filename)
    slacData = hdu[1].data
    slacMeans = slacData['%s_MEAN'%slacAmp]
    slacVars = slacData['%s_VAR'%slacAmp]
    slacGain = ptc_results['ptc_gain'][RAFT][SENSOR][slacNum]
    slacA00 = - ptc_results['ptc_a00'][RAFT][SENSOR][slacNum]
    slacNoise = ptc_results['ptc_noise'][RAFT][SENSOR][slacNum]
    if [RAFT,SENSOR,amp] in [['R30', 'S00', 'C10'], ['R03', 'S11', 'C00']]:
        ax1.text(10000,80000,"Gain Difference = DEAD", fontsize=18)
    elif dmGain > 1.0E-9:
        gainDiff = abs((dmGain - slacGain) / dmGain) * 100.0
        ax1.text(10000,80000,"Gain Difference = %.2f %%"%gainDiff, fontsize=18)
    else:
        gainDiff = 100.0
        ax1.text(10000,80000,"Gain Difference = %.2f %%"%gainDiff, fontsize=18)
    xplot = np.linspace(0,120000,200)
    
    ax1.set_title("PTC, %s_%s_%s_Det_%d_%s"%(RAFT,SENSOR, amp, DETECTOR, VENDOR), fontsize = 18)
    ax1.scatter(rawMeans, rawVars, marker = 'x', s=200, color = 'red', label = 'DM')
    ax1.scatter(slacMeans, slacVars, marker = '+', s=200, color = 'green', label = 'SLAC')
    ax1.text(10000,120000,"DM Max ADU = %.1f"%maxDM, fontsize=18)
    ax1.text(10000,110000,"EO PTC Turnoff = %.1f"%eo_PtcTurnoff, fontsize=18)
    ax1.text(10000,100000,"DM Gain = %.6f"%dmGain, fontsize=18)
    ax1.text(10000,90000,"EO Gain = %.6f"%slacGain, fontsize=18)
    #ax1.text(10000,100000,"Problem = %s"%fault, fontsize=18)
    yplot = ExpApprox(xplot,dmGain,dmA00,dmNoise)
    ax1.plot(xplot, yplot, ls = '--', color = 'red')
    yplot = ExpApprox(xplot,slacGain,slacA00,slacNoise)
    ax1.plot(xplot, yplot, ls = '--', color = 'green')
    ax1.set_xlabel("Mean(ADU)", fontsize=18)
    ax1.set_ylabel("Variance(ADU)", fontsize=18)
    #plt.xscale('log')
    #plt.yscale('log')
    ax1.set_xlim(0,150000)
    ax1.set_ylim(0,250000)
    ax1.legend(fontsize = 18)
    # Now the percent difference plots
    meanMatches, meanDiff, varDiff = CalcDifferences(rawMeans, rawVars, slacMeans, slacVars)
    ax2 = plt.axes([0.1,0.1,0.35,0.25])
    ax2.set_title("Mean Percent Difference", fontsize=12)
    ax2.scatter(meanMatches, meanDiff, label = "%s_%s"%(RAFT,SENSOR))
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(0.001, 1000)
    ax2.set_xlabel("Mean(ADU)", fontsize=12)
    ax2.set_ylabel("Mean Diff (percent)", fontsize=12)
    #ax2.legend()
    ax3 = plt.axes([0.55,0.1,0.35,0.25])
    ax3.set_title("Var Percent Difference", fontsize=12)
    ax3.scatter(meanMatches, varDiff, label = "%s_%s"%(RAFT,SENSOR))
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylim(0.001, 1000)
    ax3.set_xlabel("Mean(ADU)", fontsize=12)
    ax3.set_ylabel("Var Diff (percent)", fontsize=12)
    #ax3.legend()
    pdf.savefig(fig)  # saves the current figure into a pdf page
pdf.close()

plt.close()
```

```python
pdf = PdfPages("/repo/main/u/cslage/bps_13144C/plots/Amp_Det72_C17_13144C_21Dec21.pdf")
for ii in range(1):
    [RAFT,SENSOR,amp] = ['R20', 'S00', 'C17']
    slacAmp = slacAmps[amp]
    slacNum = int(slacAmp.strip('AMP')) - 1
    VENDOR, DETECTOR  = getDetector(RAFT, SENSOR)
    fig = plt.figure(figsize=(16,16))
    ax1 = plt.axes([0.1,0.4,0.8,0.5])
    try:
        ptcDataset = butler.get('ptc', detector=DETECTOR, exposure=exposure, instrument='LSSTCam')
        dmGain = ptcDataset.gain[amp]
        if ptcDataset.ptcFitType == 'EXPAPPROXIMATION':
            dmA00 = ptcDataset.ptcFitPars[amp][0]
        if ptcDataset.ptcFitType == 'FULLCOVARIANCE':
            dmA00 = ptcDataset.aMatrix[amp][0][0]
        dmNoise = ptcDataset.noise[amp]
        rawMeans = ptcDataset.rawMeans[amp]
        rawVars = ptcDataset.rawVars[amp]
        dmMeans = np.array(ptcDataset.finalMeans[amp])
        dmMeans = dmMeans[~np.isnan(dmMeans)]
        if len(dmMeans > 0):
            maxDM = dmMeans.max()
        else:
            maxDM = 0.0
    except:
        pass
    eo_PtcTurnoff = ptc_results['ptc_turnoff'][RAFT][SENSOR][slacNum]
    filename = "/project/cslage/BOT_LSSTCam/eotest/%s_%s_13144_ptc.fits"%(RAFT,SENSOR)
    hdu = pf.open(filename)
    slacData = hdu[1].data
    slacMeans = slacData['%s_MEAN'%slacAmp]
    slacVars = slacData['%s_VAR'%slacAmp]
    slacGain = ptc_results['ptc_gain'][RAFT][SENSOR][slacNum]
    slacA00 = - ptc_results['ptc_a00'][RAFT][SENSOR][slacNum]
    slacNoise = ptc_results['ptc_noise'][RAFT][SENSOR][slacNum]
    """
    if [RAFT,SENSOR,amp] in [['R30', 'S00', 'C10'], ['R03', 'S11', 'C00']]:
        ax1.text(10000,80000,"Gain Difference = DEAD", fontsize=18)
    elif dmGain > 1.0E-9:
        gainDiff = abs((dmGain - slacGain) / dmGain) * 100.0
        ax1.text(10000,80000,"Gain Difference = %.2f %%"%gainDiff, fontsize=18)
    else:
        gainDiff = 100.0
        ax1.text(10000,80000,"Gain Difference = %.2f %%"%gainDiff, fontsize=18)
    """
    xplot = np.linspace(0,120000,200)
    
    ax1.set_title("PTC, %s_%s_%s_Det_%d_%s"%(RAFT,SENSOR, amp, DETECTOR, VENDOR), fontsize = 18)
    #ax1.scatter(rawMeans, rawVars, marker = 'x', s=200, color = 'red', label = 'DM')
    ax1.scatter(slacMeans, slacVars, marker = '+', s=200, color = 'green', label = 'SLAC')
    #ax1.text(10000,120000,"DM Max ADU = %.1f"%maxDM, fontsize=18)
    ax1.text(10000,110000,"EO PTC Turnoff = %.1f"%eo_PtcTurnoff, fontsize=18)
    #ax1.text(10000,100000,"DM Gain = %.6f"%dmGain, fontsize=18)
    ax1.text(10000,90000,"EO Gain = %.6f"%slacGain, fontsize=18)
    #ax1.text(10000,100000,"Problem = %s"%fault, fontsize=18)
    #yplot = ExpApprox(xplot,dmGain,dmA00,dmNoise)
    #ax1.plot(xplot, yplot, ls = '--', color = 'red')
    yplot = ExpApprox(xplot,slacGain,slacA00,slacNoise)
    ax1.plot(xplot, yplot, ls = '--', color = 'green')
    ax1.set_xlabel("Mean(ADU)", fontsize=18)
    ax1.set_ylabel("Variance(ADU)", fontsize=18)
    #plt.xscale('log')
    #plt.yscale('log')
    ax1.set_xlim(0,150000)
    ax1.set_ylim(0,250000)
    ax1.legend(fontsize = 18)
    """
    # Now the percent difference plots
    meanMatches, meanDiff, varDiff = CalcDifferences(rawMeans, rawVars, slacMeans, slacVars)
    ax2 = plt.axes([0.1,0.1,0.35,0.25])
    ax2.set_title("Mean Percent Difference", fontsize=12)
    ax2.scatter(meanMatches, meanDiff, label = "%s_%s"%(RAFT,SENSOR))
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_ylim(0.001, 1000)
    ax2.set_xlabel("Mean(ADU)", fontsize=12)
    ax2.set_ylabel("Mean Diff (percent)", fontsize=12)
    #ax2.legend()
    ax3 = plt.axes([0.55,0.1,0.35,0.25])
    ax3.set_title("Var Percent Difference", fontsize=12)
    ax3.scatter(meanMatches, varDiff, label = "%s_%s"%(RAFT,SENSOR))
    ax3.set_xscale('log')
    ax3.set_yscale('log')
    ax3.set_ylim(0.001, 1000)
    ax3.set_xlabel("Mean(ADU)", fontsize=12)
    ax3.set_ylabel("Var Diff (percent)", fontsize=12)
    #ax3.legend()
    """
    pdf.savefig(fig)  # saves the current figure into a pdf page
pdf.close()

plt.close()
```

```python

```
