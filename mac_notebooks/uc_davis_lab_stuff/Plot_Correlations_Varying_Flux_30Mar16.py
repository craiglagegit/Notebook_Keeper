
# coding: utf-8

# In[2]:

# %load bias_stability.py
#!/usr/bin/python

# Plots stability of a set of bias frames vs time

import matplotlib
matplotlib.use("PDF")
import pyfits as pf
from pylab import *
import sys, glob, time
from scipy import stats

thedir='/Users/cslage/Research/LSST/code/GUI/flats/20160107_varying_flux/'
get_ipython().magic(u'cd $thedir')

keys = ['SEGMENT10','SEGMENT11','SEGMENT12','SEGMENT13','SEGMENT14','SEGMENT15',
        'SEGMENT16','SEGMENT17','SEGMENT07','SEGMENT06','SEGMENT05','SEGMENT04',
        'SEGMENT03','SEGMENT02','SEGMENT01','SEGMENT00']
values = range(16)
segdict = dict(zip(keys, values))
#%matplotlib inline


# In[3]:

outfilename = thedir+'Correlations_Varying_Flux_5Apr16.pdf'
covsteps = 6
numsegments = 16
numfiles = 50
#fluxes = [24000.0,48000.0,72000.0,96000.0,120000.0]  # Get this right!!!
fluxes = [20000.0,40000.0,60000.0,80000.0,100000.0]  # Get this right!!!
flux_value = 80000.0
seqnos = [100,200,300,400,500]
numfluxes = len(fluxes)
covariance = zeros([covsteps, covsteps, numsegments, numfiles, numfluxes])
reduced_cov = zeros([covsteps,covsteps])
for i, seqno in enumerate(seqnos):
    infilename = "correlations_%d_20160107.txt"%seqno
    file = open(infilename,'r')
    lines = file.readlines()
    file.close
    for line in lines:
        items = line.split()
        if items[0] == 'ii':
            continue
        try:
            ii = int(items[0])
            jj = int(items[1])
            n = int(items[2])
            segment = segdict[items[3]]
            covariance[ii,jj,segment,n,i] = float(items[4]) 
        except:
            break
xvals = []
yvals = []
xfit = []
yfit = []
yerr = []
for ii in range(covsteps):
    for jj in range(covsteps):
        y = []
        for k in range(numfluxes):
            y.append(covariance[ii,jj,0:15,:,k].mean() / covariance[0,0,0:15,:,k].mean())
        """
        if ii == 2 and jj == 0:
            scatter(fluxes, y)
            slope, intercept, r_value, p_value, std_err = stats.linregress(fluxes,y)
            xplot=linspace(0.0, 120000.0, 100)
            yplot = slope * xplot + intercept
            plot(xplot,yplot,color='red')
            xlim(0,120000)
            xlabel("Flux(e-)")
            ylabel("Covariance (%d, %d)"%(jj,ii))
            show()
            sys.exit()
        """
        slope, intercept, r_value, p_value, std_err = stats.linregress(fluxes,y)
        reduced_cov[ii,jj] = slope * flux_value
        rsquared = float(ii*ii + jj*jj)
        if rsquared > 0.1:
            xvals.append(rsquared)
            yvals.append(reduced_cov[ii,jj])
            yerr.append((1.0 - sqrt(r_value)) * reduced_cov[ii,jj])
        if rsquared > 1.1 and reduced_cov[ii,jj] > 0.0:
            xfit.append(rsquared)
            yfit.append(reduced_cov[ii,jj])
        print "ii = %d, jj = %d, cov = %.4f, intercept = %.5f, R value = %.3f"%(ii,jj,reduced_cov[ii,jj], intercept, r_value)
yvals = array(yvals)
yerr = array(yerr)
ylower = np.maximum(1.1E-5, yvals - yerr)
yerr_lower = yvals - ylower
figure()
title("Correlation Coefficient %d Pairs of Flats - %d Electrons"%(numfiles,flux_value))
xscale('log')
yscale('log')
xlim(0.8,100.0)
ylim(1.0E-5,1.0E-1)
#scatter(xvals,yvals)
errorbar(xvals,yvals, yerr = [yerr_lower, 2.0*yerr] , ls = 'None',marker = '.', ms = 10, color = 'blue')
slope, intercept, r_value, p_value, std_err = stats.linregress(log10(xfit),log10(yfit))
xplot=linspace(0.0, 2.0, 100)
yplot = slope * xplot + intercept
plot(10**xplot, 10**yplot, color='red', lw = 2, ls = '--')
text(10.0, 0.005, "Slope = %.3f"%slope)
text(10.0, 0.00316, "C10 = %.5g"%reduced_cov[0,1])
text(10.0, 0.002, "C01 = %.5g"%reduced_cov[1,0])
xlabel("$i^2 + j^2$")
ylabel("Covariance")
#legend()
#show()
savefig(outfilename)
#close("All")


# In[ ]:



