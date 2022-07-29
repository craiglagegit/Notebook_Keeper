
# coding: utf-8

# In[1]:

#import matplotlib
#matplotlib.use("PDF")
from pylab import *
import sys, time


# In[2]:

import numpy as np
electrons_per_pixel = 10000
NumTrials = 100
Nx = 1000
Ny = 1000
pixels = np.random.poisson(lam = electrons_per_pixel, size=(Nx,Ny,NumTrials))
print(pixels.shape)
        


# In[3]:

covsteps = 6
cov = zeros([covsteps, covsteps, NumTrials/2])
for trial in range(0, NumTrials, 2):
    fdiff = pixels[:,:,trial] - pixels[:,:,trial+1]
    (nrows,ncols)=fdiff.shape

    for k in range(0,covsteps):
        for l in range(0,covsteps):
            npixused1=0
            data1=fdiff[l:nrows  ,   k:ncols]
            data2=fdiff[0:nrows-l  , 0:ncols-k]
            npixused1=data1.size
            sum11=data1.sum()
            sum21=data2.sum()
            sum121=(data1*data2).sum()
            corr = (sum121 - sum11*sum21/float(npixused1))/float(npixused1)
            cov[k,l,trial/2] = corr
            #print("For trial %d, NumDataPoints = %d, ii = %d, jj = %d,Cij = %f"%(trial, npixused1, l, k, corr))


# In[4]:

print(cov[0,0,:].mean(),cov[0,0,:].std())
print(cov[1:-1,1:-1,:].mean(),cov[1:-1,1:-1,:].std())


# In[5]:

figure(figsize=(16,8))
subplot(1,2,1)
hist(cov[0,0,:])
subplot(1,2,2)
hist(cov[1:-1,1:-1,:].flatten())
show()


# In[ ]:

# Now add some correlations


# In[10]:

import numpy as np
electrons_per_pixel = 10000
NumTrials = 1000
Nx = 1000
Ny = 1000
pixels = np.random.poisson(lam = electrons_per_pixel, size=(Nx,Ny,NumTrials))
print(pixels.shape)
        


# In[11]:

fraction = 0.01
shift = fraction * (pixels - electrons_per_pixel)
delta =  - shift + np.roll(shift,1,axis=0)
print delta.sum()
print(shift[0:10,0,0])
print(np.roll(shift,1,axis=0)[0:10,0,0])
pixels = pixels + delta


# In[12]:

covsteps = 6
cov = zeros([covsteps, covsteps, NumTrials/2])
for trial in range(0, NumTrials, 2):
    fdiff = pixels[:,:,trial] - pixels[:,:,trial+1]
    (nrows,ncols)=fdiff.shape

    for k in range(0,covsteps):
        for l in range(0,covsteps):
            npixused1=0
            data1=fdiff[l:nrows  ,   k:ncols]
            data2=fdiff[0:nrows-l  , 0:ncols-k]
            npixused1=data1.size
            sum11=data1.sum()
            sum21=data2.sum()
            sum121=(data1*data2).sum()
            corr = (sum121 - sum11*sum21/float(npixused1))/float(npixused1)
            cov[k,l,trial/2] = corr
            #print("For trial %d, NumDataPoints = %d, ii = %d, jj = %d,Cij = %f"%(trial, npixused1, l, k, corr))


# In[13]:

print(cov[0,0,:].mean(),cov[0,0,:].std())
print(cov[0,1].mean(),cov[0,1].std())
print(cov[0,2].mean(),cov[0,2].std())
print(cov[0,3].mean(),cov[0,3].std())
print(cov[0,4].mean(),cov[0,4].std())
var = 0.0
for trial in range(NumTrials/2):
    for ii in range(-covsteps+1,covsteps):
        for jj in range(-covsteps+1,covsteps):
            var += cov[abs(ii),abs(jj),trial]

print(var/(NumTrials/2))


# In[ ]:

# Now add in some noise


# In[ ]:

import numpy as np
electrons_per_pixel = 10000
NumTrials = 1000
Nx = 1000
Ny = 1000
NoiseMean = 0.0
NoiseSigma = 10.0
pixels = np.random.poisson(lam = electrons_per_pixel, size=(Nx,Ny,NumTrials))+ np.random.normal(loc=NoiseMean,scale=NoiseSigma,size=(Nx,Ny,NumTrials))
print(pixels.shape)
        


# In[ ]:

fraction = 0.01
shift = fraction * (pixels - electrons_per_pixel)
delta =  - shift + np.roll(shift,1,axis=0)
print delta.sum()
print(shift[0:10,0,0])
print(np.roll(shift,1,axis=0)[0:10,0,0])
pixels = pixels + delta


# In[ ]:

covsteps = 6
cov = zeros([covsteps, covsteps, NumTrials/2])
for trial in range(0, NumTrials, 2):
    fdiff = pixels[:,:,trial] - pixels[:,:,trial+1]
    (nrows,ncols)=fdiff.shape

    for k in range(0,covsteps):
        for l in range(0,covsteps):
            npixused1=0
            data1=fdiff[l:nrows  ,   k:ncols]
            data2=fdiff[0:nrows-l  , 0:ncols-k]
            npixused1=data1.size
            sum11=data1.sum()
            sum21=data2.sum()
            sum121=(data1*data2).sum()
            corr = (sum121 - sum11*sum21/float(npixused1))/float(npixused1)
            cov[k,l,trial/2] = corr
            #print("For trial %d, NumDataPoints = %d, ii = %d, jj = %d,Cij = %f"%(trial, npixused1, l, k, corr))


# In[ ]:

print(cov[0,0,:].mean(),cov[0,0,:].std())
print(cov[0,1].mean(),cov[0,1].std())
print(cov[0,2].mean(),cov[0,2].std())
print(cov[0,3].mean(),cov[0,3].std())
print(cov[0,4].mean(),cov[0,4].std())
var = 0.0
for trial in range(NumTrials/2):
    for ii in range(-covsteps+1,covsteps):
        for jj in range(-covsteps+1,covsteps):
            var += cov[abs(ii),abs(jj),trial]

print(var/(NumTrials/2))


# In[ ]:



