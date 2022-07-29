
# This creates a text file with the measured covariances
from pylab import *
import os, sys, glob, time, scipy
from scipy import stats
import eups
from lsst.daf.persistence import Butler


REPO_DIR = '/home/cslage/Research/LSST/code/notebooks/notebooks_2019_02_06/stand_alone_13feb19/flats_repo_bip_quad_zero_large_bias'
bfk_butler = Butler(REPO_DIR+'/rerun/test')
# Note ITL uses raftName R02 and detector 2
# E2V uses raftName R00 and detector 0
bf_kernel = bfk_butler.get('brighterFatterKernel', dataId={'raftName': 'R00', 'detectorName': 'S00',
                                                              'detector': 0})
means = bf_kernel.means
xcorrs = bf_kernel.xcorrs
meanXcorrs = bf_kernel.meanXcorrs

amp = 'C14'
flux = 18

NumPairs = 5
NumFluxes = 20
amp = 'C14'

filename = "/home/cslage/Research/LSST/code/notebooks/notebooks_2019_02_06/stand_alone_13feb19/e2v_covariances_21oct19.txt"
file = open(filename, 'w')
line = "i \t j \t Cov mean \t Cov std \n"
file.write(line)
for jj in range(17):
    for ii in range(17):
        r2 = (ii-8)**2 + (jj-8)**2
        value = meanXcorrs[amp][ii,jj]
        n = NumFluxes - 1
        cov = []
        for m in range(NumPairs):
            i = n * NumPairs + m
            xcorr = xcorrs[amp][i][abs(ii-8),abs(jj-8)]
            cov.append(xcorr)
        cov = array(cov)
        if ii >= 8 and jj >= 8:
            line = f"{jj-8} \t {ii-8} \t {-value:#6.4g} \t {cov.std():#6.4g} \n"
            file.write(line)
            #print(line)
file.close()
