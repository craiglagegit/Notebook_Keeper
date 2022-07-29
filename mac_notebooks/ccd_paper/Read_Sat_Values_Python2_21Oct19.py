# This creates a text file with the measured covariances
from pylab import *
import os, sys, glob, time, scipy
import pickle


filename = "/Users/cslage/Research/LSST/code/poisson/Poisson_CCD_Hole20_Test11E/measurements/sat_data_20170509.pkl"
infile = open(filename, 'rb')
myList = pickle.load(infile)
infile.close()

outfilename = 'sat_data_21Oct19.txt'
outfile = open(outfilename, 'w')
line = "Vpl(V) \t Vph(V) \t Sat level (e-) \n"
outfile.write(line)

Vpls = [-8.0, -6.0, -4.0, -2.0]

for i in range(4):
    Vpl = Vpls[i]
    [Vphs, SatLevels] = myList[i + 2]
    for j, Vph in enumerate(Vphs):
        line = "%.1f \t %.1f \t %.1f \n"%(Vpl, Vph, SatLevels[j])
        outfile.write(line)
outfile.close()
