#!/usr/bin/env python
#Author: Craig Lage
#Date: 24-Apr-15
# Constants

from pylab import *
class CGS:
    def __init__(self):
        self.c = 2.99792E10
        self.pc = 3.0857E18
        self.kpc = self.pc * 1000.0
        self.Mpc = self.kpc * 1000.0
        self.yr = 365.25 * 86400
        self.Gyr = self.yr * 1.0E9
        self.ly = self.c * self.yr
        self.H0 = 67.8 * 1.0E5 / self.Mpc
        self.h =  6.62607015E-27
        self.hbar = self.h / (2.0 * pi)
        self.G = 6.67428E-8
        self.kb = 1.38065E-16
        self.Msun = 1.9891E33
        self.mp = 1.673E-24
        self.re = 2.818E-13
        self.alpha = 1.0 / 137.0
        self.e = 4.80320425E-10
        self.q = 1.602176565E-19
        self.me = 9.11E-28
        self.AU = 1.496E13
        self.Jy = 1.0E-23
        self.Rsun = 6.955E10
        self.eV = 1.602176565E-12
        self.keV = self.eV * 1000.0
        self.eps0 = 8.85419E-10        
        return

class MKS:
    def __init__(self):
        self.c = 2.99792E8
        self.pc = 3.0857E16
        self.kpc = self.pc * 1000.0
        self.Mpc = self.kpc * 1000.0
        self.yr = 365.25 * 86400
        self.Gyr = self.yr * 1.0E9
        self.ly = self.c * self.yr
        self.H0 = 67.8 * 1.0E3 / self.Mpc
        self.h =  6.62607015E-34
        self.hbar = self.h / (2.0 * pi)
        self.G = 6.67428E-11
        self.kb = 1.38065E-23
        self.Msun = 1.9891E30
        self.mp = 1.673E-27
        self.re = 2.818E-15
        self.alpha = 1.0 / 137.0
        self.q = 1.602176565E-19
        self.me = 9.11E-31
        self.AU = 1.496E11
        self.Jy = 1.0E-26
        self.Rsun = 6.955E8
        self.eV = 1.602176565E-19
        self.keV = self.eV * 1000.0
        self.mu0 = 4.0 * pi * 1.0E-7
        self.eps0 = 8.85419E-12
        return


CGS = CGS()
MKS = MKS()
        
