
# coding: utf-8

# In[2]:

get_ipython().magic(u'cd ~')
import matplotlib
matplotlib.use("Agg")
import pyfits,glob,time,scipy
import scipy.interpolate
from scipy.special import erf
from pylab import *
from subprocess import call
from IPython import parallel
from scipy.optimize import fmin_powell

topdir='/Users/cslage/Research/LSST/code/GUI/'
thedir=topdir+'profiles/'
get_ipython().magic(u'cd $thedir')

configfile=topdir+'sextractor/default-array_dither.sex'
paramfile=topdir+'sextractor/default-array_dither.param'
maskdir=topdir+'sextractor/masks/'



# In[3]:

def remove_overscan_xy(image,x_overscan_start,x_overscan_end,y_overscan_start,y_overscan_end):
    overscan=image[:y_overscan_start,x_overscan_start+1:x_overscan_end]
    image=image[:y_overscan_start,:x_overscan_start]
    finalimg=image-matrix(median(overscan,axis=1)).T*np.ones(shape(image)[1])
    return array(finalimg)

def make_reg_from_ldac(cat_ldac_file,text_tag):
    thecat=pyfits.getdata(cat_ldac_file,'LDAC_OBJECTS')
    f = open(cat_ldac_file+'.reg','w')
    for i in range(len(thecat)):
        xcoo,ycoo=thecat['XWIN_IMAGE'][i],thecat['YWIN_IMAGE'][i]
        r=thecat['A_IMAGE'][i]
        thetext=thecat[text_tag][i]
        f.write('circle '+str(xcoo)+' '+str(ycoo)+' '+str(r)+'#text="'+str(thetext)+'"\n')
    f.close()
    
def Area(xl, xh, yl, yh, sigmax, sigmay, Imax):
    # Calculates how much of a 2D Gaussian falls within a rectangular box
    ssigx = sqrt(2) * sigmax
    ssigy = sqrt(2) * sigmay    
    I = (erf(xh/ssigx)-erf(xl/ssigx))*(erf(yh/ssigy)-erf(yl/ssigy))
    return Imax * I / 4.0


# In[5]:

# This definition runs sextractor in parallel and can be given to IPython's parallel map function along with
#   the file list
def ovsub_runsex_makereg(fitsfilename):
    import pyfits
    import numpy as np
    from subprocess import call
    topdir='/Users/cslage/Research/LSST/code/GUI/'
    thedir=topdir+'profiles'
    get_ipython().magic(u'cd $thedir')

    configfile=topdir+'sextractor/default-array_dither.sex'
    paramfile=topdir+'sextractor/default-array_dither.param'
    maskdir=topdir+'sextractor/masks/'

    def remove_overscan_xy(image,x_overscan_start,x_overscan_end,y_overscan_start,y_overscan_end):
        overscan=image[:y_overscan_start,x_overscan_start+1:x_overscan_end]
        image=image[:y_overscan_start,:x_overscan_start]
        finalimg=image-np.matrix(np.median(overscan,axis=1)).T*np.ones(np.shape(image)[1])
        return finalimg

    def make_reg_from_ldac(cat_ldac_file,text_tag):
        thecat=pyfits.getdata(cat_ldac_file,'LDAC_OBJECTS')
        f = open(cat_ldac_file+'.reg','w')
        for i in range(len(thecat)):
            xcoo,ycoo=thecat['XWIN_IMAGE'][i],thecat['YWIN_IMAGE'][i]
            r=thecat['A_IMAGE'][i]
            thetext=thecat[text_tag][i]
            f.write('circle '+str(xcoo)+' '+str(ycoo)+' '+str(r)+'#text="'+str(thetext)+'"\n')
        f.close()
    for i in range(1,17):
        extname=pyfits.getheader(fitsfilename,i)['EXTNAME']
        img=pyfits.getdata(fitsfilename,extname)
        overscansubimg=remove_overscan_xy(img,509,542,2000,2022)   # cut off the overscan
        outname=fitsfilename[:-5]+extname+'.fits'
        pyfits.writeto(outname,overscansubimg,clobber=True)
        test=call(["sex",outname,"-c",configfile,"-CATALOG_NAME",outname+'.cat'])
        make_reg_from_ldac(outname+'.cat','NUMBER')


# In[4]:

zfilelist=sort(glob.glob(thedir+'114-04_spot-30um_light_3??_20150709??????.fits')+
               glob.glob(thedir+'114-04_spot-30um_light_4??_201507????????.fits'))
print len(zfilelist)


# In[110]:

for fitsfilename in zfilelist: 
    tfile1=time.time() 
    for i in [4,5,6,12,14]:
        tstart=time.time() 
        extname=pyfits.getheader(fitsfilename,i)['EXTNAME'] 
        img=pyfits.getdata(fitsfilename,extname) 
        overscansubimg=remove_overscan_xy(img,509,542,2000,2022) 
        # cut off the overscan 
        outname=fitsfilename[:-5]+extname+'.fits' 
        pyfits.writeto(outname,overscansubimg,clobber=True) 
        test=call(["sex",outname,"-c",configfile,"-CATALOG_NAME",outname+'.cat']) 
        make_reg_from_ldac(outname+'.cat','NUMBER') 
        tend=time.time() 
        print extname+" time: "+str(tend-tstart)[:4] 
        print "Time taken for file "+str(fitsfilename[-26:-23])+": "+str(time.time()-tfile1)


# #### Using the sextractor catalogs produced above to make a map of the sextractor measurement named below

# In[37]:

class Array2dSet:
    def __init__(self,xmin,xmax,nx,ymin,ymax,ny,nstamps):
        # This packages up a set of nstamps postage stamp images,
        # each image of which is nx * ny pixels
        self.nx=nx
        self.ny=ny
        self.nstamps=nstamps

        self.xmin=xmin
        self.ymin=ymin
        
        self.xmax=xmax
        self.ymax=ymax
        
        self.dx=(xmax-xmin)/nx
        self.dy=(ymax-ymin)/ny
        
        self.x=linspace(xmin+self.dx/2,xmax-self.dx/2,nx)
        self.y=linspace(ymin+self.dy/2,ymax-self.dy/2,ny)

        self.data=zeros([nx,ny,nstamps])
        self.xoffset=zeros([nstamps])
        self.yoffset=zeros([nstamps])
        self.imax=zeros([nstamps])


def BuildSpotList(fitsfilename, segmentnumber, numspots, nx, ny, minsize, maxsize):
    stampxmin = -(int(nx/2)+0.5)
    stampxmax = -stampxmin
    stampymin = -(int(ny/2)+0.5)
    stampymax = -stampymin
    xcoomin = 50
    xcoomax = 450
    ycoomin = 1400
    ycoomax = 1900
    spotlist = Array2dSet(stampxmin,stampxmax,nx,stampymin,stampymax,ny,numspots)
    hdr=pyfits.getheader(fitsfilename,segmentnumber)
    extname = hdr['EXTNAME']
    img=pyfits.getdata(fitsfilename,extname) 
    overscansubimg=remove_overscan_xy(img,509,542,2000,2022) 
    catname=fitsfilename[:-5]+extname+'.fits.cat.reg' 
    catfile = open(catname,'r')
    catlines = catfile.readlines()
    catfile.close()
    n=0
    for line in catlines:
        try:
            size = float(line.split()[3].split('#')[0])
            if size < minsize or size > maxsize:
                continue
            xcoord = float(line.split()[1])
            ycoord = float(line.split()[2])
            if xcoord < xcoomin or xcoord > xcoomax or ycoord < ycoomin or ycoord > ycoomax:
                continue
            xint = int(xcoord-0.5)
            yint = int(ycoord-0.5)
            xmin = xint - int(nx/2)
            xmax = xint + int(nx/2) + 1
            ymin = yint - int(ny/2)
            ymax = yint + int(ny/2) + 1
            stamp = overscansubimg[ymin:ymax, xmin:xmax]
           
            xsum = 0.0
            ysum = 0.0
            datasum = 0.0
             
            for i in range(nx):
                for j in range(ny):
                    spotlist.data[i,j,n] = float(stamp[j,i])                    
                    ysum += spotlist.y[j] * spotlist.data[i,j,n]
                    xsum += spotlist.x[i] * spotlist.data[i,j,n]
                    datasum += spotlist.data[i,j,n]
            xoff = xsum / datasum
            yoff = ysum / datasum
            spotlist.xoffset[n] = xoff#xcoord - xint - 1.0
            spotlist.yoffset[n] = yoff#ycoord - yint - 1.0     
                    
            n += 1
            if n == numspots:
                return spotlist
        except:
            continue
    # Reaching this point means we found less spots than requested.
    newspotlist = Array2dSet(stampxmin,stampxmax,nx,stampymin,stampymax,ny,n)
    newspotlist.xoffset = spotlist.xoffset[0:n]
    newspotlist.yoffset = spotlist.yoffset[0:n]
    newspotlist.data = spotlist.data[:,:,0:n]
    del spotlist
    return newspotlist

import forward

def FOM(params):
    [sigmax, sigmay] = params
    result = forward.forward(spotlist,sigmax,sigmay)
    return result

def PyFOM(params):
    fom = 0.0
    [Imax, sigmax, sigmay] = params
    area=zeros([spotlist.nx,spotlist.ny])
   
    for spot in range(spotlist.nstamps):
        for ii in range(spotlist.nx):
            for jj in range(spotlist.ny):
                xl = spotlist.x[ii] - spotlist.xoffset[spot] - 0.5
                xh = xl + 1.0
                yl = spotlist.y[jj] - spotlist.yoffset[spot] - 0.5
                yh = yl + 1.0
                #if spot == 78 and ii == 4 and jj == 4 or spot==0:
                    #print "ii = %d, jj = %d, img = %.4f"%(ii,jj,spotlist.data[ii,jj,spot])
                #print "ii = %d, jj = %d,xl = %.2f, xh = %.2f, yl = %.2f, yh = %.2f"%(ii,jj,xl,xh,yl,yh)
                area[ii,jj] = Area(xl, xh, yl, yh, sigmax, sigmay, Imax)
                fom += square(area[ii,jj]-spotlist.data[ii,jj,spot])
    #print "Imax = %.1f, sigmax = %.2f, sigmay = %.2f, fom = %.1f"%(Imax, sigmax, sigmay, fom)
    return fom


# In[8]:

for i, segmentnumber in enumerate([4,5,6,12,14]):
    hdr=pyfits.getheader(zfilelist[0],segmentnumber)
    extname = hdr['EXTNAME']
    print "Segment # %d is %s"%(segmentnumber,extname)


# In[111]:

nx = 9
ny = 9
numspots =400
segmentnumber = 5

fitsfilename = zfilelist[10]
print fitsfilename
    
spotlist = BuildSpotList(fitsfilename, segmentnumber, numspots, nx, ny, 0.7, 1.4)

print "nstamps = %d"%spotlist.nstamps

param0 = [1.00, 1.00]

args = ()#(spotlist)
Result = fmin_powell(FOM, param0, args)

print Result



# In[112]:

import gc
nx = 9
ny = 9
numspots = 400
numfiles = len(zfilelist)
sigmaxs = zeros([5,numfiles])
sigmays = zeros([5,numfiles])
imaxs = zeros([5,numfiles])
totelectrons = zeros([5,numfiles])
gains = [4.89,5.03,5.04,4.94,5.03]

for i, segmentnumber in enumerate([4,5,6,12,14]):
    gain = gains[i]

    for j,fitsfilename in enumerate(zfilelist):
        param0 = [1.0, 1.0]
    
        spotlist = BuildSpotList(fitsfilename, segmentnumber, numspots, nx, ny,0.7,2.0)
        print "nstamps = %d"%spotlist.nstamps
        args = ()#(spotlist)
        Result = fmin_powell(FOM, param0, args)

        print fitsfilename
        imax = spotlist.imax.mean()
        ADU_correction = Area(-0.5,0.5,-0.5,0.5,Result[0],Result[1],1.0)
        tote = spotlist.data.sum() * gain / spotlist.nstamps # Total electrons in the stamp
        print Result, imax, tote
        imaxs[i,j] = imax * ADU_correction * gain
        sigmaxs[i,j] = abs(Result[0])
        sigmays[i,j] = abs(Result[1])
        totelectrons[i,j] = tote
        del spotlist
        gc.collect()


# In[113]:

print imaxs.shape, sigmaxs.shape, sigmays.shape
print imaxs[4,:]


# In[114]:

from scipy import stats
rcParams.update({'font.size':18})
# The following masks out bad fits files of undetermiend cause, probably due to the shutter sticking
#mask = [19, 20, 22] # Use this one for the 30V Vbb
#mask = [15,19,22,23,24,25,26,28,29,30,31] # Use this one for the 45V Vbb
mask = []#[18,19,21,22,24,25] # Use this one for the 60V Vbb

for i, segmentnumber in enumerate([4,5,6,12,14]):
    hdr=pyfits.getheader(zfilelist[0],segmentnumber)
    extname = hdr['EXTNAME']

    figure()
    title("Brighter-Fatter - 30 micron Spots - %s"%extname)
    #title("Brighter-Fatter - 3 micron Spots")
    scatter(delete(imaxs[i,:],mask), delete(sigmaxs[i,:],mask), color = 'green', lw = 2, label = 'Sigma-x')
    scatter(delete(imaxs[i,:],mask), delete(sigmays[i,:],mask), color = 'red', lw = 2, label = 'Sigma-y')

    slope, intercept, r_value, p_value, std_err = stats.linregress(imaxs[i,10:40],sigmaxs[i,10:40])

    xplot=linspace(-5000.0,300000.0,100)
    yplot = slope * xplot + intercept
    plot(xplot, yplot, color='blue', lw = 2, ls = '--')
    tslope = slope * 100.0 * 50000.0
    #text(10000.0,1.03,"X Slope = %.2f %% per 50K e-, Intercept = %.3f"%(tslope,intercept))
    text(2000.0,1.18,"X Slope = %.2f %% per 50K e-"%tslope, fontsize = 24)

    slope, intercept, r_value, p_value, std_err = stats.linregress(imaxs[i,10:40],sigmays[i,10:40])

    xplot=linspace(-5000.0,300000.0,100)
    yplot = slope * xplot + intercept
    plot(xplot, yplot, color='black', lw = 2, ls = '--')
    tslope = slope * 100.0 * 50000.0
    #text(10000.0,1.02,"Y Slope = %.2f %% per 50K e-, Intercept = %.3f"%(tslope,intercept))
    text(2000.0,1.16,"Y Slope = %.2f %% per 50K e-"%tslope, fontsize = 24)

    xlim(0.0,300000.0)
    xticks([0,100000,200000])
    ylim(0.90,1.20)
    xlabel('Central Peak(electrons)')
    ylabel('Sigma (Pixels)')
    legend(loc= 'lower right')
    savefig(datadir+"Forward_16Mar17_%s.png"%extname)


# In[115]:

segmentnumber = 5
hdr=pyfits.getheader(zfilelist[0],segmentnumber)
extname = hdr['EXTNAME']
gain = gains[0]

figure()
testn = 8
spot = 228
numspots = 400
fitsfilename = zfilelist[testn]
spotlist = BuildSpotList(fitsfilename, segmentnumber, numspots, nx, ny,0.2,1.4)
param0 = [1.0, 1.0]
args = ()#(spotlist)
Result = fmin_powell(FOM, param0, args)
sigmax = abs(Result[0])
sigmay = abs(Result[1])
ADU_correction = Area(-0.5,0.5,-0.5,0.5,sigmax, sigmay,1.0)
Imax = spotlist.imax[spot] 
seqno = int(zfilelist[testn].split("_")[4])
suptitle=("Sequence # %d, Spot # %d"%(seqno,spot))
area=zeros([nx,ny])
for ii in range(nx):
    for jj in range(ny):
        xl = ii - int(nx/2) - spotlist.xoffset[spot] - 0.5
        xh = xl + 1.0
        yl = jj - int(ny/2) - spotlist.yoffset[spot] - 0.5
        yh = yl + 1.0
        area[ii,jj] = Area(xl, xh, yl, yh, sigmax, sigmay, Imax) 

        
#loste = (area.sum() - spotlist.data.sum() / spotlist.nstamps) * gain
#loste = (area.sum() - spotlist.data[:,:,spot].sum() ) * gain
#print "Loste = %.0f"%loste
#outer_row = (spotlist.data[0,:,spot].sum() + spotlist.data[14,:,spot].sum() + spotlist.data[:,0,spot].sum() + 
#            spotlist.data[:,14,spot].sum()) * gain
#print "Electrons in outer row = %d"%outer_row
subplots_adjust(wspace = 2.0)
subplot(1,3,1)
imshow(spotlist.data[:,:,spot] ,interpolation="None")
subplot(1,3,2)
plot(spotlist.data[int(nx/2),:,spot] , lw = 2, label="Data")
plot(area[int(nx/2),:], lw = 2, label="Model")
xlabel("X (Pixels)")
ylabel("Signal(ADU)")
xticks([0,2,4,6,8])
subplot(1,3,3)
plot(spotlist.data[:,int(ny/2),spot], lw = 2, label="Data")
plot(area[:,int(ny/2)], lw = 2, label="Model")
xlabel("Y (Pixels)")
ylabel("Signal(ADU)")
xticks([0,2,4,6,8])
legend(loc = (-6.5,0.8))
savefig(datadir+"Typical_Fit_16Mar17_%s_%s"%(testn,extname))


# In[66]:

mask = [18,19,21,22,24,25]
masked_zfilelist = delete(zfilelist,mask)
print masked_zfilelist.shape


# In[116]:

#mask = []#[15,19,22,23,24,25,26,28,29,30,31]#[19,20,22]
from scipy import stats
figure()
exp_current=[]
masked_zfilelist = delete(zfilelist,mask)
for i,fitsfilename in enumerate(masked_zfilelist):
    try:
        exptime=pyfits.getheader(fitsfilename,0)['EXPTIME'] 
        mondiode=pyfits.getheader(fitsfilename,0)['MONDIODE'] 
        exp_current.append(float(exptime) * mondiode / (1.0E-9))
    except:
        continue
slope, intercept, r_value, p_value, std_err = stats.linregress(exp_current[0:12],imaxs[0,0:12])

xplot=linspace(0.0,3.0,100)
yplot = slope * xplot + intercept
scatter(exp_current,delete(imaxs[0,:],mask)/1.0E3)
plot(xplot, yplot/1.0E3, color='red')
text(5.0,60,"Linear Fit R^2 = %.5f"%r_value)
xlabel("Exp Time * Monitor Diode (nA-sec)")
ylabel("Modeled Peak Intensity (10^3 e-)")
#xlim(0.0, 3.0)
#ylim(0.0,350.0)
savefig(datadir+"Intensity_Check_16Mar17.png")


# In[95]:

from scipy import stats
mask=[]
figure()
title("Total Electrons vs Intensity * Exposure Time")
exp_current=[]
masked_zfilelist = delete(zfilelist,mask)
for i,fitsfilename in enumerate(masked_zfilelist):
    try:
        exptime=pyfits.getheader(fitsfilename,0)['EXPTIME'] 
        mondiode=pyfits.getheader(fitsfilename,0)['MONDIODE'] 
        exp_current.append(float(exptime) * mondiode / (1.0E-9))
    except:
        continue
slope, intercept, r_value, p_value, std_err = stats.linregress(exp_current[0:10],totelectrons[0,0:10])

xplot=linspace(0.0,exp_current[-1],100)
yplot = (slope * xplot + intercept)/1.0E6
scatter(exp_current,totelectrons[0,:]/1.0E6)
plot(xplot, yplot, color='red')
loste = yplot[-1] - (totelectrons[0,-1])/1.0E6
text(10.0,6.0,"%.2f Million e- lost"%loste)
text(5.0,1.0,"Linear Fit R^2 = %.5f"%r_value)
xlabel("Exp Time * Monitor Diode (nA-sec)")
ylabel("Total Electrons in Spot (10^6 e-)")
xlim(0.0, 25.0)
ylim(0.0,10.0)
savefig(datadir+"Total_Electrons_15.png")

#print totelectrons[0,29], imaxs[0,29], exp_current[29]
#print xplot[20:70], yplot[20:70]


InputElec = array(exp_current) * slope / 1.0E6
CollectedElec = array(totelectrons[0,:])/1.0E6
slope1, intercept1, r_value, p_value, std_err = stats.linregress(InputElec[2:8],CollectedElec[2:8])
xplot1=linspace(0.0,3.0,100)
yplot1 = slope1 * xplot1 + intercept1

slope2, intercept2, r_value, p_value, std_err = stats.linregress(InputElec[-4:-1],CollectedElec[-4:-1])
xplot2=linspace(1.0,4.0,100)
yplot2 = slope2 * xplot2 + intercept2

Input_Breakpoint = (intercept2 - intercept1) / (slope1 - slope2) * 1.0E6
Output_Breakpoint = intercept1 + slope1 * Input_Breakpoint
Slope1E = slope1
Slope2E = slope2

figure()
title("Electron Loss", fontsize = 16)
scatter(InputElec, CollectedElec, color = 'blue')
plot(xplot1, yplot1, color='red', lw = 2, ls = '--')
plot(xplot2, yplot2, color='green', lw = 2, ls = '--')
text(1.5,1.0,"Slope1 = %.2f"%slope1,fontsize=12)
text(1.5,0.8,"Slope2 = %.2f"%slope2,fontsize=12)
text(1.5,0.6,"Breakpoint = %.0f input electrons"%Input_Breakpoint,fontsize=12)
text(1.5,0.4,"Breakpoint = %.0f collected electrons"%Output_Breakpoint,fontsize=12)
text(1.5,0.2,"%.2f Million e- lost"%loste, fontsize = 12)
xlim(0.0,4.0)
ylim(0.0,4.0)
xlabel("Input Electrons (10^6 e-)")
ylabel("Collected Electrons (10^6 e-)")
savefig(datadir+"Collected_Electrons.png")

ExpCurrent = array(exp_current)
CollectedElec = array(totelectrons[0,:])/1.0E6
slope1, intercept1, r_value, p_value, std_err = stats.linregress(ExpCurrent[2:8],CollectedElec[2:8])
xplot1=linspace(0.0,15.0,100)
yplot1 = slope1 * xplot1 + intercept1

slope2, intercept2, r_value, p_value, std_err = stats.linregress(ExpCurrent[-4:-1],CollectedElec[-4:-1])
xplot2=linspace(1.0,15.0,100)
yplot2 = slope2 * xplot2 + intercept2

figure()
title("Electron Loss", fontsize = 16)
scatter(ExpCurrent, CollectedElec, color = 'blue')
plot(xplot1, yplot1, color='red', lw = 2, ls = '--')
plot(xplot2, yplot2, color='green', lw = 2, ls = '--')
text(2.5,1.0,"Slope1 = %.2f"%Slope1E,fontsize=12)
text(2.5,0.8,"Slope2 = %.2f"%Slope2E,fontsize=12)
text(2.5,0.4,"Breakpoint = %.0f collected electrons"%Output_Breakpoint,fontsize=12)
text(2.5,0.2,"%.2f Million e- lost"%loste, fontsize = 12)
xlim(0.0,8.0)
ylim(0.0,4.0)
xlabel("Exp Time * Monitor Diode (nA-sec)")
ylabel("Collected Electrons (10^6 e-)")
savefig(datadir+"Collected_Electrons_2.png")


# In[ ]:



