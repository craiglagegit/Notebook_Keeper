
# coding: utf-8

# In[9]:

import matplotlib
matplotlib.use("PDF")
from pylab import *
import sys, glob, xlrd, datetime
from scipy import stats
thedir='/Users/cslage/Research/LSST/code/poisson/Poisson_CCD_Hole17/'
get_ipython().magic(u'cd $thedir')


# In[4]:

def mu_si(E, T):
    vm = 1.53E9 * pow(T, -0.87)
    Ec = 1.01 * pow(T, 1.55)
    beta = 2.57E-2 * pow(T, 0.66)
    return ((vm/Ec) / pow(1.0 + pow(abs(E)/Ec,beta), 1/beta))


# In[5]:

def ReadConfigFile(filename):
    # This reads the Poisson simulator config file for
    # the settings that were run
    # and returns a dictionary with the values
    ConfigData = {}
    try:
        file = open(filename,'r')
        lines=file.readlines()
        file.close()
    except IOError:
        print "Configuration file %s not found"%filename
        return False, ConfigData 

    try:
        for line in lines:
            ThisLine=line.strip().split()
            ThisLineLength=len(ThisLine)
            if ThisLineLength < 3:
                continue
            if list(ThisLine[0])[0]=='#' or ThisLine[0]=='\n':
                continue
            try:
                ParamName = ThisLine[0]
                ThisLine.remove(ThisLine[0])
                for counter,item in enumerate(ThisLine):
                    if list(item)[0] == '#':
                        del ThisLine[counter:] # Strip the rest of the line as a comment
                        continue
                    if item == '=':
                        ThisLine.remove(item)
                        continue
                if len(ThisLine) == 0:
                    continue
                elif len(ThisLine) == 1:
                    ThisParam = ThisLine[0]
                    try: ConfigData[ParamName] = int(ThisParam)
                    except ValueError:
                        try:
                            ConfigData[ParamName] = float(ThisParam)
                        except ValueError:
                            try:
                                ConfigData[ParamName] = ThisParam
                            except ValueError:
                                return False, ConfigData 
                else:
                    ThisParam = []
                    for item in ThisLine:
                        try: ThisParam.append(int(item))
                        except ValueError:
                            try: ThisParam.append(float(item))
                            except ValueError:
                                ThisParam.append(item)
                    ConfigData[ParamName] = ThisParam
            except (IOError, ValueError):
                continue
    except Exception as e:
        print "Error reading configuration file %s. Exception of type %s and args = \n"%(filename,type(e).__name__), e.args 
        return False, ConfigData 

    return True, ConfigData



# In[10]:

q = 1.6E-19   # MKS
me = 9.11E-31   # MKS
k = 1.38E-23   # MKS
T = 273.0
E=1000.0
Tox = 1.0E-5
Eps_ox = 8.85E-14 * 4.2
mu = mu_si(E, T)
print "Mobility = %f"%mu
Vds = 0.5
L = 5.0E-4
Vgates = []
Is = []
#dirs = ['A','B','C','D','G','E','H','F','I','J']
dirs = ['AF','BF','CF','KF','DF','GF','EF','HF','FF','IF','JF']
#dirs = ['HS','FS']
#dirs = ['EF', 'JF']
for dir in dirs:
    cfg_success, ConfigData = ReadConfigFile('data/transrun1%s/trans.cfg'%dir)
    if not cfg_success:
        print "Configuration file issue. Quitting"
        sys.exit()

    Vgate = ConfigData["FixedRegionVoltage_11"]
    file = open('data/transrun1%s/charge.txt'%dir,'r')
    lines = file.readlines()
    file.close
    Ne = float(lines[0].split()[5].strip(','))
    print dir, Vgate, Ne
    I = Vds * q * mu * Ne / L**2
    Vgates.append(Vgate)
    Is.append(I)
Vgates0 = []
Is0 = []
dirs = ['AG','BG','CG','KG','DG','GG','EG','HG','FG','IG','JG']
for dir in dirs:
    cfg_success, ConfigData = ReadConfigFile('data/transrun1%s/trans.cfg'%dir)
    if not cfg_success:
        print "Configuration file issue. Quitting"
        sys.exit()

    Vgate = ConfigData["FixedRegionVoltage_11"]
    file = open('data/transrun1%s/charge.txt'%dir,'r')
    lines = file.readlines()
    file.close
    Ne = float(lines[0].split()[5].strip(','))
    print dir, Vgate, Ne
    I = Vds * q * mu * Ne / L**2
    Vgates0.append(Vgate)
    Is0.append(I)


# In[11]:

thedir='/Users/cslage/Research/LSST/optical_simulator/timing/ltspice/sta3800_output_device'
get_ipython().magic(u'cd $thedir')
data_wb = xlrd.open_workbook('STA3800_meas.xls')
idvg1_data = data_wb.sheet_by_name('8-1011.TXT')
idvg3_data = data_wb.sheet_by_name('8-1014.TXT')


# In[20]:

# ID-VG 1
W= 28.0
L = 5.0
Mu = 1000.0
Vbs = 0.0
Vds = 0.50
Vgs = []
Ids = []
for i in range(idvg1_data.nrows):
    try:
        if type(idvg1_data.row(i)[0].value) is float:
            #if idvg1_data.row(i)[1].value > 0.0:
            Vgs.append(idvg1_data.row(i)[1].value)
            Ids.append(idvg1_data.row(i)[3].value)
    except:
        continue
#Qss = 1.5E12
#Qss = 0.0
DeltaVt = - 5.0
figure()
ax1=axes([0.2,0.1,0.6,0.6])
ax1.set_title("STA3800C Id-Vg")
ax1.plot(Vgs, array(Ids)*1000.0, color = 'green', lw = 2, label='Measured')
ax1.scatter(array(Vgates), array(Is)*1000.0, marker = 'x', color='red', label='Sim - QS=1.4E12')
ax1.scatter(array(Vgates0), array(Is0)*1000.0, marker = '+', color='blue', label='Sim - QS=0.0')
#ax1.scatter(array(Vgates_100), array(Is_100)*1000.0, marker = '*', color='blue', label='Sim:Thick=100')
#ax1.scatter(array(Vgates) + DeltaVt, array(Is)*1000.0, marker = '*', color='blue', label='$\delta V$=%.1f'%DeltaVt)
ax1.set_xlabel("Vgs (volts)")
ax1.set_ylabel("Ids(mA)")
ax1.set_ylim(0,1.0)
ax1.legend(loc='upper left')
ax1.text(-23, 0.3, 'Vds = 0.5V')
#ax1.text(-23, 0.3, 'DeltaV = %f'%DeltaVt)
#show()
savefig("IdVg_QS_no_QS_18May17.pdf")


# In[25]:

print "DVt = %.3f"%(7.5E11 * q / Eps_ox * Tox)


# In[60]:

# ID-VG 3
W= 28.0
L = 5.0
Mu = 1000.0
Vds = 5.0
Vbs = -8.0
Vgs_num_steps = 201
Vgs = zeros(Vgs_num_steps)
Ids = zeros([Vgs_num_steps])
for i in range(idvg3_data.nrows):
    try:
        if type(idvg3_data.row(i)[0].value) is float:
            vgs_index = (int(idvg3_data.row(i)[0].value) - 1) % Vgs_num_steps
            vbs_index = 1
            Vgs[vgs_index] = idvg3_data.row(i)[1].value
            Ids[vgs_index] = idvg3_data.row(i)[3].value
    except:
        continue

DeltaVt = - 5.0
fig = figure()
ax1=axes([0.2,0.1,0.6,0.6])
ax1.set_title("STA3800C Id-Vg")
ax1.plot(Vgs, array(Ids)*1000.0, color = 'green', lw = 2, label='Measured')
ax1.scatter(array(Vgates), array(Is2)*1000.0, marker = '*', color='red', label='Simulation')
ax1.scatter(array(Vgates) + DeltaVt, array(Is2)*1000.0, marker = '*', color='blue', label='$\delta V$=%.1f'%DeltaVt)
ax1.set_xlabel("Vgs (volts)")
ax1.set_ylabel("Ids(mA)")
ax1.set_ylim(0,8.0)
ax1.legend(loc='upper left')
ax1.text(-23, 0.3, 'Vds = 5.0V')
savefig("IdVg3_5May17.pdf")


# In[53]:

dq = 1.0 * Eps_ox / Tox / q

print "DQ = %g"%dq


# In[ ]:



