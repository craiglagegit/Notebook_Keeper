#!/usr/bin/env python
# coding: utf-8

# ## Reading dome monitor with a labJack T7
# Uses 2 analog input (AINs) to read the data.
# 
# Craig Lage - Aug 17, 2023

# Need to set the LD_LIBRARY_PATH variable to find the labJack .so library
import os, sys, time, datetime, argparse
import numpy as np
import pandas as pd
import pickle as pkl
import matplotlib.pyplot as plt
from astropy.time import Time, TimeDelta
sys.path.append("/scratch/cslage/labJack")
from labjack import ljm  # Needed pip install labjack-ljm

########## SUBROUTINES ##################

def readStream(handle, names, readTime, scanRate=100):
    # This reads the accelerometers for a time readTime                                                             
    # and returns a Pandas timeSeries with the data                                                                 
    # names is the list of AIN ports (9 in total)                                                          
    # handle is the handle for talking to the labJack                                                               
    # scanRate is the read frequency in Hertz                                                                       
    # readTime is the total time of read in seconds                                                                 

    numAddresses = len(names)
    aScanList = ljm.namesToAddresses(numAddresses, names)[0]
    scansPerRead = int(scanRate * readTime)
    # Configure and start stream 
    # First close any existing streams
    # There shouldn't be any existing streams, but 
    # there have been some unexplained errors
    try:
        ljm.eStreamStop(handle)
    except ljm.LJMError:
        pass
    streamStartSuccess = False
    for i in range(5):
        # This sometimes fails, so try multiple times.
        try:
            scanRate = ljm.eStreamStart(handle, scansPerRead, numAddresses, aScanList, scanRate)
            streamStartSuccess = True
            break
        except:
            print(f"StreamStart # {i+1} failed")
            sys.stdout.flush()
            time.sleep(5.0)
    if not streamStartSuccess:
        print(f"Stream start failed!")
        sys.stdout.flush()
        return
        

    start = datetime.datetime.now()
    # Stream the data                                                                                               
    # Added the try-except block because it has failed occasionally
    try:
        ret = ljm.eStreamRead(handle)
    except ljm.LJMError:
        print(f"Stream read failed!")
        sys.stdout.flush()
        return
    # Stop the stream                                                                                               
    ljm.eStreamStop(handle)
    aData = ret[0]
    # Reshape the data                                                                                              
    newData = np.resize(aData, (scansPerRead, numAddresses))
    end = start + datetime.timedelta(seconds = readTime)
    date_rng = pd.date_range(start=start, end=end, periods=scansPerRead)
    # Create the Pandas dataFrame                                                                                   
    df = pd.DataFrame(newData, index=date_rng,
                      columns=names)
    # Pickle the dataframe.  Use a temporary filename,                                                              
    # then update it after the exposure has finished.   
    # Stuff to rename file                                                                                                                 
    name = start.strftime("%Y-%m-%dT%H:%M:%S")
    name = name.replace('-','').replace(':','')
    name = name + 'Z'
    filename = f'/scratch/labJackData/Dome_{name}.pkl'
    file = open(filename, 'wb')
    pkl.dump(df, file)
    file.close()

    finish = datetime.datetime.now()
    print(f"Finishing Accel data: {finish}")
    sys.stdout.flush()
    return

########## MAIN PROGRAM ##################

parser = argparse.ArgumentParser()
parser = argparse.ArgumentParser(add_help=False)
parser.add_argument('--end_time', help='Time to finish.  Format = 2023-08-17T016:00:00')
parser.add_argument('-h', '--help', action='help', default=argparse.SUPPRESS,
                    help='To run this in the background, use the following syntax: nohup python Dome_Monitor_Test_17Aug23.py --end_time 2023-08-17T10:00:00 >& monitor_log_17aug23.txt&')
args = parser.parse_args()
end_time = args.end_time
finish_time = Time(end_time, format='isot')

# Open LabJack T7
handle = ljm.openS("T7", "wifi", "139.229.170.164") 
info = ljm.getHandleInfo(handle)
print("Opened a LabJack with Device type: %i, Connection type: %i,\n"
      "Serial number: %i, IP address: %s, Port: %i,\nMax bytes per MB: %i" %
      (info[0], info[1], info[2], ljm.numberToIP(info[3]), info[4], info[5]))

names = ["AIN0", "AIN1", "AIN2"]

while (time.time() - finish_time.unix_tai) < 0.0:
    readStream(handle, names, 300, scanRate=100)
print("Finished with loop")
sys.stdout.flush()

ljm.close(handle)

