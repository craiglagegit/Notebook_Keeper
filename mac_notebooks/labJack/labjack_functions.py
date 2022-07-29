'''
This code has a set of functions to read the labJack T7-Pro,
which is interfaced to 3 - 3-axis accelerometers
Uses 9 analog inputs (AINs) to read the data at 200 Hz.
This notebook reads the accelerometers while taking an image
To run it, you will first have to do the following:

(1) pip install labjack-ljm.  This will build the labJack Python code
    in your local directory at ~/.local/lib/python3.8/site-packages/labjack \
(2) In your ~/notebooks/.user_setups file, add the following line: \
      export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/scratch/lib

To use from an observing notebook:
from labjack_functions import LabJackFunctions
# Other stuff to create a latiss instance
lj = LabJackFunctions(latiss)
await lj.take_object_with_accel_data(expTime, filter='empty_1', grating='empty_1')
lj.close()

Craig Lage - Mar 10, 2022
'''
import sys, time, datetime, asyncio
from subprocess import Popen
import numpy as np
import pandas as pd
import pickle as pkl
from labjack import ljm

# Structures to hold the calibration data and connection information
class zeroOffset:
    def __init__(self, off_x, off_y, off_z):
        # This is the reading in Volts when accel = 0
        self.x = off_x
        self.y = off_y
        self.z = off_z
        
class gMult:
    def __init__(self, mult_x, mult_y, mult_z):
        # This is the conversion to acceleration in V/g
        self.x = mult_x
        self.y = mult_y
        self.z = mult_z
        
class AIN_name:
    def __init__(self, ain_x, ain_y, ain_z):
        # This is where the sensors are connected to the labJack
        self.x = ain_x
        self.y = ain_y
        self.z = ain_z

class calData:
    def __init__(self, serial="", ain_x="", ain_y="", ain_z="", off_x=0.0, off_y=0.0, off_z=0.0, mult_x=1.0, mult_y=1.0, mult_z=1.0):
        # The serial number is imprinted on the accelerometer
        self.serial = serial
        self.AIN_name = AIN_name(ain_x, ain_y, ain_z)
        self.zeroOffset = zeroOffset(off_x, off_y, off_z)
        self.gMult = gMult(mult_x, mult_y, mult_z)

# This holds the necessary functions
class LabJackFunctions:
    def __init__(self, latiss):
        # This sets up the labJack for reading data
        self.latiss = latiss

        # First, the calibration info
        calDict = {}
        calDict["1"] = calData(serial="A395429", ain_x="AIN1", ain_y="AIN2", ain_z="AIN3", off_x=2.49017, off_y=2.44424, off_z=2.44589, mult_x=0.98959, mult_y=0.98572, mult_z=0.99946)
        calDict["2"] = calData(serial="A395423", ain_x="AIN4", ain_y="AIN5", ain_z="AIN6", off_x=2.49874, off_y=2.49595, off_z=2.41423, mult_x=0.99740, mult_y=1.00142, mult_z=0.99595)
        calDict["3"] = calData(serial="A395446", ain_x="AIN7", ain_y="AIN8", ain_z="AIN9", off_x=2.47830, off_y=2.48088, off_z=2.41385, mult_x=0.97957, mult_y=0.98699, mult_z=1.00376)

        # Now open communications with the labJack
        self.handle = ljm.openS("T7", "wifi", "139.229.164.249")  
        # Verify that you are communicating
        info = ljm.getHandleInfo(self.handle)
        print("Opened a LabJack with Device type: %i, Connection type: %i,\n"
              "Serial number: %i, IP address: %s, Port: %i,\nMax bytes per MB: %i" %
              (info[0], info[1], info[2], ljm.numberToIP(info[3]), info[4], info[5]))

        # Ensure triggered stream is disabled.
        ljm.eWriteName(self.handle, "STREAM_TRIGGER_INDEX", 0)
        # Enabling internally-clocked stream.
        ljm.eWriteName(self.handle, "STREAM_CLOCK_SOURCE", 0)

        # Set up analog parameters
        aRange = 10.0 # +/- 10.0 Volts
        aSettle = 0 # 0 microsecond settling time
        resIndex = 0
        aNames = ["AIN_ALL_NEGATIVE_CH", "STREAM_SETTLING_US", "STREAM_RESOLUTION_INDEX"] # List of set-up parameters
        aValues = [ljm.constants.GND, aSettle, resIndex] # List of set-up values

        self.aScanListNames = [] # List of AIN names which will be read
        self.offsets = []
        self.gMults = []
        for name in ["1", "2", "3"]:
            for axis in ["x", "y", "z"]:
                print(f"aName = calDict['{name}'].AIN_name.{axis}")
                print(calDict['1'].AIN_name.x)
                exec(f"aName = calDict['{name}'].AIN_name.{axis}", locals()) 
                
                print(aName)
 
                self.aScanListNames.append(aName)
                aNames.append(aName+"_RANGE")
                aValues.append(aRange)
                exec(f"off = calDict['{name}'].zeroOffset.{axis}")
                self.offsets.append(off)
                exec(f"gMult = calDict['{name}'].gMult.{axis}")
                self.gMults.append(gMult)

        self.offsets = np.array(self.offsets)
        self.gMults = np.array(self.gMults)

        # Write the analog inputs' negative channels (when applicable), ranges,
        # stream settling time and stream resolution configuration.
        numFrames = len(aNames)
        ljm.eWriteNames(self.handle, numFrames, aNames, aValues)
        return

    def __del__(self):
        ljm.close(self.handle)        
        return

    def close(self):
        ljm.close(self.handle)        
        return

    def readStream(self, readTime, scanRate=200):
        # This reads the accelerometers for a time readTime
        # and returns a Pandas timeSeries with the data
        # aScanListNames is the list of AIN ports (9 in total)
        # handle is the handle for talking to the labJack
        # scanRate is the read frequency in Hertz
        # readTime is the total time of read in seconds
        # calDict is the dictionary with the calibration data
        # The function returns a Pandas dataframe with the three 
        # accelerometers times three axes results

        numAddresses = len(self.aScanListNames)
        aScanList = ljm.namesToAddresses(numAddresses, self.aScanListNames)[0]
        scansPerRead = int(scanRate * readTime)
        # Configure and start stream
        scanRate = ljm.eStreamStart(self.handle, scansPerRead, numAddresses, aScanList, scanRate)
        start = datetime.datetime.now()
        # Stream the data
        ret = ljm.eStreamRead(self.handle)
        # Stop the stream
        ljm.eStreamStop(self.handle)
        aData = ret[0]
        # Reshape the data
        newData = np.resize(aData, (scansPerRead, numAddresses))
        # Convert to g
        accelData = (newData - self.offsets) / self.gMults
        # Create the timestamps
        end = start + datetime.timedelta(seconds = readTime)
        date_rng = pd.date_range(start=start, end=end, periods=scansPerRead)
        # Create the Pandas dataFrame
        df = pd.DataFrame(accelData, index=date_rng, 
                          columns=['ELM2', 'AZM2', 'ZM2', 'ELT', 'ZT', 'AZT', 'ELM1', 'AZM1', 'ZM1'])
        # Pickle the dataframe.  Use a temporary filename,
        # then update it after the exposure has finished.
        file = open(f'/scratch/labJackData/temp.pkl', 'wb')
        pkl.dump(df, file)
        file.close()

        finish = datetime.datetime.now()
        print(f"Finishing Accel data: {finish}")
        return 

    async def waitForAccelData(self, readTime):
        # This function wraps the readStream in an async function
        start = datetime.datetime.now()
        print(f"Starting Accel data: {start}")
        loop = asyncio.get_event_loop()
        await loop.run_in_executor(None, lambda:readStream(readTime))
        return

    async def objectExposure(self, expTime, filter='empty_1', grating='empty_1'):
        # Takes an exposure
        start = datetime.datetime.now()
        print(f"Starting exposure: {start}")
        exp_ids = await self.latiss.take_object(expTime, 1, filter=filter, grating=grating)
        finish = datetime.datetime.now()
        print(f"Exposure done: {finish}")
        return exp_ids

    async def take_object_with_accel_data(expTime, filter='empty_1', grating='empty_1'):
        # This puts them together
        readTime = expTime + 8.0 # Add some buffer for the accel data
        # This buffer can be reduced after astropy tai problem is fixed
        result = await asyncio.gather(waitForAccelData(readTime, expId), objectExposure(expTime, filter=filter, grating=grating))
        print(result)
        # Need to extract expId from result
        """
        # Stuff to rename file
        command = f"mv /scratch/labJackData/temp.pkl /scratch/labJackData/{expId}.pkl"
	changeName = Popen(command, shell=True)
	Popen.wait(changeName)
        print(f"Wrote file /scratch/labJackData/{expId}.pkl")
        """
        return

