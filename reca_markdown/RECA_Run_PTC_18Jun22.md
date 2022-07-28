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

## Running the PTC task from a Notebook

Craig Lage - 18-Jun-22

```python
import time
import numpy as np
import matplotlib.pyplot as plt
import astropy.io.fits as pf
import lsst.daf.butler as dafButler
from lsst.cp.pipe.ptc import PhotonTransferCurveExtractTask, PhotonTransferCurveExtractConfig
from lsst.cp.pipe.ptc import PhotonTransferCurveSolveTask, PhotonTransferCurveSolveConfig
```

Ignore these warnings!

Below we get a list of the exposures that are part of run 13144.  We can see the flat pairs.\
Note that the first group of flat pairs incorporate the neutral density filter ND_OD0.5.\
Then we build a list of flat pairs to be used in the PTC.

```python jupyter={"outputs_hidden": true} tags=[]
REPO_DIR = '/repo/main'
butler = dafButler.Butler(REPO_DIR, collections=["LSSTCam/raw/all","LSSTCam/calib", "u/cslage/bps_13144M"])
 
exposureList = []
for record in butler.registry.queryDimensionRecords("exposure", where="exposure.science_program='13144'"):
    exposureList.append([record.id, record.observation_type, record.exposure_time, record.physical_filter])
exposureList = sorted(exposureList, key=lambda x: x[0])  
flatList = []
print("expId \t Type \t Exposure Time \t Filter\n")    
for items in exposureList:
    print(f"{items[0]} \t {items[1]} \t {items[2]} \t\t {items[3]}")
    if items[1] == 'flat':
        flatList.append(items)

```

Now we build a list of objects needed for the PTC.  This is pretty deep into the working of the butler, so don't worry if you don't understand all of this.  We're just getting a list of stuff needed by the PtcExtractTask.  The "cpPtcProc" data is the data from the CCDs that have already been through ISR (Instrument Signature Removal).

```python tags=[]
detector = 94 # Just choosing an arbitrary detector

# This let's you skip some of the flat pairs.  A value of 1 will do all
# of the pairs.  A value of 10 will only do every 10th pair, etc.
skipNPairs = 1 

expIds = []
expDict = {}
metadata = []
i = 0
while i < (len(flatList) - 2):
    expTime = flatList[i][2]
    nextExpTime = flatList[i + 1][2]
    if abs(expTime - nextExpTime) < 1E-4:
        expId = flatList[i][0]
        nextExpId = flatList[i+1][0]
        #print(expId, nextExpId)
        expIds.append(expId)
        expIds.append(nextExpId)
        dataId = {'exposure':expId, 'detector':detector, 'instrument':'LSSTCam'}
        nextDataId = {'exposure':nextExpId, 'detector':detector, 'instrument':'LSSTCam'}
        ref1 = butler.getDeferred(datasetRefOrType='cpPtcProc', dataId=dataId)
        ref2 = butler.getDeferred(datasetRefOrType='cpPtcProc', dataId=nextDataId)
        expDict[str(expTime)] = ((ref1, expId), (ref2, nextExpId))
        meta1 = butler.get('isr_metadata', dataId=dataId)
        meta2 = butler.get('isr_metadata', dataId=nextDataId)
        metadata.append(meta1)
        metadata.append(meta2)
        i += 2 * skipNPairs
    else:
        i += 1
        
    
print(f"There are {len(expDict)} pairs in the list.")    
```

```python
PTCExtractConfig = PhotonTransferCurveExtractConfig()
```

We can examine the parameters in the PTCExtractConfig object as follows:
This will list parameters that we might want to adjust\
In this case we are happy with the defaults so we don't modify any of the parameters.

```python jupyter={"outputs_hidden": true} tags=[]
PTCExtractConfig?
```

Now we create the object for the PTC extract task

```python
PTCExtractTask = PhotonTransferCurveExtractTask(config=PTCExtractConfig)
```

Again, we can examine the parameters in the PTCExtract object:
This will also explain the inputs and outputs.

```python jupyter={"outputs_hidden": true} tags=[]
PTCExtractTask?
```

Now we run the PTCExtract task.  This will take a while, especially if we use all of the flat pairs (skipNPairs = 1).  The reason it takes a long time is that it is calculating the variances and covariances across all of the CCD.  When I ran all 337 flat pairs on one CCD, it took about 1-1/2 hours.

```python jupyter={"outputs_hidden": true} tags=[]
start = time.time()
PTCExtractResult = PTCExtractTask.run(inputExp=expDict, inputDims=expIds,
                                         taskMetadata=metadata)
finish = time.time()
print(f"Took {finish - start} seconds.")
```

Now we create and run the PTCSolve task.  This runs much faster.

```python
PTCSolveConfig = PhotonTransferCurveSolveConfig()
PTCSolveConfig.ptcFitType = "EXPAPPROXIMATION"
```

```python
PTCSolveTask = PhotonTransferCurveSolveTask()
```

```python jupyter={"outputs_hidden": true} tags=[]
PTCSolveResult = PTCSolveTask.run(PTCExtractResult.outputCovariances)
```

Now we can extract the result.  This will have objects like gain and noise that you are familiar with.

```python
ptcDataset = PTCSolveResult.outputPtcDataset
```

```python tags=[]
ptcDataset.gain
```

Now we can plot the result.

```python

```
