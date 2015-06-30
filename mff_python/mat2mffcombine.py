# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 13:55:51 2015
@author: jesong1126
"""

#import os 
import scipy.io as sio
import numpy as np   
import matplotlib.pyplot as plt 

TempData = sio.loadmat('./data/MEP108RT.mat') 
RTData = TempData['data']
RTData = RTData.astype('float32')
RTData.dtype
nX1 = RTData.shape[2]

TempData = sio.loadmat('./data/MEP108LT.mat') 
LTData = TempData['data']
LTData = LTData.astype('float32')
LTData.dtype
nX2 = LTData.shape[2]


nC = 257
nE = nC - 1
nSamplesPre = 200
nSamples = 300
srate = 1000
nTrials = nX1 + nX2
baseline = (nSamplesPre * 1000 / srate)  
msSamples = np.arange(0, nSamples,1) * 1000/srate  - baseline 

#fig, ax= plt.subplots(2,2)
#ax[0,0].plot(msSamples, LTData[:,:,0].T)
#ax[0,1].plot(msSamples, RTData[:,:,0].T)
#plt.show()
# 
#RL = np.concatenate((RT227, LT259), axis=2) 
RL = np.concatenate((RTData, LTData), axis=2) 
RL = RL.astype('float32')
RL.dtype
 
from mff_python import libMFF 
filePath ='/Users/jesong1126/Python27/GeoPy/mff_python/test_data/empty_seg_486.mff'

rr = libMFF.MFFReader()
rr.openFile(filePath)  
numInfoN = rr.readInfoN()
print "Num InfoN: " + str(numInfoN)
numSignalBlock = rr.prepareSignalBlocks()
numCat = rr.readCategories()
numEpochs = rr.readEpochs()
numData = rr.readInfoN()
numPNS = rr.readPNSSets()
numEventTracks = rr.readEventTracks()

ww = libMFF.MFFWriter()
isOpen = ww.openFile('/Users/jesong1126/Python27/GeoPy/mff_python/test_data/SEP_107_0046_seg_486.mff')

# First we create an Info object.
info = libMFF.PyInfo()
info.setMFFVersion(3)
info.setRecordTime("2015-03-13T00:00:00.000000-07:00")
info.setAmpSerialNumber("A08080117")
info.setAmpFirmwareVersion("0")
info.setMovieDeltaPresent(False)
ww.setInfo(info)

# Copy the categories.
for i in range(numCat):
    ww.addCategory(rr.getCategory(i))

# Copy the Epochs
for i in range(numEpochs):
    ww.addEpoch(rr.getEpoch(i))

# Copy the InfoN data.
for i in range(numInfoN):
    print "adding infoN: " + str(i) 
    ww.addInfoN(rr.getInfoN(i))

# Copy the EventTracks
for i in range(numEventTracks):
    ww.addEventTrack(rr.getEventTrack(i))

# Write the meta data we copied now.
infoN=rr.getInfoN(0)
rr.readSignalBlocksFromInfoN('/Users/jesong1126/Python27/GeoPy/mff_python/test_data/empty_seg_486.mff/signal1.bin', infoN, 0, -1)
ww.writeMetaData(ww.CoordinateLayout_HydroCel256(), ww.SensorLayout_HydroCel256()) 
#ww.writeMetaData(ww.CoordinateLayout_HydroCel128(), ww.SensorLayout_HydroCel128())
numblocks = rr.numSignalBlocks(infoN)
#
#for bid in range(numblocks):
#    block = rr.getSignalBlock(infoN, 0) 
#    ar = block.getDataAsNumpyArray() #(257 300)
#    block.setDataFromNumpyArray(ar)
#    ww.writeSignalBlockToFile(infoN, block)

for bid in range(numblocks):
    block = rr.getSignalBlock(infoN, 0) 
    ar = RL[:,:,bid].flatten().reshape((nSamples, nC)).T 
    print ar.shape
    ar = ar.astype('float32')
    block.setDataFromNumpyArray(ar)
    ww.writeSignalBlockToFile(infoN, block)

#clean up after ourselves.
rr.releaseSignalBlocks()
ww.closeFile()

#    ar = block.getDataAsNumpyArray()
#    print ar
#    print ar.shape    
#    block.setDataFromNumpyArray(ar)
