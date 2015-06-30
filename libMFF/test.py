#!/usr/bin/python
# cd ~/Python27/GeoPy/libMFF/build_anaconda

import os 
os.chdir("./libMFF/build_anaconda")

import libMFF
#import numpy as np
#import matplotlib.pyplot as plt

#eegMean0 = np.concatenate((eegMean, np.zeros((1, nSamples))), axis=0)
#eegStd0 = np.concatenate((eegStd, np.zeros((1, nSamples))), axis=0)
#eegT0 = np.concatenate((eegT, np.zeros((1, nSamples))), axis=0) 
#eee0 = np.concatenate((eee, np.zeros((1, nSamples))), axis=0)
#
#eegMean =  eegMean.astype(np.float32)
#eegStd=  eegStd.astype(np.float32)
#eegT =  eegT.astype(np.float32)
#eee =  eee.astype(np.float32)
#
#eegMean0 =  eegMean0.astype(np.float32)
#eegStd0 =  eegStd0.astype(np.float32)
#eegT0 =  eegT0.astype(np.float32)
#eee0 =  eee0.astype(np.float32)

filePath ='/Users/jesong1126/Python27/GeoPy/libMFF/Empty.mff'

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
isOpen = ww.openFile('/Users/jesong1126/Python27/GeoPy/mff_python/test_data/Wavelet.mff')

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
infoN = rr.getInfoN(0)
rr.readSignalBlocksFromInfoN('/Users/jesong1126/Python27/GeoPy/libMFF/Empty.mff/signal1.bin', infoN, 0, -1)
#r = ww.writeMetaData(ww.CoordinateLayout_HydroCel256(), ww.SensorLayout_HydroCel256())
numblocks = rr.numSignalBlocks(infoN)

block = rr.getSignalBlock(infoN, 0)
block.setDataFromNumpyArray(eegMean0)
ww.writeSignalBlockToFile(infoN, block)

block = rr.getSignalBlock(infoN, 0)
block.setDataFromNumpyArray(eegStd0)
ww.writeSignalBlockToFile(infoN, block)

block = rr.getSignalBlock(infoN, 0)
block.setDataFromNumpyArray(eegT0)
ww.writeSignalBlockToFile(infoN, block)

block = rr.getSignalBlock(infoN, 0)
block.setDataFromNumpyArray(eee0)
ww.writeSignalBlockToFile(infoN, block)

#for bid in range(numblocks):
#    block = rr.getSignalBlock(infoN, 0)
#    ar = block.getDataAsNumpyArray()
#    block.setDataFromNumpyArray(ar)
#    ww.writeSignalBlockToFile(infoN, block)

#type(ar) #numpy.narray
#ar.shape #(257, 1600)
#ar.dtype #dtype('float32')

#clean up after ourselves.
rr.releaseSignalBlocks()
ww.closeFile()

## Copy the PNS Data
#if numPNS > 0: 
#    info.setPNSSet(rr.getPNSSet(i))

##-----------------
#filePath ='/Users/jesong1126/Python27/GeoPy/mff_python/test_data/testOutput.mff'
#rr = libMFF.MFFReader()
#rr.openFile(filePath)  
#numInfoN = rr.readInfoN()
#print "Num InfoN: " + str(numInfoN)
#numSignalBlock = rr.prepareSignalBlocks()
#numCat = rr.readCategories()
#numEpochs = rr.readEpochs()
#numData = rr.readInfoN()
#numPNS = rr.readPNSSets()
#numEventTracks = rr.readEventTracks()
#infoN = rr.getInfoN(0)
#rr.readSignalBlocksFromInfoN('/Users/jesong1126/Python27/GeoPy/mff_python/test_data/testOutput.mff/signal1.bin', infoN, 0, -1)
#numblocks = rr.numSignalBlocks(infoN)
#for bid in range(numblocks):
#    block = rr.getSignalBlock(infoN, 0)
#    ar = block.getDataAsNumpyArray()
#    print ar
#


