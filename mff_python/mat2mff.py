# -*- coding: utf-8 -*-
"""
Created on Mon Mar 16 13:55:51 2015
@author: jesong1126
"""

#import os 
import scipy.io as sio
import numpy as np   
import matplotlib.pyplot as plt 

TempData = sio.loadmat('./data/grandave.mat') 
MatData = TempData['data']
MatData = MatData.astype('float32')
MatData.dtype

plt.plot(MatData[:,:,0].T)
plt.show()

nC = 129
nE = nC - 1
nSamplesPre = 50
nSamples = 425
srate = 250
nTrials = 53 
baseline = (nSamplesPre * 1000 / srate)  
msSamples = np.arange(0, nSamples,1) * 1000/srate  - baseline 

MatData0 = np.zeros((nC, nSamples, nTrials))
for i in range(nTrials):
    a = np.concatenate((MatData[:,:,i], np.zeros((1, nSamples)))) 
    MatData0[:,:,i] = a

MatData = MatData0.astype('float32')
MatData.dtype
#
#plt.plot(MatData[:,:,0].T)
#plt.show()

#MatDataT = np.zeros((nE, nSamples, nTrials))
#for j in range(nTrials):
#    scalp = MatData[:,:,j]
#    scalpflat = scalp.flatten()   
#    scalp1 = scalpflat.reshape((nSamples, nE))
#    scalp0 = scalp1.T
#    MatDataT[:,:,j] = scalp0
#
#MatData0 = np.zeros((nC, nSamples, nTrials))
#for i in range(nTrials):
#    a = np.concatenate((MatDataT[:,:,i], np.zeros((1, nSamples)))) 
#    MatData0[:,:,i] = a

#    
#MatData0 = MatData0.astype('float32')
#MatData0.dtype
#plt.plot(MatDataT[:,:,0])
#plt.show()

#MatDataT = np.zeros((nC, nSamples, nTrials))
#for j in range(nTrials):
#    scalp = MatData[:,:,j]
#    scalpflat = scalp.flatten()   
#    scalp1 = scalpflat.reshape((nSamples, nC))
#    scalp0 = scalp1.T
#    MatDataT[:,:,j] = scalp0
##    
#MatDataT= MatDataT.astype('float32')
#MatDataT.dtype
#plt.plot(MatDataT[:,:,0])
#plt.show()


from mff_python import libMFF 
filePath ='/Users/jesong1126/Python27/GeoPy/mff_python/test_data/granderp.mff'

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
isOpen = ww.openFile('/Users/jesong1126/Python27/GeoPy/mff_python/test_data/MatData.mff')

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
rr.readSignalBlocksFromInfoN('/Users/jesong1126/Python27/GeoPy/mff_python/test_data/granderp.mff/signal1.bin', infoN, 0, -1)
#ww.writeMetaData(ww.CoordinateLayout_HydroCel256(), ww.SensorLayout_HydroCel256())
ww.writeMetaData(ww.CoordinateLayout_HydroCel128(), ww.SensorLayout_HydroCel128())
numblocks = rr.numSignalBlocks(infoN)


for bid in range(numblocks):
    block = rr.getSignalBlock(infoN, 0) 
    ar = MatData0[:,:,bid].flatten().reshape((nSamples, nC)).T
    #ar = np.kron(np.ones((1, nSamples)), np.arange(nC).reshape((nC,1)))
    ar = ar.astype('float32')
    #ar0 = ar.flatten().reshape((nSamples, nC)).T
    #block.setDataFromNumpyArray(MatData0[:,:,bid]) 
    block.setDataFromNumpyArray(ar) 
    ww.writeSignalBlockToFile(infoN, block)

#clean up after ourselves.
rr.releaseSignalBlocks()
ww.closeFile()

#    ar = block.getDataAsNumpyArray()
#    print ar
#    print ar.shape    
#    block.setDataFromNumpyArray(ar) 
