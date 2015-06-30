#!/usr/bin/python
# cd ~/Python27/GeoPy/libMFF/build_anaconda

#import os 
#os.chdir("./libMFF/mff_python")
import numpy as np
#from braink import read_lfm, read_bk   
#import matplotlib.pyplot as plt 

HMdir = '/Users/jesong1126/Python27/GeoPy/data/AMale_40_256'
lfmfilename = HMdir+'/AMale_40_256_trip.lfm'
fd = open(lfmfilename, 'r')
Ktri = np.fromfile(file=fd, dtype=np.dtype('d')) #big endian ? 
fd.close()

nE = 256 
nC = 257
nV3 = Ktri.shape[0]/nE  #13257   
nV = nV3/3 #4419 

Ktri = Ktri.reshape(nE, nV3)

#H = np.eye((nC))- np.ones((nC, nC))/nC
Ktri2 = np.concatenate((Ktri, np.zeros((1, nV3)))) 
K = Ktri2 
#K = np.dot(H, Ktri2) 

Dipole3 = np.eye(nV3)
ii = 0 

from RabbitMQ import Connection
c = Connection.Connection('localhost')
c.connect()
#c.sendOrientedData(abs(K[:, ii]) , "Oriented", "3:00pm")
c.sendTriplesData(Dipole3[:,ii] , "Triples", "3:00pm")
c.sendEEGData(Ktri[:, ii], "Electrodes", "3:00pm")
c.disconnect()       

##
Dipole3141 = [41, 89, 187]
Dipole3151 = [127, 89, 187]

mepL = K[:, 3141*3] ; mepL =  mepL.reshape((nC, 1))
mepR = K[:, 3151*3] ; mepR =  mepR.reshape((nC, 1))

freqL = 10.0
freqR = 20.0

nSamples = 1600 
tp = np.arange(0, nSamples)/1000.0

signalL = np.sin(2 * np.pi *freqL * tp)
signalR = np.sin(2 * np.pi *freqR * tp)
signalL =  signalL.reshape((1, nSamples))
signalR =  signalR.reshape((1, nSamples))

scalpL = np.kron(mepL, signalL) 
scalpR = np.kron(mepR, signalR) 
scalpL = scalpL.astype('float32')
scalpR = scalpR.astype('float32')
 
# 1% noise  
nTrialsL = 33
Left3141 = np.zeros((nC, nSamples, nTrialsL))
for j in range(nTrialsL):
    noise = np.random.normal(0, 0.01*np.std(mepL), (nC, nSamples)) 
    scalpNoise = scalpL + noise    
    scalpflat = scalpNoise.flatten()   
    scalp1 = scalpflat.reshape((nSamples, nC))
    scalp0 = scalp1.T
    Left3141[:,:,j] = scalp0

nTrialsR = 33 
Right3151 = np.zeros((nC, nSamples, nTrialsR))
for j in range(nTrialsR):
    noise = np.random.normal(0, 0.01*np.std(mepR), (nC, nSamples)) 
    scalpNoise = scalpR + noise          
    scalpflat = scalpNoise.flatten()   
    scalp1 = scalpflat.reshape((nSamples, nC))
    scalp0 = scalp1.T
    Right3151[:,:,j] = scalp0
 
LR = np.concatenate((Left3141, Right3151), axis=2) 
LR = LR.astype('float32')
LR.dtype

#type(ar) #numpy.narray
#ar.shape #(257, 1600)
#ar.dtype #dtype('float32')

from mff_python import libMFF 
filePath ='/Users/jesong1126/Python27/GeoPy/mff_python/test_data/Empty_1600_66.mff'
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
rr.readSignalBlocksFromInfoN('/Users/jesong1126/Python27/GeoPy/mff_python/test_data/Empty_1600_66.mff/signal1.bin', infoN, 0, -1)
ww.writeMetaData(ww.CoordinateLayout_HydroCel256(), ww.SensorLayout_HydroCel256())
numblocks = rr.numSignalBlocks(infoN)
#import numpy as np
#LR0 = np.zeros((257, 1600)) 
#LR0 = LR0.astype('float32')
#LR0.dtype
for bid in range(numblocks):
    block = rr.getSignalBlock(infoN, 0) 
    ar = block.getDataAsNumpyArray()
#    block.setDataFromNumpyArray(ar)   #    
    block.setDataFromNumpyArray(LR[:,:,bid]) 
#    block.setDataFromNumpyArray(LR0) 
    ww.writeSignalBlockToFile(infoN, block)

#clean up after ourselves.
rr.releaseSignalBlocks()
ww.closeFile()

#    ar = block.getDataAsNumpyArray()
#    ar.shape
#    ar.dtype
#    scalpL.shape
#    scalpL.dtype        
#    block.setDataFromNumpyArray(ar)  
#    block.setDataFromNumpyArray(scalp0) 
#    block.setDataFromNumpyArray(scalp2) 

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
#
#block = rr.getSignalBlock(infoN, 0)
#block.setDataFromNumpyArray(eegMean0)
#ww.writeSignalBlockToFile(infoN, block)
#
#block = rr.getSignalBlock(infoN, 0)
#block.setDataFromNumpyArray(eegStd0)
#ww.writeSignalBlockToFile(infoN, block)
#
#block = rr.getSignalBlock(infoN, 0)
#block.setDataFromNumpyArray(eegT0)
#ww.writeSignalBlockToFile(infoN, block)
#
#block = rr.getSignalBlock(infoN, 0)
#block.setDataFromNumpyArray(eee0)
#ww.writeSignalBlockToFile(infoN, block)
#
#type(ar) #numpy.narray
#ar.shape #(257, 1600)
#ar.dtype #dtype('float32')

#nC = nE + 1
#H = np.eye((nC))- np.ones((nC, nC))/nC
#Ktri2 = np.concatenate((Ktri, np.zeros((1, nV)))) 
#K = np.dot(H, Ktri2) 
    
# 0-2399 0-1199 NumPoints=149332 233344
#if os.path.exists(HMdir+'/fdmForwardMatrixOriented'):          
#    K  = read_lfm.forward(HMdir+'/fdmForwardMatrixOriented', nE)  
#    print K 
#    print("fdmForwardMatrixOriented is loaded. ")            
#elif os.path.exists(HMdir+'/Leadfield.lfm'): 
#    K  = read_lfm.lfm(HMdir+'/Leadfield.lfm',  nE)                          
#    print("Leadfield.lfm is loaded. ")      
#    print K 
#else:
#    print("LFM is not loaded. ")

#nV = K.shape[1] # = 2383

##
#lbkd = read_bk.bkd('/Users/jesong1126/Python27/GeoPy/data/108_MIE/108_3.Dipoles_left_MRI.bkd')
#rbkd = read_bk.bkd('/Users/jesong1126/Python27/GeoPy/data/108_MIE/108_3.Dipoles_right_MRI.bkd')
##lbkd['ndipoles'] + rbkd['ndipoles']
#
#lcoord = np.array((lbkd['dipole_location_index']['x'], lbkd['dipole_location_index']['y'], lbkd['dipole_location_index']['z'])).T
#rcoord = np.array((rbkd['dipole_location_index']['x'], rbkd['dipole_location_index']['y'], rbkd['dipole_location_index']['z'])).T
#VoxCoord = np.concatenate((lcoord, rcoord), axis=0)
#nVL = lbkd['ndipoles'] 
#nVR = rbkd['ndipoles'] 

#VoxCoord[[1034, 2215],:]

#os.system("open /Applications/EAV/EAV.app")
#sKtri = np.cumsum(Ktri)
#nV3 = np.where(sKtri!=0)[0][0] 
#nE = Ktri.shape[0]/nV3 

