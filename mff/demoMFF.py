
import numpy as np
from xml.dom.minidom import parse
import os.path
import glob 
import matplotlib.pyplot as plt
from mff.read_mff_header import read_mff_header
from mff.read_mff_data import read_mff_data
from mff.read_mff_event import read_mff_event

from mff.mff_getSummaryInfo import mff_getSummaryInfo
from mff.getSignalBlocks import getSignalBlocks
from mff.getSignalNbin import getSignalNbin 
from mff.mff_getSignalFilename import mff_getSignalFilename
from mff.mff_micros2Sample import mff_micros2Sample
from mff.blockSample2BlockAndSample import blockSample2BlockAndSample
from mff.read_signalN import read_signalN

from mff.ns2pyTime import ns2pyTime

mfffilepath ='/Users/jesong1126/Python3/gs3/data' 
filePath = mfffilepath+'/SEP_107_0046_fil_seg_bcr.mff'

signalNbin = 'signal1.bin'
SignalBlocks = getSignalNbin(filePath, signalNbin) 

 
###----------------------------------------------------------------------------- 
# aaa = getSignalNbin(filePath, 'signal1.bin') 
# getEpochInfos(filePath, 500)
# SignalBlocks = getSignalBlocks(filePath)
# summaryInfo = mff_getSummaryInfo(filePath)  

# filePath = '/Users/jesong1126/Python27/Data/mffdata_test/MEP_107_LT.mff'
hdr = read_mff_header(filePath)
hdr.keys()
summaryInfo = hdr['orig']

#event1 = read_mff_event(filePath, hdr)

hdr['nTrials']
dataEpoch = read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr)
print(dataEpoch.shape)

ts = (np.arange(0,  hdr['nSamples']) - hdr['nSamplesPre']) * (1000 / hdr['Fs'])

plt.figure(1)
plt.plot(ts, dataEpoch[:,:,0].T, ':')
plt.xlim((ts[0], ts[-1]))


hdr['nSamples']
startInd = 4001; lastInd = 40000 ;
dataSample = read_mff_data(filePath, 'sample', startInd, lastInd, hdr)
print dataSample.shape

ts = (np.arange(0,  hdr['nSamples']) - hdr['nSamplesPre']) * (1000 / hdr['Fs'])

plt.figure(1)
#plt.plot(ts, data1.T, ':')
#plt.plot(dataSample.T, ':')
plt.plot(ts, dataSample[:,:,0].T, ':')
plt.xlim((ts[0], ts[-1]))

###-----------------------------------------------------------------------------
#    
#    if hdr == None :
#        summaryInfo = mff_getSummaryInfo(filePath)
#    else:
#        summaryInfo = hdr['orig']
#
#    # Pull the information about the epochs out of the summary info
#    epochBeginSamps = summaryInfo['epochBeginSamps'] 
#    epochNumSamps = summaryInfo['epochNumSamps'] 
#    epochFirstBlocks = summaryInfo['epochFirstBlocks'] 
#    epochLastBlocks = summaryInfo['epochLastBlocks']  
#    
#    infofile = filePath+'/info.xml'
#    infotime = parse(infofile)    
#    beginTime = infotime.getElementsByTagName('recordTime')[0].firstChild.data.encode()
#    bTime = ns2pyTime(beginTime)
#
#    numEpochs = len(epochBeginSamps) 
#    events = {'type':np.zeros((numEpochs), dtype='S50'),'sample':np.zeros((numEpochs), dtype='i4'), 'value':np.zeros((numEpochs), dtype='S50'), 'offset':np.zeros((numEpochs), dtype='i4'),'duration':np.zeros((numEpochs), dtype='i4'), 'timestamp':np.zeros((numEpochs), dtype='i4'),'sampleRemainder':np.zeros((numEpochs), dtype='i4'), 'durationRemainder':np.zeros((numEpochs), dtype='i4')}
#    eventInds ={'0':np.zeros((numEpochs), dtype='i4'), '1':np.zeros((numEpochs), dtype='i4') }  #dict() 
#    eventInd = 0
#    for p in range(numEpochs):
#        events['type'][eventInd] = 'break '+ summaryInfo['epochType']
#        events['sample'][eventInd] = summaryInfo['epochBeginSamps'][p]
#        events['value'][eventInd] = summaryInfo['epochLabels'][p]
#        events['offset'][eventInd] = 0
#        events['duration'][eventInd] = summaryInfo['epochNumSamps'][p] 
#        events['timestamp'][eventInd] = 0
#        events['sampleRemainder'][eventInd] = 0 
#        events['durationRemainder'][eventInd] = 0
#        eventInds['0'][eventInd] = events['sample'][eventInd] 
#        eventInds['1'][eventInd] = eventInd 
#        eventInd = eventInd+1 
#
#    # os.chdir(filePath)
#    for file in glob.glob(filePath+"/Events_*.xml"):
#        eventsfile = file
#
#    eventsListObj = parse(eventsfile)
#    eventsList = eventsListObj.getElementsByTagName('event')
#    numEvents= eventsList.length
#    for p in range(numEvents):
#        theEvent= eventsList[p]
#        eventTime = theEvent.getElementsByTagName('beginTime')[0].firstChild.data.encode()
#        eTime = ns2pyTime(eventTime)
#        eventTimeInMicros = (eTime-bTime).microseconds 
#        events['type'][eventInd] = 'break '+ summaryInfo['epochType']
#        events['sample'][eventInd] = summaryInfo['epochBeginSamps'][p]
#        events['value'][eventInd] = summaryInfo['epochLabels'][p]
#        events['offset'][eventInd] = 0
#        events['duration'][eventInd] = summaryInfo['epochNumSamps'][p] 
#        events['timestamp'][eventInd] = 0
#        events['sampleRemainder'][eventInd] = 0 
#        events['durationRemainder'][eventInd] = 0
#        eventInds['0'][eventInd] = events['sample'][eventInd] 
#        eventInds['1'][eventInd] = eventInd 
#        eventInd = eventInd+1 

##--------------------------------------------------------------------
           
## from datetime import datetime
## ## dt_str = '9/24/2010 5:03:29 PM'
## ## dt_obj = datetime.strptime(dt_str, '%m/%d/%Y %I:%M:%S %p')

## beginTime0 = '2013-07-10 12:21:48.000000'
## eventTime0 = '2013-07-10 12:21:48.067000'

## b0 = datetime.strptime(beginTime0, '%Y-%m-%d %H:%M:%S.%f')
## e0 = datetime.strptime(eventTime0, '%Y-%m-%d %H:%M:%S.%f')

## (e0-b0).microseconds

## beginDate = beginTime[0:10]
## beginTime0 = beginTime[12:26]
## beginTime00 = beginDate + " " +beginTime0
## #--------------------------------------------------------------------------        
 
## def ns2pyTime(nsTime):
##     from datetime import datetime
##     nsDate = nsTime[0:10]
##     nsTime0 = nsTime[11:26]
##     nsTime00 = nsDate + " " + nsTime0
##     pyTime = datetime.strptime(nsTime00, '%Y-%m-%d %H:%M:%S.%f')
##     return pyTime

## beginTime0 = ns2pyTime(beginTime)
## eventTime0 = ns2pyTime(eventTime)
## (eventTime0-beginTime0).microseconds 



