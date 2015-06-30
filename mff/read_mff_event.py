## read_mff_event.py ----------------------------------------------------------------------
#  Python File
#  author Jasmine Song 
#  date 8/20/2014
#  Copyright 2014 EGI. All rights reserved.
##------------------------------------------------------------------------------------------
def read_mff_event(filePath, hdr): 
    import numpy as np
    from xml.dom.minidom import parse
    import os.path
    from mff.mff_micros2Sample import mff_micros2Sample
    import glob
    import os
    from mff.mff_getSummaryInfo import mff_getSummaryInfo
    from mff.ns2pyTime import ns2pyTime
    from mff.mff_micros2Sample import mff_micros2Sample
    from mff.samples2EpochSample import samples2EpochSample
    
    if hdr == None :
        summaryInfo = mff_getSummaryInfo(filePath)
    else:
        summaryInfo = hdr['orig']

    # Pull the information about the epochs out of the summary info
    epochBeginSamps = summaryInfo['epochBeginSamps'] 
    epochNumSamps = summaryInfo['epochNumSamps'] 
    epochFirstBlocks = summaryInfo['epochFirstBlocks'] 
    epochLastBlocks = summaryInfo['epochLastBlocks']  
    
    infofile = filePath+'/info.xml'
    infotime = parse(infofile)    
    beginTime = infotime.getElementsByTagName('recordTime')[0].firstChild.data.encode()
    bTime = ns2pyTime(beginTime)

    numEpochs = len(epochBeginSamps) 
    events = []
    eventInds = np.zeros((numEpochs,2), dtype='i')  
    eventInd = 0
    
    for p in range(numEpochs):
        origp = [0, 0, 'metadata',{}]
        eventp = ['break '+ summaryInfo['epochType'],
                  summaryInfo['epochBeginSamps'][p],
                  summaryInfo['epochLabels'][p], [],
                  summaryInfo['epochNumSamps'][p], [], origp ]
        events.append(eventp)
        eventInds[eventInd,:] = [eventp[1], eventInd]
        eventInd = eventInd+1 
           
    for file in glob.glob(filePath+"/Events_*.xml"):
        eventsfile = file

    eventsListObj = parse(eventsfile)
    trackname = eventsListObj.getElementsByTagName('name')[0].firstChild.data.encode() 
    eventsList = eventsListObj.getElementsByTagName('event')
    numEvents = eventsList.length
    
    for p in range(numEvents):
        theEvent= eventsList[p]
        eventTime = theEvent.getElementsByTagName('beginTime')[0].firstChild.data.encode()
        eTime = ns2pyTime(eventTime)
        eventTimeInMicros = (eTime-bTime).microseconds
        out = mff_micros2Sample(eventTimeInMicros,summaryInfo['sampRate'])
        eventTimeInSamples = int(out[0])
        sampleRemainder = int(out[1])
        eventTimeInEpochSamples= samples2EpochSample(eventTimeInSamples, epochBeginSamps, epochNumSamps)

        if eventTimeInEpochSamples >= 0:
            keylist = theEvent.getElementsByTagName('key') 
            eventkeycount = len(keylist) 
            keyp = []
            for q in range(eventkeycount):
                theKey = keylist[q] 
                theKeyCode = theKey.getElementsByTagName('keyCode')[0].firstChild.data.encode()
                theKeyData = theKey.getElementsByTagName('data')[0].firstChild.data.encode()
                theKeyData1 = theKey.getElementsByTagName('data')
                theKeyData0 = theKeyData1[0] # print theKeyData0.toxml()
                a = theKeyData0.attributes["dataType"]
                theKeyDatatype = a.value.encode()
                theKeyq = [theKeyCode, int(theKeyData), theKeyDatatype,'']
                keyp.append(theKeyq)   
            
            duration = int(theEvent.getElementsByTagName('duration')[0].firstChild.data.encode())
            out = mff_micros2Sample(duration, summaryInfo['sampRate'])   
                
            origp =[sampleRemainder, int(out[1]), trackname, keyp]
            eventp=[theEvent.getElementsByTagName('code')[0].firstChild.data.encode(),eventTimeInEpochSamples,0, [], int(out[0]), [], origp]
            
            events.append(eventp)
            eventIndp = np.array([eventp[1], eventInd], dtype='i').reshape((1,2))
            eventInds = np.concatenate((eventInds,eventIndp), axis=0)
            eventInd = eventInd+1 

    eventInds0 = eventInds[eventInds[:,0].argsort()]
    sortedEvents = events 
    for p in range(eventInd): 
        nextEventInd = eventInds0[p,1]
        sortedEvents[p] = events[nextEventInd] 
    events = sortedEvents
                                
    return events

 
## ------------------------------------------------------------------------------------
#        origp = {'sampleRemainder':0,'durationRemainder':0,'trackname':'metadata','keys':{}}
#        eventp = {'type':'break '+ summaryInfo['epochType'],
#        'sample':summaryInfo['epochBeginSamps'][p],
#        'value':summaryInfo['epochLabels'][p],
#        'offset':[],
#        'duration':summaryInfo['epochNumSamps'][p],
#        'timestamp':[],
#        'orig': origp }
#        events = {'type','sample','value','offset','duration','timestamp','orig'}
#        eventInds.append(eventIndp)
#        eventInds ={'0','1'}
## --------------------------------------------------------------------------------------







