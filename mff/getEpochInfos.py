
def getEpochInfos(filePath, sampRate):
    import numpy as np
    from xml.dom.minidom import parse
    import os.path
    from mff.mff_micros2Sample import mff_micros2Sample

    ##------------------------------------------------------------------------------
    epochfile = filePath+'/epochs.xml'
    epochList = parse(epochfile)    
    epochs = epochList.getElementsByTagName('epoch')
    numEpochs = epochs.length 

    epochBeginSamps = np.zeros((numEpochs), dtype='i8')
    epochNumSamps = np.zeros((numEpochs), dtype='i8')
    epochFirstBlocks = np.zeros((numEpochs), dtype='i8')
    epochLastBlocks = np.zeros((numEpochs), dtype='i8')
    epochTime0 = np.zeros((numEpochs), dtype='i8')
    epochLabels = [None]*numEpochs #np.zeros((numEpochs), dtype='S50')
    epochSubjects = [] 
    epochFilenames = [] 
    epochSegStatus = [None]*numEpochs #np.zeros((numEpochs), dtype='S50')
    multiSubj = False 

    for p in range(numEpochs):
        anEpoch = epochs[p]
        epochBegin = int(anEpoch.getElementsByTagName('beginTime')[0].firstChild.data) 
        epochEnd = int(anEpoch.getElementsByTagName('endTime')[0].firstChild.data)
        epochBeginSamps[p] = mff_micros2Sample(epochBegin, sampRate)[0]
        epochTime0[p] = epochBeginSamps[p]
        epochNumSamps[p] = mff_micros2Sample(epochEnd, sampRate)[0]- epochBeginSamps[p]    
        epochFirstBlocks[p] = int(anEpoch.getElementsByTagName('firstBlock')[0].firstChild.data)
        epochLastBlocks[p] = int(anEpoch.getElementsByTagName('lastBlock')[0].firstChild.data)
        epochLabels[p] = 'epoch' 

    epochType = 'cnt'
    totalNumSegs = 0

    ##------------------------------------------------------------------------------------
    categfile = filePath+'/categories.xml' 
    if os.path.isfile(categfile):
        epochType = 'seg'
        categList = parse(categfile)    
        cats = categList.getElementsByTagName('cat')
        numCategs = cats.length 

        multiSubj = False
        if os.path.exists(filePath+'/subjects'):
            multiSubj = True
            for p in range(numCategs):
                aCateg = cats[p]
                segList = aCateg.getElementsByTagName('segments')
                numSegs = segList.length
                for q in range(numSegs):
                    aSeg = segList[q]
                    keyListArray = aSeg.getElementsByTagName('keyCode')
                    numKeys = keyListArray.length
                    numSubjs = 0
                    for r in range(numKeys):
                        aKey = keyListArray[r]
                        if aKey.firstChild.data.encode() == 'subj':
                            numSubjs = numSubjs + 1

        if multiSubj:            
            epochSubjects = [None]*numEpochs #np.zeros((numEpochs), dtype='S50')
            epochFilenames = [None]*numEpochs #np.zeros((numEpochs), dtype='S50')          
                    
        for p in range(numCategs):
            aCateg = cats[p]
            categLabel = aCateg.getElementsByTagName('name')[0].firstChild.data #.encode() 
            segList = aCateg.getElementsByTagName('seg')
            numSegs = segList.length
            totalNumSegs = totalNumSegs + numSegs
            for q in range(numSegs):
                aSeg = segList[q]
                segBegin = int(aSeg.getElementsByTagName('beginTime')[0].firstChild.data)            
                segBeginSamp = int(mff_micros2Sample(segBegin, sampRate)[0])
                segInd = np.where(epochBeginSamps == segBeginSamp)[0]
                epochSegStatus[segInd[0]] = aSeg.getAttribute('status') 
                epochLabels[segInd[0]] = categLabel
                time0 = int(aSeg.getElementsByTagName('evtBegin')[0].firstChild.data)
                time0Samp = mff_micros2Sample(time0, sampRate)[0]
                time0Samp = time0Samp - segBeginSamp
                epochTime0[segInd[0]] = time0Samp
                if multiSubj:
                    keyListArray = aSeg.getElementsByTagName('keyCode')
                    numKeys = keyListArray.length                     
                    dataListArray = aSeg.getElementsByTagName('data')
                    for r in range(numKeys):
                        aKey = keyListArray[r]
                        aData = dataListArray[r] 
                        subject = None
                        filename = None
                        if aData.firstChild != None: 
                            if aKey.firstChild.data.encode() =='subj': 
                                subject = aData.firstChild.data #.encode()
                            elif aKey.firstChild.data.encode() == 'FILE': 
                                filename = aData.firstChild.data.encode() 
                        epochSubjects[segInd[0]] = subject
                        epochFilenames[segInd[0]] = filename
     
    if (multiSubj and len(set(epochFilenames))==1 and len(set(epochSubjects)) == 1):
        epochSubjects = [] 
        epochFilenames = []
        multiSubj = False

    ##------------------------------------------------------------------------------------

    epochInfo = {'epochType':epochType, 'epochBeginSamps':epochBeginSamps, 
    'epochNumSamps':epochNumSamps, 'epochFirstBlocks':epochFirstBlocks, 
    'epochLastBlocks':epochLastBlocks, 'epochLabels':epochLabels, 
    'epochTime0':epochTime0, 'multiSubj':multiSubj, 'epochSubjects':epochSubjects, 
    'epochFilenames':epochFilenames, 'epochSegStatus':epochSegStatus}

    
    return epochInfo

 
