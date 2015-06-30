
def mff_getSummaryInfo(filePath):
    import numpy as np
    from xml.dom.minidom import parse
    from mff.getEpochInfos import getEpochInfos
    from mff.getSignalBlocks import getSignalBlocks

    SignalBlocks = getSignalBlocks(filePath)
    sampRate = SignalBlocks['sampRate']
    numblocks = SignalBlocks['blocks']
    blockNumSamps = np.array(SignalBlocks['binObj'])  

    ##------------------------------------------------------------------------------------
    pibHasRef = False
    pibNChans = 0 
    if SignalBlocks['pibSignalFile'] != []:
        pnsSetFile = filePath+'/pnsSet.xml'
        pnsSetObj = parse(pnsSetFile)
        pnsSensors = pnsSetObj.getElementsByTagName('sensor')
        pibNChans = pnsSensors.length
        if SignalBlocks['npibChan'] - pibNChans == 1 :
            pibHasRef = True        
 
    ##------------------------------------------------------------------------------------
    epochInfo = getEpochInfos(filePath, sampRate)

    ##------------------------------------------------------------------------------------
    blockBeginSamps = np.zeros((numblocks), dtype='i8')
    for x in range(0, (numblocks-1)):
        blockBeginSamps[x+1] = blockBeginSamps[x] + blockNumSamps[x]
        
    ##------------------------------------------------------------------------------------
    summaryInfo = {'blocks':SignalBlocks['blocks'],'eegFilename':SignalBlocks['eegFile'],
    'sampRate':SignalBlocks['sampRate'],'nChans':SignalBlocks['nChan'],
    'pibBinObj':SignalBlocks['pibBinObj'],'pibBlocks':SignalBlocks['pibBlocks'],
    'pibNChans':pibNChans,'pibFilename':SignalBlocks['pibSignalFile'],
    'pibHasRef':pibHasRef,'epochType':epochInfo['epochType'],
    'epochBeginSamps':epochInfo['epochBeginSamps'],'epochNumSamps':epochInfo['epochNumSamps'],
    'epochFirstBlocks':epochInfo['epochFirstBlocks'],'epochLastBlocks':epochInfo['epochLastBlocks'],
    'epochLabels':epochInfo['epochLabels'],'epochTime0':epochInfo['epochTime0'],
    'multiSubj':epochInfo['multiSubj'],'epochSubjects':epochInfo['epochSubjects'],
    'epochFilenames':epochInfo['epochFilenames'],'epochSegStatus':epochInfo['epochSegStatus'],
    'blockBeginSamps':blockBeginSamps,'blockNumSamps':blockNumSamps}
                
    return summaryInfo 

 
