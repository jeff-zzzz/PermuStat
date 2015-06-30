
def read_mff_data(filePath, indType, startInd, lastInd, hdr):
    import numpy as np
    from mff.mff_getSummaryInfo import mff_getSummaryInfo
    from mff.blockSample2BlockAndSample import blockSample2BlockAndSample
    from mff.read_signalN import read_signalN

    if hdr != None:
        summaryInfo = hdr['orig']
    else: 
        summaryInfo = mff_getSummaryInfo(filePath) 

    ##------------------------------------------------------------------------------
    blockNumSamps = summaryInfo['blockNumSamps']
    if indType == 'sample':
        [startBlockNum, startSamp] = blockSample2BlockAndSample(startInd-1, blockNumSamps) 
        [lastBlockNum, lastSamp] = blockSample2BlockAndSample(lastInd, blockNumSamps)
        blocksInEpoch = range(startBlockNum, lastBlockNum+1)
    elif indType == 'epoch':
        blocksInEpoch = [] 
        for i in range(startInd-1, lastInd): 
            blocksInd = range(summaryInfo['epochFirstBlocks'][i]-1,summaryInfo['epochLastBlocks'][i]) 
            blocksInEpoch.extend(blocksInd)
    else:
        print("Error: indType must to be either 'epoch' or 'sample'") 

    ##--------------------------------------------------------------------------------        
    dataeeg = read_signalN(filePath, 'signal1.bin', blocksInEpoch)
    if summaryInfo['pibFilename'] != []: 
        datapib = read_signalN(filePath, summaryInfo['pibFilename'], blocksInEpoch)
        datalist = [] 
        for i in range(len(dataeeg)):
            dataI = np.concatenate((dataeeg[i], datapib[i]), axis=0)
            datalist.append(dataI)   
    else:
        datalist = dataeeg 
    
    ##--------------------------------------------------------------------------------        
    if indType == 'sample': 
        if startBlockNum == lastBlockNum : 
            data = datalist[0][:, startSamp:lastSamp] 
        elif lastBlockNum-startBlockNum == 1: 
            datalistS = datalist[0][:,startSamp:] 
            datalistL = datalist[1][:,:lastSamp] 
            data = np.concatenate((datalistS, datalistL), axis=1) 
        else:
            datalistS = datalist[0][:,startSamp:]
            for j in range(1, len(datalist)-1):
                datalistJ = datalist[j] 
                datalistS = np.concatenate((datalistS, datalistJ), axis=1)
            datalistL = datalist[-1][:,:lastSamp]  
            data = np.concatenate((datalistS, datalistL), axis=1) 

    elif indType == 'epoch':    
        if len(datalist) > 1:
            if len(set(blockNumSamps)) > 1 : 
                data = datalist[0]
                for j in range(1,len(datalist)): 
                    datalistj = datalist[j]
                    data = np.concatenate((data, datalistj), axis=1) 
            else :
                data = np.zeros((datalist[0].shape[0], datalist[0].shape[1], len(datalist)))  
                for j in range(len(datalist)):
                    data[:,:,j] = datalist[j]
        else:
            data = datalist[0]
    else:
        print("Error: indType must to be either 'epoch' or 'sample'") 

    return data  

    ##-------------------------------------------------------------------------------- 
    ## binfile = filePath+'/signal1.bin' 
    ## fid = open(binfile, 'rb')  
    ## s = fid.read()
    ## l = len(s)
    ## fid.seek(0,0)
    ## position = fid.tell()
    ## datalist = []
    ## blockI = 0 
    ## while position < l:  
    ##     version = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
    ##     if version == 0: 
    ##         if blockI in blocksInEpoch: 
    ##             data1 = np.fromfile(fid, dtype= np.dtype('f4'), count=hl)
    ##             data1 = data1.reshape((nC, nSamples)) 
    ##             datalist.append(data1)
    ##         else: 
    ##             position = fid.tell()
    ##             fid.seek(position+blocksize)
    ##     else:
    ##         headersize = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
    ##         blocksize = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
    ##         hl = blocksize/4 
    ##         nC = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
    ##         nSamples = hl/nC
    ##         sigOffset = np.fromfile(fid, dtype= np.dtype('i4'), count=nC)
    ##         sigFreq = np.fromfile(fid, dtype= np.dtype('i4'), count=nC)
    ##         srate = (sigFreq[0]-32)/(nC-1)
    ##         optionHead = np.fromfile(fid, dtype= np.dtype('i4'), count=(headersize/4-(4+2*nC)))
    ##         if blockI in blocksInEpoch: 
    ##             data1 = np.fromfile(fid, dtype= np.dtype('f4'), count=hl)
    ##             data1 = data1.reshape((nC, nSamples)) 
    ##             datalist.append(data1)
    ##         else:
    ##             position = fid.tell()
    ##             fid.seek(position+blocksize)
    ##     position = fid.tell()
    ##     blockI = blockI + 1 
    ## else:
    ##     fid.close()
    ##
    ##-------------------------------------------------------------------------------- 

