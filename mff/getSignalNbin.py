
def getSignalNbin(filePath, signalNbin): 
    import numpy as np

    binfile = filePath +'/' + signalNbin  
    fid = open(binfile, 'rb')  
    s = fid.read()
    l = len(s)
    fid.seek(0,0)
    position = fid.tell()
    blockNumSamps = []
    numblocks = 0
    while position < l:
        version = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
        numblocks = numblocks + 1
        if version == 0:
            blockNumSamps.append(nSamples)
            position = fid.tell()
            fid.seek(position+blocksize)
            position = fid.tell()
        else:
            headersize = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
            blocksize = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
            hl = int(blocksize/4) 
            nC = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
            nSamples = int(hl/nC)
            blockNumSamps.append(nSamples)
            sigOffset = np.fromfile(fid, dtype= np.dtype('i4'), count=nC)
            sigFreq = np.fromfile(fid, dtype= np.dtype('i4'), count=nC)
            sampRate = (sigFreq[0]-32)/(nC-1)
            optionHead = np.fromfile(fid, dtype= np.dtype('i4'), count=int(headersize/4-(4+2*nC)))
            position = fid.tell()
            fid.seek(position+blocksize)
            position = fid.tell() 
    else:
        fid.close()

    blockNumSamps = np.array(blockNumSamps)        
    SignalBlocks ={'nC':nC,'sampRate':sampRate,'blocks':numblocks,'blockNumSamps':blockNumSamps,'signalFile':signalNbin}

    return SignalBlocks



