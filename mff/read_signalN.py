
def read_signalN(filePath, Nbinfile, blocksInEpoch):
    import numpy as np
 
    binfile = filePath+'/'+ Nbinfile    
    fid = open(binfile, 'rb')  
    s = fid.read()
    l = len(s)
    fid.seek(0,0)
    position = fid.tell()
    datalist = []
    blockI = 0 
    while position < l:  
        version = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
        if version == 0: 
            if blockI in blocksInEpoch: 
                data1 = np.fromfile(fid, dtype= np.dtype('f4'), count=hl)
                data1 = data1.reshape((nC, nSamples)) 
                datalist.append(data1)
            else: 
                position = fid.tell()
                fid.seek(position+blocksize)
        else:
            headersize = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
            blocksize = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
            hl = int(blocksize/4) 
            nC = np.fromfile(fid, dtype= np.dtype('i4'), count=1)[0]
            nSamples = int(hl/nC)
            sigOffset = np.fromfile(fid, dtype= np.dtype('i4'), count=nC)
            sigFreq = np.fromfile(fid, dtype= np.dtype('i4'), count=nC)
            srate = int((sigFreq[0]-32)/(nC-1))
            optionHead = np.fromfile(fid, dtype= np.dtype('i4'), count=int(headersize/4-(4+2*nC)))
            if blockI in blocksInEpoch: 
                data1 = np.fromfile(fid, dtype= np.dtype('f4'), count=hl)
                data1 = data1.reshape((nC, nSamples)) 
                datalist.append(data1)
            else:
                position = fid.tell()
                fid.seek(position+blocksize)
        position = fid.tell()
        blockI = blockI + 1 
    else:
        fid.close()

    return datalist 

 
