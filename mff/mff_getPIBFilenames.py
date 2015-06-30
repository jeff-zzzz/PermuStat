  
  
##---------------------------------------------------------------------------------------
def mff_getPIBFilenames(filePath):
    import glob 
    import numpy as np
    from xml.dom.minidom import parse
    
    binfiles = []  
    for ff in glob.glob(filePath+'/signal*.bin'):
        binfiles.append(ff[len(filePath)+1:])

    signalFile = []
    infoFile = []
    if len(binfiles) > 1:
        for p in range(1, len(binfiles)): # p is from 1 to * to find PIBfile
            binFilename = binfiles[p]
            # All this to strip the number (binNumStr) from the signal file in order to apply
            # it to the info file. 
            prefix = 'signal' 
            prefixLen = len(prefix)
            extension = '.bin'
            extensionLen = len(extension)
            binNumDigits = len(binFilename) - (prefixLen + extensionLen)
            binNumStr = binFilename[np.asarray(range(prefixLen,  prefixLen + binNumDigits))]

            infoObjFile = filePath+'/info'+binNumStr+'.xml'
            infoObj = parse(infoObjFile) 
            # 'PNSData' 'Spectral' 'sourceData' 'JTF' 'TValues' 'Filter'
            if infoObj.getElementsByTagName('PNSData') != None :
                infoFile = 'info'+binNumStr+'.xml'
                signalFile = 'signal'+binNumStr+'.bin'

    pibInfo = {'pibSignalFile':signalFile, 'pibInfoFile': infoFile}
    return pibInfo

##--------------------------------------------------------------------------------------- 
