
##---------------------------------------------------------------------------------------
def mff_getSignalFilename(filePath, infoNType):
    import glob 
    import numpy as np
    from xml.dom.minidom import parse
    
    binfiles = []  
    for ff in glob.glob(filePath+'/signal*.bin'):
        binfiles.append(ff[len(filePath)+1:])

    signalFile = []
    infoFile = []
    if len(binfiles) > 1:
        for p in range(len(binfiles)): # p is from 1 to * to find PIBfile
            binFilename = binfiles[p]
            prefix = 'signal' 
            prefixLen = len(prefix)
            extension = '.bin'
            extensionLen = len(extension)
            binNumDigits = len(binFilename) - (prefixLen + extensionLen)
            binNumStr = binFilename[np.asarray(range(prefixLen,  prefixLen + binNumDigits))]

            infoObjFile = filePath+'/info'+binNumStr+'.xml'
            infoObj = parse(infoObjFile) 
            if infoObj.getElementsByTagName(infoNType) != None :
                infoFile = 'info'+binNumStr+'.xml'
                signalFile = 'signal'+binNumStr+'.bin'

    return signalFile

##--------------------------------------------------------------------------------------- 
## infoNType ='EEG' 'PNSData' 'Spectral' 'sourceData' 'JTF' 'TValues' 'Filter'  
 
