
##---------------------------------------------------------------------------------------
# Returns the bin objects corresponding to the EEG, and PIB if it exists,
# and the blocks arrays associated associated with them.

def getSignalBlocks(filePath): 
    from mff.mff_getSignalFilename import mff_getSignalFilename            
    from mff.getSignalNbin import getSignalNbin            

#    eegSignalFile = mff_getSignalFilename(filePath,'EEG')
    eegSignalFile = 'signal1.bin' 
    eegInfo = getSignalNbin(filePath, eegSignalFile)
    pibSignalFile = mff_getSignalFilename(filePath, 'PNSData')
    if pibSignalFile == []:
        pibInfo ={'nC':0,'sampRate':0,'blocks':0,'blockNumSamps':[],'signalFile':[]}
    else:
        pibInfo = getSignalNbin(filePath, pibSignalFile)

    SignalBlocks = {'npibChan':pibInfo['nC'],'pibBinObj':pibInfo['blockNumSamps'], 'pibBlocks':pibInfo['blocks'], 
    'pibSignalFile':pibSignalFile,'sampRate':eegInfo['sampRate'], 'nChan':eegInfo['nC'],'binObj':eegInfo['blockNumSamps'],'blocks':eegInfo['blocks'],'eegFile':eegSignalFile}

    return SignalBlocks

   
