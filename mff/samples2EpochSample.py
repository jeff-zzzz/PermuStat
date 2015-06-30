## samples2EpochSample.py ----------------------------------------------------------------------
#  Python File
#  author Jasmine Song 
#  date 8/20/2014
#  Copyright 2014 EGI. All rights reserved.
##------------------------------------------------------------------------------------------  
# Converts from samples since start of recording (as if there were no
# breaks) to samples in file.   
def samples2EpochSample(sampleNum, epochBeginSamps, epochNumSamps):
    numEpochs = len(epochBeginSamps) 
    epochSampleNum = 0 
    p = 0 
    while (sampleNum > (epochBeginSamps[p] + epochNumSamps[p] - 1)) and (p < numEpochs):
        epochSampleNum = epochSampleNum + epochNumSamps[p] 
        p = p+1 

    if p <= numEpochs:
        if sampleNum >= epochBeginSamps[p]: 
            epochSampleNum = epochSampleNum + ((sampleNum - epochBeginSamps[p])+1) 
        else:
            epochSampleNum = -1
            # Error: sample falls between epochs
    else:
        epochSampleNum = -2 
        # Error: sample is after last epoch

    return epochSampleNum

## --------------------------------------------------------------------------------------
