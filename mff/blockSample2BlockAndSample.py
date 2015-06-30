
#sampleNum = 700000
#blockNumSamps = summaryInfo['blockNumSamps']
def blockSample2BlockAndSample(sampleNum, blockNumSamps):
    blockNum = 0
    numSamps = blockNumSamps[blockNum] #blockNumSamps(blockNum)
    while sampleNum > numSamps :
        blockNum = blockNum + 1
        numSamps = numSamps + blockNumSamps[blockNum] #blockNumSamps(blockNum)

    numSamps = numSamps - blockNumSamps[blockNum] #blockNumSamps(epochNum)    
    sample = sampleNum - numSamps             
    return [blockNum, sample]


 
