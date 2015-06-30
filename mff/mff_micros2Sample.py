## mff_micros2Sample.py
#  Python File
#  author Jasmine Song
#  date 7/21/2014
#  Copyright 2014 EGI. All rights reserved.
#  Support routine for MFF Python code. Not intended to be called directly.
#
#  Converts from nanoseconds to samples, given the sampling rate. 
##

def mff_micros2Sample(microsecs, sampRate):
    import numpy as np
    
    sampDuration = 1000000/sampRate  
    sampleNum = microsecs/sampDuration 
    remainder = microsecs % sampDuration
    sampleNum = np.fix(sampleNum) 
#    out = {'sampleNum':sampleNum, 'remainder':remainder}
    out = [sampleNum, remainder]
    return out 

## function [sampleNum, remainder] = mff_micros2Sample(microsecs, sampRate)
## microsecs = double(microsecs);
## sampDuration = 1000000/sampRate;
## sampleNum = microsecs/sampDuration;
## remainder = uint64(rem(microsecs, sampDuration));
## sampleNum = fix(sampleNum);
