# -*- coding: utf-8 -*-
"""
Created on Fri Dec  5 11:33:11 2014
@author: jesong1126
"""

# cd '/Users/jesong1126/Python27/GeoPy'
from mff import read_mff_header, read_mff_data  #, getEpochInfos, mff_getSummaryInfo

filePath ='/Users/jesong1126/Python27/mff_anaconda/test_data/testInput.mff'
filePath ='/Users/jesong1126/Python27/mff_anaconda/test_data/MEP_108_gav_blc_bcr_ref.mff'
filePath ='/Users/jesong1126/Python27/mff_anaconda/test_data/SEP_108_0691_blc_ave_aref.mff'
filePath ='/Users/jesong1126/Python27/mff_anaconda/test_data/VGT_8subj_bcr_blc_ave.mff'


hdr = read_mff_header.read_mff_header(filePath)
        
nC = hdr['nChans']
nE = nC - 1
nSamples = hdr['nSamples']
nSamplesPre = hdr['nSamplesPre']
nTrials = hdr['nTrials']
srate = hdr['Fs']
summaryInfo = hdr['orig'] 
trialsName = summaryInfo['epochLabels']   
categoryName = list(set(trialsName))
nCategory = len(categoryName)

data = read_mff_data.read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr)     

len(hdr['orig']['epochLabels'])

