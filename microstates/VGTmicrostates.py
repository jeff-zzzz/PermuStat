# -*- coding: utf-8 -*-
"""
Created on Mon Apr 27 11:18:27 2015 @author: jesong1126
"""
 
import numpy as np   
import matplotlib.pyplot as plt
from matplotlib import gridspec
from mff import read_mff_header, read_mff_data  #, getEpochInfos, mff_getSummaryInfo
from microstates import Microstates #,Kmeans  
from pandas import DataFrame 
import pandas as pd 

gsn257 = pd.read_table('/Users/jesong1126/Python27/GeoPy/nscolor_hgsn/GSN257ToPy2.dat')  
frame = DataFrame(gsn257 , columns=['ChLabel', 'X3', 'Y3', 'Z3','X2','Y2'])
x2 = frame.values[:,4] 
y2 = frame.values[:,5] 

##-----------------------------------------------------------------------------
## VGT_105_20140307_095017_A800_bcr_ref_blc.mff 
## VGT_121_A800_bcr_ref_blc.mff 
## VGT_125_20141111_022603_A800_bcr_ref_blc.mff 
## VGT_130_20150126_030959_fil_segA_blc_bcr.mff 

filePath ='/Users/jesong1126/Python27/data_GeoPy/VGT/VGT_105_20140307_095017_A800_bcr_ref_blc.mff'
filePath ='/Users/jesong1126/Python27/data_GeoPy/VGT/VGT_121_A800_bcr_ref_blc.mff'
filePath ='/Users/jesong1126/Python27/data_GeoPy/VGT/VGT_125_20141111_022603_A800_bcr_ref_blc.mff'
filePath ='/Users/jesong1126/Python27/data_GeoPy/VGT/VGT_130_20150126_030959_fil_segA_blc_bcr.mff'

hdr = read_mff_header.read_mff_header(filePath)
nC = hdr['nChans']
nSamples = hdr['nSamples']
nSamplesPre = hdr['nSamplesPre']
nTrials = hdr['nTrials']
srate = hdr['Fs']
summaryInfo = hdr['orig'] 
trialsName = summaryInfo['epochLabels']   
categoryName = list(set(trialsName))
nCategory = len(categoryName)

data = read_mff_data.read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr)    
baseline = (nSamplesPre * 1000 / srate)  
msSamples = np.arange(0, nSamples,1) * 1000/srate  - baseline 
xlimMin = msSamples[0]
xlimMax = msSamples[-1]

len(hdr['orig']['epochLabels'])
SegStatus = hdr['orig']['epochSegStatus']  
GoodSeg = []; BadSeg = []; epochLabels = []; 
for i in range(len(SegStatus)):
    if SegStatus[i] == 'bad' :  
        BadSeg.append(i)
    else :  
        GoodSeg.append(i)
nTrials = len(GoodSeg)

## average reference 
H = np.identity(nC) - np.ones((nC, nC))/nC  
s = np.zeros((data.shape[0], data.shape[1], nTrials))  
sGFP = np.zeros((nSamples, nTrials))
if nTrials > 1:
    for i in range(nTrials):
        s[:,:,i] = np.dot(H, data[:,:,GoodSeg[i]])
        sGFP[:,i] = np.std(s[:,:,i], axis=0)    
else :
    s = np.dot(H, data)  
    sGFP = np.std(s, axis=0)  

X = s ; nX = nTrials    
Xmean = np.mean(X, axis=2)
Xstd = np.std(X, axis=2)
Xvar = np.var(X, axis=2)
XT = np.sqrt(nX) * Xmean / Xstd  

##-----------------------------------------------------------------------------
k=8
nSim = 1000
Data = XT[:, nSamplesPre:]
ChangeOut = Microstates.Change(Data, k, nSim) 
Tmaps = ChangeOut['Tmaps'] 
KmeanId = np.array(ChangeOut['KmeanId'], 'i4')
GEVs= ChangeOut['GEVs'] 
CutId =ChangeOut['CutId']

GFP = np.std(XT, axis=0) 
nT = Data.shape[1] 

Mycolorkeys = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','purple','firebrick','darkgoldenrod','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00','b', 'g', 'r', 'c', 'm', 'y', 'k') 
KmeanColorId = [] # [None] * nT 
for i in range(nT):
    KmeanColorId.append(Mycolorkeys[KmeanId[i]])  

##-----------------------------------------------------------------------------
fig = plt.figure(figsize=(16, 8)) 
gs = gridspec.GridSpec(2, 1, height_ratios=[3,1]) 
ax0 = plt.subplot(gs[0])
ax0.plot(msSamples, XT.T)
ax0.set_title('XT')
ax0.set_xlim([msSamples[0], msSamples[-1]])
ax1 = plt.subplot(gs[1])
ax1.plot(msSamples, GFP)
ax1.set_xlim([msSamples[0], msSamples[-1]])
ax1.set_title('gfp ')
for i in range(k):
    ax1.axvline(msSamples[nSamplesPre+CutId[i]], color=Mycolorkeys[i] ,  linewidth=2)
    ax1.text(msSamples[nSamplesPre+(CutId[i]+CutId[i+1])/2], 1, ('%d' % i), fontsize=15)    
#    ax1.text(msSamples[nSamplesPre+(CutId[i]+CutId[i+1])/2],  GFP[nSamplesPre+(CutId[i]+CutId[i+1])/2]-1, ('%d' % i), fontsize=15)    
ax1.vlines(msSamples[nSamplesPre:], [0], GFP[nSamplesPre:], color=KmeanColorId, alpha= 0.2)    
plt.show()
plt.savefig('VGT130MicroStateGFP.png')


fig,axes = plt.subplots(1, k, figsize=(3*k, 3)) #, sharex=True)
for ii in range(k):
    axes[ii].scatter(x2, y2, c=Tmaps[:,ii], s=30, cmap=plt.get_cmap('seismic'), alpha= .5)
    axes[ii].set_alpha(0.75)
    axes[ii].set_title(('MicroState %d' % ii))
    axes[ii].set_xticks([])
    axes[ii].set_yticks([])
plt.savefig('VGT130MicroStateTMaps.png')




##-----------------------------------------------------------------------------
from braink import read_lfm
from swcm import Imatrix      
nE = nC -1  
K  = read_lfm.lfm('/Users/jesong1126/Python27/data_GeoPy/108_HM/Leadfield.lfm', nE)                                  
nV = K.shape[1]

alpha = 0.001 
Imat = Imatrix.MN(alpha, K)

sdenTmaps = np.zeros((nV, Tmaps.shape[1])) 
for i in range( Tmaps.shape[1]):             
    sdenTmaps[:,i] = np.dot(Imat, Tmaps[:,i]) 

import os
os.system("open /Applications/EAV/EAV.app")

from RabbitMQ import Connection 

c = Connection.Connection('localhost')
c.connect()

i=0
c.sendEEGData(Tmaps[:, i], "Tmap", "")
c.sendOrientedData(abs(sdenTmaps[:, i]), "sdenTmap", "")

c.disconnect()       

##-----------------------------------------------------------------------------

#    plt.savefig('VGT%sMicroStateTMaps.png' %subjectName[sIdx])
#    plt.close() 

#ax2 = plt.subplot(gs[2])
#ax2.scatter(msSamples[nSamplesPre:], KmeanId) #ax2.set_title('mean gdiss(lap=%d)'%lap)
#ax2.set_xlim([msSamples[0], msSamples[-1]])
#ax3 = plt.subplot(gs[3])
#ax3.plot(GEVs) #, color='k') 
 #Data = XT  

#nSim = 1000
#CV = np.zeros((20,2))
#for k in range(2, 20):
#    ChangeOut = Microstates.Change(Data, k, nSim)
#    GEVs= ChangeOut['GEVs'] 
#    CV[k, 0]= k 
#    CV[k, 1]= GEVs[-1]
#
#CVslope = (np.append(np.diff(CV[:,1]),0 ) + np.append(0, np.diff(CV[:,1])))/2
#
#fig, ax= plt.subplots(2,1)
#ax[0].scatter(CV[2:,0], CV[2:,1])
#ax[0].set_xlabel('# of cluster')
#ax[0].set_ylabel('GEV')
#ax[1].plot(CV[2:,0], CVslope[2:])
#ax[1].set_xlim([0, 20])
#ax[1].set_xlabel('# of cluster')
#ax[1].set_ylabel('sclope')
#plt.show()
#
