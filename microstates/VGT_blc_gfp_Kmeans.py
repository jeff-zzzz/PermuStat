# -*- coding: utf-8 -*-
#import config 

import numpy as np   
import matplotlib.pyplot as plt

from mff import read_mff_header, read_mff_data  #, getEpochInfos, mff_getSummaryInfo
#filePath ='/Volumes/Jasmine 1/Work/Data/VGT/VGT_8subj_bcr_blc_ave.mff' #this has too much bad channels 
#filePath ='/Users/jesong1126/Python27/GeoPy/data/VGT/VGT_105_20140307_095017_A800_bcr_ref_blc.mff'

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
epochName = summaryInfo['epochFilenames']
subjectName = []  
for i in range(nTrials):
    subjectName.append(epochName[i][:7])

data = read_mff_data.read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr)    
baseline = (nSamplesPre * 1000 / srate)  
msSamples = np.arange(0, nSamples,1) * 1000/srate  - baseline 
xlimMin = msSamples[0]
xlimMax = msSamples[-1]

## average reference 
H = np.identity(nC) - np.ones((nC, nC))/nC 
s = np.zeros(data.shape)  
sGFP = np.zeros((nSamples, nTrials))
if nTrials > 1:
    for i in range(nTrials):
        s[:,:,i] = np.dot(H, data[:,:,i])
        sGFP[:,i] = np.std(s[:,:,i], axis=0)    
else :
    s = np.dot(H, data)  
    sGFP = np.std(s, axis=0)  
    
ax = plt.subplots(1,2)    
ax[0].plot(msSamples,sGFP)
ax[0].xlim([msSamples[0], msSamples[-1]])
ax[0].legend(subjectName, loc='upper left')
#plt.savefig('VGT8subjGFP.png')
plt.show()
      
## baseline correction 
s_blc = np.zeros(data.shape)    
s_blc_GFP = np.zeros(sGFP.shape)
if nTrials > 1:
    for i in range(nTrials):
        Data = s[:,:,i] 
        BaselineData = Data[:, 0:nSamplesPre]
        BaselineData0 = np.mean(BaselineData, axis=1)
        s0 = Data - BaselineData0.reshape(BaselineData0.shape[0],1) 
        s_blc[:,:,i] = s0
        s_blc_GFP[:,i] = np.std(s0, axis=0) 
else :
    s_blc = np.dot(H, s)  
    s_blc_GFP = np.std(s_blc, axis=0)  
 
fig = plt.figure() 
plt.plot(msSamples,s_blc_GFP)
plt.xlim([msSamples[0], msSamples[-1]])
plt.legend(subjectName, loc='upper left')
#plt.savefig('VGT8subjGFP.png')
plt.show()

sIdx = 0
#Data = s_blc[:,:,sIdx] 
#s_blc_GDISS = np.zeros((Data.shape[1],Data.shape[1])) 
#for j in range(Data.shape[1]-1):
#    for i in range(j+1, Data.shape[1]):
#        u = Data[:,j] 
#        v = Data[:,i]
#        u_modulus = np.std(u)  
#        v_modulus = np.std(v)  
#        cos_angle = np.dot(u,v) / (u_modulus * v_modulus) /nC #-> cosine of the angle
#        s_blc_GDISS[i,j] = (1 - cos_angle)
##        s_blc_GDISS[i,j] = cos_angle


#s_blc_GDISS = np.zeros((Data.shape[1] ,nTrials)) 
#for i in range(nTrials):
#    Data = s_blc[:,:,i] 
#    for j in range(Data.shape[1]-1):
#        u = Data[:,j] 
#        v = Data[:,j+1]
#        u_modulus = np.std(u)  
#        v_modulus = np.std(v)  
#        cos_angle = np.dot(u,v) / (u_modulus * v_modulus) /nC #-> cosine of the angle
#        s_blc_GDISS[j,i] = (1 - cos_angle)
#
#fig, ax = plt.subplots(2,1)    
#ax[0].plot(s_blc_GFP[:1000,:])            
#ax[1].plot(s_blc_GDISS[:1000,:])            
#plt.show() 


Data = s_blc[:,:,sIdx] 
s_blc_GDISS = np.zeros((Data.shape[1],Data.shape[1])) 
for j in range(Data.shape[1]-1):
    for i in range(j+1, Data.shape[1]):
        u = Data[:,j] 
        v = Data[:,i]
        u_modulus = np.std(u)  
        v_modulus = np.std(v)  
        cos_angle = np.dot(u,v) / (u_modulus * v_modulus) /nC #-> cosine of the angle
        s_blc_GDISS[i,j] = (1 - cos_angle)

Similarity = s_blc_GDISS + s_blc_GDISS.T
plt.plot(np.mean(Similarity, axis=0))


fig, ax = plt.subplots(8,1)
for i in range(nTrials):
    sss = s_blc[:,:,i].T
    ax[i].plot(sss)

plt.show()


from matplotlib import gridspec
fig = plt.figure(figsize=(8, 12)) 
gs = gridspec.GridSpec(3, 1, height_ratios=[1,1,4]) 
ax0 = plt.subplot(gs[0])
ax0.plot(s_blc_GFP[:,sIdx])
ax0.set_title('gfp')
ax1 = plt.subplot(gs[1])
ax1.plot(np.mean(Similarity, axis=0))
ax1.set_title('ms vs mean(gdiss)')
ax2 = plt.subplot(gs[2])
#im2 = ax2.imshow(Similarity, cmap='jet')# , aspect=0)
im2 = ax2.imshow(Similarity, cmap='hot')
#plt.colorbar(im1)
plt.tight_layout()
ax2.set_title('dissimilarity')
plt.show()



Tmaps = Data[:, [300,500]]
k = 2
nT = Data.shape[1] 
    # Correlation between Tmaps and Samples 
Cuv = np.zeros((k, nT))  
for j in range(nT) :
    for i in range(k) :
        u = Tmaps[:,i] 
        v = Data[:,j]
        u_modulus = np.sqrt((u*u).sum()) 
        v_modulus = np.sqrt((v*v).sum())
        cos_angle= np.dot(u,v) / u_modulus / v_modulus  
        Cuv[i, j] = cos_angle 
            
plt.figure()
plt.plot(Cuv.T)            
plt.show() 

KmeanId = np.zeros((nT,1)) 
for j in range(nT):
    KmeanId[j] = np.where((Cuv[:,j])== max(Cuv[:,j]))[0]
KmeanId = np.array(KmeanId, dtype='i4')


#im4 = ax4.imshow(tot2, norm=LogNorm(vmin=0.001, vmax=1), aspect='auto')
#divider4 = make_axes_locatable(ax4)
#cax4 = divider4.append_axes("right", size="20%", pad=0.05)
#cbar4 = plt.colorbar(im4, cax=cax4

#    u = Data[:,200] 
#    v = Data[:,700]
#    u_modulus = np.std(u)  
#    v_modulus = np.std(v)  
#    cos_angle = np.dot(u,v) / (u_modulus * v_modulus) /nC #-> cosine of the angle
#    cos_angle

                             
fig, ax = plt.subplots(3,1) 
ax[0].plot(msSamples,sGFP)
ax[0].set_xlim([msSamples[0], msSamples[-1]])
ax[0].set_title('gfp') 
# ax[0].legend(('105','106','107', '110' ,'111','112','113','117'), loc='upper left')
ax[1].plot(msSamples, s_blc_GFP)
ax[1].set_xlim([msSamples[0], msSamples[-1]])
ax[1].set_title('baseline correced gfp')
# ax[1].legend(('105','106','107', '110' ,'111','112','113','117'), loc='upper left')
#plt.savefig('VGT8subjbcr_GFP.png')
ax[2].plot(msSamples, 100*s_blc_GDISS)
ax[2].set_xlim([msSamples[0], msSamples[-1]])
ax[2].set_title('global dissimilarity')
plt.show() 


from pandas import DataFrame 
import pandas as pd 

gsn257 = pd.read_table('/Users/jesong1126/Python3/gs3/nscolor_hgsn/GSN257ToPy.sfp') 
frame = DataFrame(gsn257, columns=['ChLabel', 'X3', 'Y3', 'Z3','X2','Y2'])
x2 = frame.values[:,4] 
y2 = frame.values[:,5] 

colorkeys = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','firebrick','darkgoldenrod','purple','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00')
Kmeans(Data, k, nStable) 

from microstates import Kmeans  
#SubjId =('105','106','107', '110' ,'111','112','113','117') 
k0 = 5 
nSim = 1000 
StableLegnth = 50 # ms 

for sIdx in range(nTrials): 
    Data = s_blc[:,:,sIdx] 
    blc_gfp = s_blc_GFP[:,sIdx] 
    blc_diss = s_blc_GDISS[:,sIdx] 

    plt.plot(msSamples,Data.T) 
    plt.xlim([msSamples[0], msSamples[-1]])
    plt.savefig('VGT%sButterfly.png' %subjectName[sIdx])
    plt.close() 

    KmeansOut = Kmeans.Kmeans(Data, k0, nSim)  
    Tmaps = KmeansOut['Tmaps'] 
    KmeanId = np.array(KmeansOut['KmeanId'], 'i4')
    GEV = KmeansOut['GEV']
    GEVs= KmeansOut['GEVs'] 
    C_UTmaps = KmeansOut['C_UTmaps'] 
    KmeanColorId = KmeansOut['KmeanColorId']

    ChangePt = np.concatenate([[0], 1+np.where(np.diff(KmeanId)!=0)[0], [nSamples-1]])
    StateDuration = np.diff(ChangePt) 

    plt.plot(msSamples,blc_gfp, color='k')
    plt.vlines(msSamples, [0], blc_gfp, color=KmeanColorId, alpha= 0.2)    
    plt.xlim([msSamples[0], msSamples[-1]])    
    plt.title(('GFP and %d groups (GEV=%.2f)' % (k0, GEV)))
    plt.savefig('VGT%sKmeanInitial.png' %subjectName[sIdx])
    plt.close() 
    
    aaa = np.where(StateDuration > StableLegnth)[0]
    nMS = len(aaa) 
    BeginMS = ChangePt[aaa]
    EndMS = ChangePt[aaa+1]
    Microstate = np.zeros((nC, nMS))
    for iii in range(nMS): 
        DataI = Data[:,BeginMS[iii]:EndMS[iii]] 
        DataI = np.reshape(DataI, (nC,-1)) 
        Microstate[:,iii] = np.mean(DataI, axis=1) 
    

    plt.plot(msSamples,blc_gfp, color='k')
    plt.xlim([msSamples[0], msSamples[-1]])    
    for ii in range(nMS):
        tempii = range(BeginMS[ii], EndMS[ii])     
        plt.vlines(msSamples[tempii], [0], blc_gfp[tempii], color=colorkeys[ii], alpha= 0.2)    
        plt.text(msSamples[tempii[round(len(tempii)/2)]], 1, ('%d' % ii), fontsize=15)    
    
    plt.savefig('VGT%sGfpMicroState.png' %subjectName[sIdx])
    plt.close() 

             
    fig,axes = plt.subplots(1, nMS, figsize=(3*nMS,3)) #, sharex=True)
    for ii in range(nMS):
        axes[ii].scatter(x2, y2, c=Microstate[:,ii], s=30, cmap=plt.get_cmap('seismic'), alpha= .5)
        axes[ii].set_alpha(0.75)
        axes[ii].set_title(('MicroState %d' % ii))
        axes[ii].set_xticks([])
        axes[ii].set_yticks([])

    plt.savefig('VGT%sMicroStateTMaps.png' %subjectName[sIdx])
    plt.close() 





  
    
#    ##
#from mff import KmeansRJ 
#
#Data = s_blc[:,:,0] 
#blc_gfp = s_blc_GFP[:,0] 
#
#plt.plot(msSamples,Data.T) 
#plt.xlim([msSamples[0], msSamples[-1]])
#plt.savefig('VGT105Butterfly.png')
#plt.show()
#
#k0 = 5 
#nSim = 10000
#RJout = KmeansRJ.KmeansRJ(Data, k0, nSim) 
##RJout.keys()
#KmeanId = np.array(RJout['KmeanId'], 'i4')
#GEV = RJout['GEV']
#GEVs = RJout['GEVs'] 
#PartitionId = RJout['PartitionId'] 
#Cuv = RJout['Cuv'] 
#Tmaps = RJout['Tmaps'] 
#KmeanColorId= RJout['KmeanColorId']
#
#plt.plot(msSamples,blc_gfp, color= 'k')
#plt.vlines(msSamples, [0], blc_gfp, color=KmeanColorId)    
#plt.title(('GFP and %d groups (GEV=%.2f)' % (k0, GEV)))
#plt.show()
#
#
#from mff import KmeansMixed 
#
#k0 = 5 
#nSim = 1000
#Mixedout = KmeansMixed.KmeansMixed(Data, k0, nSim) 
#KmeanIdMixed = np.array(Mixedout['KmeanId'], 'i4')
#GEVMixed = Mixedout['GEV']
#GEVsMixed = Mixedout['GEVs'] 
#PartitionIdMixed = Mixedout['PartitionId'] 
#CuvMixed = Mixedout['Cuv'] 
#TmapsMixed = Mixedout['Tmaps'] 
#KmeanColorIdMixed = Mixedout['KmeanColorId']
#
#plt.plot(msSamples,blc_gfp, color= 'k')
#plt.vlines(msSamples, [0], blc_gfp, color=KmeanColorIdMixed)    
#plt.title(('GFP and %d groups (GEV=%.2f)' % (k0, GEVMixed)))
#plt.show()
#

    
    
#plt.scatter(msSamples, KmeanId)
#plt.show()
#nT =1600 
#KmeanId[:200]

#
#
##sss = pd.read_table('/Users/jesong1126/Python3/gs3/nscolor_hgsn/Seismic.clr', header=None) 
##sssframe = DataFrame(sss)
##sssframe
##aaa = list(open('/Users/jesong1126/Python3/gs3/nscolor_hgsn/Seismic.clr'))
##sesmic = sssframe.values
# 
#Tmaps0 = Tmaps 
#Tmaps = Tmaps0 - np.mean(Tmaps0)
# 
#fig,axes = plt.subplots(1, kOpt, figsize=(15,3))#, sharex=True)
#axes[0].scatter(x2, y2, c=Tmaps[:,0], s=20, cmap=plt.get_cmap('seismic'))
#axes[0].set_alpha(0.75)
#axes[0].set_title('Cluster 1')
#axes[0].set_xticks([])
#axes[0].set_yticks([])
#
#axes[1].scatter(x2, y2, c=Tmaps[:,1],s=20, cmap=plt.get_cmap('seismic'))
#axes[1].set_title('Cluster 2')
#axes[1].set_xticks([])
#axes[1].set_yticks([])
#
#axes[2].scatter(x2, y2, c=Tmaps[:,2],s=20, cmap=plt.get_cmap('seismic'))
#axes[2].set_title('Cluster 3')
#axes[2].set_xticks([])
#axes[2].set_yticks([])
#
#axes[3].scatter(x2, y2, c=Tmaps[:,3],s=20, cmap=plt.get_cmap('seismic'))
#axes[3].set_title('Cluster 4')
#axes[3].set_xticks([])
#axes[3].set_yticks([])
#
#axes[4].scatter(x2, y2, c=Tmaps[:,4],s=20, cmap=plt.get_cmap('seismic'))
#axes[4].set_title('Cluster 5') 
#axes[4].set_xticks([])
#axes[4].set_yticks([])
##plt.savefig('VGT106kTMaps.png')
#plt.show()
#
#
#fig = plt.figure()
#ax1 = plt.subplot2grid((3,5), (0,0), colspan=5)
#ax1.plot(msSamples, Data.T)  
#ax1.set_ylabel(' ')
#ax1.set_xlim([msSamples[0], msSamples[-1]])
#ax2 = plt.subplot2grid((3,5), (1,0), colspan=5)
#ax2.vlines(msSamples, [0], GFP, color=KmeanIDcolor)  
#ax2.set_xlim([msSamples[0], msSamples[-1]])
#ax2.set_xlabel('ms')
#ax2.set_ylabel('GFP ')
#ax2.text(0, 4, '1', color='blue',fontsize=15)
#ax2.text(500, 3.5, '2', color='green',fontsize=15)
#ax2.text(750, 3.5, '3', color='red',fontsize=15)
#ax2.text(1000, 3, '4', color='c',fontsize=15)
#ax2.text(1300, 6, '5', color='m',fontsize=15)
#ax3 = plt.subplot2grid((3,5), (2,0)) #, rowspan=2)
#ax3.scatter(x2, y2, c=Tmaps[:,0], s=20, cmap=plt.get_cmap('seismic'))
#ax3.set_alpha(0.75)
#ax3.set_xlabel('Cluster 1')
#ax3.spines['top'].set_color('blue')
#ax3.spines['bottom'].set_color('blue')
#ax3.spines['left'].set_color('blue')
#ax3.spines['right'].set_color('blue')
#ax3.set_xticks([])
#ax3.set_yticks([])
#ax4 = plt.subplot2grid((3,5), (2,1))
#ax4.scatter(x2, y2, c=Tmaps[:,1], s=20, cmap=plt.get_cmap('seismic'))
#ax4.set_alpha(0.75)
#ax4.set_xlabel('Cluster 2')
#ax4.spines['top'].set_color('green')
#ax4.spines['bottom'].set_color('green')
#ax4.spines['left'].set_color('green')
#ax4.spines['right'].set_color('green')
#ax4.set_xticks([])
#ax4.set_yticks([])
#ax5 = plt.subplot2grid((3,5), (2,2))
#ax5.scatter(x2, y2, c=Tmaps[:,2], s=20, cmap=plt.get_cmap('seismic'))
#ax5.set_alpha(0.75)
#ax5.set_xlabel('Cluster 3')
#ax5.spines['top'].set_color('red')
#ax5.spines['bottom'].set_color('red')
#ax5.spines['left'].set_color('red')
#ax5.spines['right'].set_color('red')
#ax5.set_xticks([])
#ax5.set_yticks([])
#ax6 = plt.subplot2grid((3,5), (2,3))
#ax6.scatter(x2, y2, c=Tmaps[:,3], s=20, cmap=plt.get_cmap('seismic'))
#ax6.set_alpha(0.75)
#ax6.set_xlabel('Cluster 4')
#ax6.spines['top'].set_color('c')
#ax6.spines['bottom'].set_color('c')
#ax6.spines['left'].set_color('c')
#ax6.spines['right'].set_color('c')
#ax6.set_xticks([])
#ax6.set_yticks([])
#ax7 = plt.subplot2grid((3,5), (2,4))
#ax7.scatter(x2, y2, c=Tmaps[:,4], s=20, cmap=plt.get_cmap('seismic'))
#ax7.set_alpha(0.75)
#ax7.set_xlabel('Cluster 5')
#ax7.spines['top'].set_color('m')
#ax7.spines['bottom'].set_color('m')
#ax7.spines['left'].set_color('m')
#ax7.spines['right'].set_color('m')
#ax7.set_xticks([])
#ax7.set_yticks([])
##plt.savefig('VGT106kButterGfpTMaps.png')
#plt.show()
# 
# 
 
 
###------------------------------------------------------------------------------
#colorkeys = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','firebrick',
#'darkgoldenrod','purple','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00') 
#  
#
#import random 
#PartitionId = np.array([   0,  171,  180,  329,  448, 1600])
#
#
##MaxGEVAt = np.where(GEVs == max(GEVs))[0] 
##plt.plot(range(nSim), GEVs, color = 'k')
##plt.vlines(MaxGEVAt[0], [0], GEVs[MaxGEVAt[0]], color='r')    
##plt.title(('GEVs converge at %d with %d groups (GEV=%.2f)' % (MaxGEVAt[0], k0, GEV)))
##plt.show()
##
#plt.plot(KmeanId) 
#plt.ylim([-1, k0])
#plt.show()
#
#PartitionId    
#
#BaseMean = np.mean(blc_gfp[0:nSamplesPre]) 
#BaseStd = np.std(blc_gfp[0:nSamplesPre])  
#BaseCI99 =(BaseMean+2.575*BaseStd) 
#
#plt.plot(msSamples,blc_gfp)
#plt.xlim([msSamples[0], msSamples[-1]])
#plt.hlines(BaseMean, msSamples[0], msSamples[-1], color='r')
#plt.hlines((BaseMean+2.575*BaseStd), msSamples[0], msSamples[-1], color='k')
#plt.hlines((BaseMean-2.575*BaseStd), msSamples[0], msSamples[-1], color='k')
#plt.title('VGT 105')
##plt.legend(('gfp','rmse'), loc='best')
##plt.savefig('VGT_blc_gfp_105.png')
#plt.show()

## Kmeans clustring 
             
## 
#for nnn in range(100): 
#    
#SelId = 2 + random.sample(range(k-1), 1)[0]   #2~5(k-1)                
#PartitionIdProp = np.delete(PartitionId, SelId)
#PartitionIdProp = np.sort(np.concatenate((PartitionIdProp, random.sample(range(FirstNonBaseline, nSamples), 1)), axis=0))   
#             
#IdProp = np.ones((PartitionIdProp.shape[0]-1, nSamples))
#KmeanIdProp = np.array([]) 
#for i in range(PartitionIdProp.shape[0]-1):
#    KmeanIdProp= np.concatenate((KmeanIdProp, (i) * IdProp[i, (PartitionIdProp[i]):(PartitionIdProp[i+1])]))
#    
#KmeanIdProp = np.array(KmeanIdProp, 'i4')
#
#TmapsProp = np.zeros((nC, PartitionIdProp.shape[0]-1))
#for i in range(PartitionIdProp.shape[0]-1):
#    TmapsProp[:,i] = np.mean(Data[:, (PartitionIdProp[i]):(PartitionIdProp[i+1])], axis=1)
#        
## Correlation between Tmaps and Samples 
#CuvProp = np.array([]) #np.zeros((1, nSamples))  
#for i in range(k+1) :
#    DataI = Data[:, range(PartitionIdProp[i], PartitionIdProp[i+1])]
#    CuvI = np.zeros((DataI.shape[1]))
#    u = Tmaps[:,i] 
#    u_modulus = np.sqrt((u*u).sum()) 
#    for j in range(DataI.shape[1]): 
#        v = DataI[:,j]
#        v_modulus = np.sqrt((v*v).sum())
#        cos_angle= np.dot(u,v) / u_modulus / v_modulus  
#        CuvI[j] = cos_angle
#    CuvProp = np.concatenate((CuvProp, CuvI))   
#       
#GEVProp = sum(blc_gfp * CuvProp * blc_gfp * CuvProp)/sum(blc_gfp*blc_gfp )  
#GEVProp
#
#if GEVProp > GEV: 
#    GEV = GEVProp
#    PartitionId = PartitionIdProp
#    KmeanId = KmeanIdProp
#    Tmaps = TmapsProp
#

#plt.plot(msSamples, abs(Cuv)); #plt.ylim([-.2, 1.2])    
#plt.vlines(msSamples, [0], abs(Cuv), color=KmeanIdColor)    
#plt.show()  

##----------------------------------------------------------------------------

#Data = DataAll[:, FirstNonBaseline:]
#VGT_Kmeans = Kmeans.Kmeans(Data, k, nStable) 
#GEV = VGT_Kmeans['GEV'] 
#GEVs = VGT_Kmeans['GEVs']
#C_UTmaps = VGT_Kmeans['C_UTmaps']
#KmeanIDPost = VGT_Kmeans['KmeanID']
#Tmaps = VGT_Kmeans['Tmaps']
#
#GEV 
#
#colorkeys = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','firebrick',
#'darkgoldenrod','purple','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00')
#
#KmeanID0 = k * np.ones((FirstNonBaseline,), 'i4') 
#KmeanID = np.concatenate((KmeanID0,  KmeanIDPost), axis=0) 
#
#KmeanIDcolor = [None] * KmeanID.shape[0]
#for i in range(KmeanID.shape[0]):
#    KmeanIDcolor[i] = colorkeys[KmeanID[i]]
#
#fig,axes = plt.subplots(2,1)#, sharex=True)
#axes[0].plot(msSamples, DataAll.T)  
#axes[0].set_ylabel(' ')
#axes[0].set_xlim([msSamples[0], msSamples[-1]])
#axes[1].vlines(msSamples, [0], blc_gfp, color=KmeanIDcolor)  
#axes[1].set_xlim([msSamples[0], msSamples[-1]])
#axes[1].set_xlabel('ms')
#axes[1].set_ylabel('GFP and 8 clusters')
##plt.savefig('VGT_blc_Kmeans105.png')
#plt.show()

#SlopeRmse= np.diff(VGTrmse) 
#Thres = np.where(VGTrmse > BaseCI99)[0] 
#PosSlope = np.where(SlopeRmse>0)[0] 
#Thres = (VGTrmse > BaseCI99)
#PosSlope = (SlopeRmse>0)
#Thres2= Thres[1:] * PosSlope
#
#PosSlopeInd = np.where(SlopeRmse>0)[0] 
#SlopeRmse[PosSlopeInd[:(-1)]]
#PosNegSlope = (SlopeRmse[PosSlopeInd[:(-1)]+1] < 0) 
#PosNegSlopeInd = np.where(SlopeRmse[PosSlopeInd[:(-1)]+1] < 0)[0] 
#
#plt.scatter(msSamples[range(1,nSamples)], Thres2,  s=2)
#plt.show()

#fig, ax= plt.subplots(2,1)
#ax[0].plot(msSamples,VGTrmse)
#ax[0].xlim([msSamples[0], msSamples[-1]])
#ax[0].hlines(BaseMean, msSamples[0], msSamples[-1], color='r')
#ax[0].hlines((BaseMean+2.575*BaseStd), msSamples[0], msSamples[-1], color='g')
#ax[0].hlines((BaseMean-2.575*BaseStd), msSamples[0], msSamples[-1], color='g')
#
#ax[1].scatter(msSamples[range(1,nSamples)], SlopeRmse,  s=2)
#ax[1].hlines(0, msSamples[1], msSamples[-1], color='r')
#ax[1].xlim([msSamples[0], msSamples[-1]])
##plt.savefig('VGT106RmseBaseline.png')
#plt.show()





#
#from pandas import DataFrame 
#import pandas as pd 
#
#gsn257 = pd.read_table('/Users/jesong1126/Python3/gs3/nscolor_hgsn/GSN257ToPy.sfp') 
#frame = DataFrame(gsn257, columns=['ChLabel', 'X3', 'Y3', 'Z3','X2','Y2'])
#x2 = frame.values[:,4] 
#y2 = frame.values[:,5] 
#
##sss = pd.read_table('/Users/jesong1126/Python3/gs3/nscolor_hgsn/Seismic.clr', header=None) 
##sssframe = DataFrame(sss)
##sssframe
##aaa = list(open('/Users/jesong1126/Python3/gs3/nscolor_hgsn/Seismic.clr'))
##sesmic = sssframe.values
# 
#Tmaps0 = Tmaps 
#Tmaps = Tmaps0 - np.mean(Tmaps0)
# 
#fig,axes = plt.subplots(1, kOpt, figsize=(15,3))#, sharex=True)
#axes[0].scatter(x2, y2, c=Tmaps[:,0], s=20, cmap=plt.get_cmap('seismic'))
#axes[0].set_alpha(0.75)
#axes[0].set_title('Cluster 1')
#axes[0].set_xticks([])
#axes[0].set_yticks([])
#
#axes[1].scatter(x2, y2, c=Tmaps[:,1],s=20, cmap=plt.get_cmap('seismic'))
#axes[1].set_title('Cluster 2')
#axes[1].set_xticks([])
#axes[1].set_yticks([])
#
#axes[2].scatter(x2, y2, c=Tmaps[:,2],s=20, cmap=plt.get_cmap('seismic'))
#axes[2].set_title('Cluster 3')
#axes[2].set_xticks([])
#axes[2].set_yticks([])
#
#axes[3].scatter(x2, y2, c=Tmaps[:,3],s=20, cmap=plt.get_cmap('seismic'))
#axes[3].set_title('Cluster 4')
#axes[3].set_xticks([])
#axes[3].set_yticks([])
#
#axes[4].scatter(x2, y2, c=Tmaps[:,4],s=20, cmap=plt.get_cmap('seismic'))
#axes[4].set_title('Cluster 5') 
#axes[4].set_xticks([])
#axes[4].set_yticks([])
##plt.savefig('VGT106kTMaps.png')
#plt.show()
#
#
#fig = plt.figure()
#ax1 = plt.subplot2grid((3,5), (0,0), colspan=5)
#ax1.plot(msSamples, Data.T)  
#ax1.set_ylabel(' ')
#ax1.set_xlim([msSamples[0], msSamples[-1]])
#ax2 = plt.subplot2grid((3,5), (1,0), colspan=5)
#ax2.vlines(msSamples, [0], GFP, color=KmeanIDcolor)  
#ax2.set_xlim([msSamples[0], msSamples[-1]])
#ax2.set_xlabel('ms')
#ax2.set_ylabel('GFP ')
#ax2.text(0, 4, '1', color='blue',fontsize=15)
#ax2.text(500, 3.5, '2', color='green',fontsize=15)
#ax2.text(750, 3.5, '3', color='red',fontsize=15)
#ax2.text(1000, 3, '4', color='c',fontsize=15)
#ax2.text(1300, 6, '5', color='m',fontsize=15)
#ax3 = plt.subplot2grid((3,5), (2,0)) #, rowspan=2)
#ax3.scatter(x2, y2, c=Tmaps[:,0], s=20, cmap=plt.get_cmap('seismic'))
#ax3.set_alpha(0.75)
#ax3.set_xlabel('Cluster 1')
#ax3.spines['top'].set_color('blue')
#ax3.spines['bottom'].set_color('blue')
#ax3.spines['left'].set_color('blue')
#ax3.spines['right'].set_color('blue')
#ax3.set_xticks([])
#ax3.set_yticks([])
#ax4 = plt.subplot2grid((3,5), (2,1))
#ax4.scatter(x2, y2, c=Tmaps[:,1], s=20, cmap=plt.get_cmap('seismic'))
#ax4.set_alpha(0.75)
#ax4.set_xlabel('Cluster 2')
#ax4.spines['top'].set_color('green')
#ax4.spines['bottom'].set_color('green')
#ax4.spines['left'].set_color('green')
#ax4.spines['right'].set_color('green')
#ax4.set_xticks([])
#ax4.set_yticks([])
#ax5 = plt.subplot2grid((3,5), (2,2))
#ax5.scatter(x2, y2, c=Tmaps[:,2], s=20, cmap=plt.get_cmap('seismic'))
#ax5.set_alpha(0.75)
#ax5.set_xlabel('Cluster 3')
#ax5.spines['top'].set_color('red')
#ax5.spines['bottom'].set_color('red')
#ax5.spines['left'].set_color('red')
#ax5.spines['right'].set_color('red')
#ax5.set_xticks([])
#ax5.set_yticks([])
#ax6 = plt.subplot2grid((3,5), (2,3))
#ax6.scatter(x2, y2, c=Tmaps[:,3], s=20, cmap=plt.get_cmap('seismic'))
#ax6.set_alpha(0.75)
#ax6.set_xlabel('Cluster 4')
#ax6.spines['top'].set_color('c')
#ax6.spines['bottom'].set_color('c')
#ax6.spines['left'].set_color('c')
#ax6.spines['right'].set_color('c')
#ax6.set_xticks([])
#ax6.set_yticks([])
#ax7 = plt.subplot2grid((3,5), (2,4))
#ax7.scatter(x2, y2, c=Tmaps[:,4], s=20, cmap=plt.get_cmap('seismic'))
#ax7.set_alpha(0.75)
#ax7.set_xlabel('Cluster 5')
#ax7.spines['top'].set_color('m')
#ax7.spines['bottom'].set_color('m')
#ax7.spines['left'].set_color('m')
#ax7.spines['right'].set_color('m')
#ax7.set_xticks([])
#ax7.set_yticks([])
##plt.savefig('VGT106kButterGfpTMaps.png')
#plt.show()
# 
# 
 
 
#def RMSE(u,v):
#    ddd = u-v
#    rmse= np.std(ddd)
#    return rmse  

#VGTrmse = np.zeros((nSamples, ))
#for i in range(nSamples):
#    VGTrmse[i]= RMSE(BaselineData0, Data[:,i])
#
#plt.plot(msSamples,sGFP[:,0])
#plt.plot(msSamples,VGTrmse)
#plt.xlim([msSamples[0], msSamples[-1]])
##plt.title('VGT 106')
#plt.legend(('gfp','rmse'), loc='best')
#plt.savefig('VGT106GfpRmse.png')
#plt.show()


 #colorcache = {u'blue': (0.0, 0.0, 1.0),
#'green': (0.0, 0.5019607843137255, 0.0), 
#'purple': (0.5019607843137255, 0.0, 0.5019607843137255),
# 'r': (1.0, 0.0, 0.0),   
#'firebrick': (0.6980392156862745, 0.13333333333333333, 0.13333333333333333), 
#'cyan': (0.0, 1.0, 1.0),
#'yellow': (1.0, 1.0, 0.0),
#u'c': (0.0, 0.75, 0.75), 
#u'k': (0.0, 0.0, 0.0),
#'y': (0.75, 0.75, 0),
#'darkgoldenrod': (0.7215686274509804, 0.5254901960784314, 0.043137254901960784),
#'magenta': (1.0, 0.0, 1.0), 
#u'0.75': (0.75, 0.75, 0.75), 
#'0.8': (0.8, 0.8, 0.8), 
#u'w': (1.0, 1.0, 1.0),
#u'#bcbcbc': (0.7372549019607844, 0.7372549019607844, 0.7372549019607844), 
#u'white': (1.0, 1.0, 1.0), 
#u'#ffed6f': (1.0, 0.9294117647058824, 0.43529411764705883), 
#u'#467821': (0.27450980392156865, 0.47058823529411764, 0.12941176470588237), 
#u'#eeeeee': (0.9333333333333333, 0.9333333333333333, 0.9333333333333333), 
#u'#F0E442': (0.9411764705882353, 0.8941176470588236, 0.25882352941176473), 
#u'0.50': (0.5, 0.5, 0.5), 
#u'#E24A33': (0.8862745098039215, 0.2901960784313726, 0.2), 
#u'#f0f0f0': (0.9411764705882353, 0.9411764705882353, 0.9411764705882353), 
#u'0.40': (0.4, 0.4, 0.4), 
#'#afeeee': (0.6862745098039216, 0.9333333333333333, 0.9333333333333333),  
#'0.5': (0.5, 0.5, 0.5), 
#u'#fc4f30': (0.9882352941176471, 0.30980392156862746, 0.18823529411764706),  
#u'0.00': (0.0, 0.0, 0.0), 
#u'#bfbbd9': (0.7490196078431373, 0.7333333333333333, 0.8509803921568627), 
#u'#ccebc4': (0.8, 0.9215686274509803, 0.7686274509803922), 
#u'#A60628': (0.6509803921568628, 0.023529411764705882, 0.1568627450980392), 
#u'#988ED5': (0.596078431372549, 0.5568627450980392, 0.8352941176470589), 
#u'#777777': (0.4666666666666667, 0.4666666666666667, 0.4666666666666667), 
#u'#EEEEEE': (0.9333333333333333, 0.9333333333333333, 0.9333333333333333), 
#u'#fdb462': (0.9921568627450981, 0.7058823529411765, 0.3843137254901961), 
# u'#FFB5B8': (1.0, 0.7098039215686275, 0.7215686274509804), 
# u'#30a2da': (0.18823529411764706, 0.6352941176470588, 0.8549019607843137), 
# u'#555555': (0.3333333333333333, 0.3333333333333333, 0.3333333333333333), 
# u'#7A68A6': (0.47843137254901963, 0.40784313725490196, 0.6509803921568628), 
# u'#8b8b8b': (0.5450980392156862, 0.5450980392156862, 0.5450980392156862), 
# u'gray': (0.5019607843137255, 0.5019607843137255, 0.5019607843137255), 
# u'#8dd3c7': (0.5529411764705883, 0.8274509803921568, 0.7803921568627451), 
# u'#bc82bd': (0.7372549019607844, 0.5098039215686274, 0.7411764705882353), 
# u'#CC79A7': (0.8, 0.4745098039215686, 0.6549019607843137),  
# u'#E5E5E5': (0.8980392156862745, 0.8980392156862745, 0.8980392156862745), 
# u'0.70': (0.7, 0.7, 0.7), 
# u'#009E73': (0.0, 0.6196078431372549, 0.45098039215686275), 
# u'#FBC15E': (0.984313725490196, 0.7568627450980392, 0.3686274509803922), 
# u'#feffb3': (0.996078431372549, 1.0, 0.7019607843137254), 
# u'#56B4E9': (0.33725490196078434, 0.7058823529411765, 0.9137254901960784), 
# u'#e5ae38': (0.8980392156862745, 0.6823529411764706, 0.2196078431372549), 
# u'#348ABD': (0.20392156862745098, 0.5411764705882353, 0.7411764705882353), 
# u'#cbcbcb': (0.796078431372549, 0.796078431372549, 0.796078431372549),  
# u'#D55E00': (0.8352941176470589, 0.3686274509803922, 0.0), 
# u'#81b1d2': (0.5058823529411764, 0.6941176470588235, 0.8235294117647058),
# u'#8EBA42': (0.5568627450980392, 0.7294117647058823, 0.25882352941176473), 
# u'#0072B2': (0.0, 0.4470588235294118, 0.6980392156862745), 
# u'#6d904f': (0.42745098039215684, 0.5647058823529412, 0.30980392156862746),   
# '#00FFCC': (0.0, 1.0, 0.8), 
# u'#fa8174': (0.9803921568627451, 0.5058823529411764, 0.4549019607843137), 
# u'#b3de69': (0.7019607843137254, 0.8705882352941177, 0.4117647058823529)}
#colorkeys = list(colorcache.keys())

