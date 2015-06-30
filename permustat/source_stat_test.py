
import numpy as np   
import matplotlib.pyplot as plt 
from mpl_toolkits.mplot3d import Axes3D
from scipy.stats import norm #, t
import matplotlib.mlab as mlab
from swcm import Imatrix 
from mff import read_mff_header, read_mff_data  #, getEpochInfos, mff_getSummaryInfo
from pandas import DataFrame 
import pandas as pd 


filePath ='/Users/jesong1126/Work/Data/VGT/VGT_105_fil_segA_mff_32_100_32_tp_bcr_blc.mff'
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
#epochLabels = hdr['orig']['epochLabels'] 
#[epochLabel.encode('utf-8') for epochLabel in epochLabels]  
SegStatus = hdr['orig']['epochSegStatus']  
GoodSeg = []
epochLabels = []
for i in range(len(SegStatus)):
    if SegStatus[i] == 'good' :  
        GoodSeg.append(i)
        epochLabels.append(trialsName[i])
        
#GoodSeg = np.array(GoodSeg, 'i4')    
data = data[:nE,:,GoodSeg]   
nTrials = len(GoodSeg)

## average reference 
#H = np.identity(nC) - np.ones((nC, nC))/nC 
H = np.identity(nE) - np.ones((nE, nE))/nE 
s = np.zeros(data.shape) 
sGFP = np.zeros((nSamples, nTrials))
if nTrials > 1 :
    for i in range(nTrials):
        s[:,:,i] = np.dot(H, data[:,:,i])
        sGFP[:,i] = np.std(s[:,:,i], axis=0)

erp = np.mean(s, axis=2)
sGFPmean = np.mean(sGFP, axis=1)
#s = s[:, :1000, :] #sGFP = sGFP[:1000, :] #nSamples = 1000 
baseline = (nSamplesPre * 1000 / srate)  
msSamples = np.arange(0, nSamples,1) * 1000/srate  - baseline 
xlimMin = msSamples[0]
xlimMax = msSamples[-1]

fig, ax = plt.subplots(4,1)
ax[0].plot(msSamples, erp.T)
ax[1].plot(msSamples, np.std(erp, axis=0))
ax[2].plot(msSamples, sGFP)
ax[3].plot(msSamples, sGFPmean)
plt.show()
#fig.savefig('EEG105.png')    

gsn257 = pd.read_table('/Users/jesong1126/Python27/GeoPy/nscolor_hgsn/GSN257ToPy.sfp') 
frame = DataFrame(gsn257, columns=['ChLabel', 'X3', 'Y3', 'Z3','X2','Y2'])
x2 = frame.values[:nE,4] 
y2 = frame.values[:nE,5] 

#for tt in [100, 300, 800, 1100, 1500]:
#    for i in range(nTrials): 
#        fig = plt.figure()
#        pltscatter = plt.scatter(x2, y2, c=s[:,tt,i], s=50, cmap=plt.get_cmap('seismic') )
#        plt.colorbar(pltscatter)#cax, ticks=[-1, 0, 1])
##        fig.savefig('Topo105_%d_%dms.png'%(i,tt))    
#        plt.close()
#
#for i in range(nTrials): 
#    fig = plt.figure()
#    pltscatter = plt.imshow(s[:,:,i], cmap='seismic', aspect="auto", extent=[-100,1500, 257,0])  
#    plt.colorbar(pltscatter) 
##    fig.savefig('Topo105_%d.png'%(i))    
#    plt.close()

for tt in [100, 300, 800, 1100, 1500]:
    fig = plt.figure()
    pltscatter = plt.scatter(x2, y2, c=erp[:, tt], s=50, cmap=plt.get_cmap('seismic') )
    plt.colorbar(pltscatter) 
#    fig.savefig('Topo105_erp_%dms.png'%(tt))    
    plt.close()

##-----------------------------------------------------------------------------
eegMean = np.mean(s, axis=2)
eegStd = np.std(s, axis=2)
eegT = eegMean / eegStd  

[np.min(eegMean), np.max(eegMean)]
[np.min(eegStd), np.max(eegStd)]
[np.min(eegT), np.max(eegT)]

#
ut = 2.58 * np.std(eegT)
#ut = 3 * np.std(eegT)
eee = (abs(eegT) > ut) 
nSigChan = np.sum(eee, axis=0) 
nSigTime = np.sum(eee, axis=1)

atTime= np.where(nSigChan>0)[0]
a = np.where(np.diff(atTime)>1)[0] + 1
MicroStart = np.concatenate((np.array(atTime[0]).reshape((1)), atTime[a]))
MicroEnd = np.concatenate((atTime[a-1], np.array(atTime[-1]).reshape((1))))
MicroMiddle = (MicroStart + MicroEnd)/2
len(MicroMiddle) #9 array([ 284,  538,  610, 1500 ])

fig = plt.figure()
im = plt.imshow(eegMean, cmap='seismic', aspect="auto", extent=[-100,1500, 257,0]) 
plt.colorbar(im)
#fig.savefig('Topo105_mean_image.png')    
plt.show()    

fig = plt.figure()
im = plt.imshow(eegStd, cmap='hot', aspect="auto", extent=[-100,1500, 257,0]) 
plt.colorbar(im) 
#fig.savefig('Topo105_std_image.png')    
plt.show()

fig = plt.figure()
im = plt.imshow(eegT, cmap='seismic', aspect="auto", extent=[-100,1500, 257,0]) 
plt.colorbar(im) 
#fig.savefig('Topo105_SN_image.png')    
plt.show()
 
fig = plt.figure()
im = plt.imshow(eee, cmap='hot', aspect="auto", extent=[-100,1500, 257,0]) 
cbar = plt.colorbar(im) 
cbar.set_ticks([0,1])
cbar.set_ticklabels(['No', 'Yes']) 
#fig.savefig('Topo105_SigVT_image.png')    
plt.show()


eegTvector = eegT.reshape((nE*nSamples,1))
mu = np.mean(eegTvector)
sigma = np.std(eegTvector)
vu = mu + 2.58 * sigma 
vl = mu - 2.58 * sigma 

fig = plt.figure()
num_bins = 50
n, bins, patches = plt.hist(eegTvector, num_bins, normed=1, facecolor='green', alpha=0.5)
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, 'r--')
plt.vlines([vl, vu], 0, 0.1, color='black', lw=2)
plt.xlim([-3, 3])
plt.xlabel('T')
plt.ylabel('Probability')
plt.title(r'Normal: $\mu=%5.2f$, $\sigma=%3.2f$'%(mu, sigma)) 
#fig.savefig('Topo105_NormalHist.png')    
plt.show()

fig, ax = plt.subplots(4,2,figsize=(20,10))
ax[0,0].plot(msSamples, eegMean.T)
ax[0,0].set_xlim([msSamples[0], msSamples[-1]])
ax[0,0].set_ylabel('mean')
ax[0,0].set_title('times')
ax[1,0].plot(msSamples, eegStd.T)
ax[1,0].set_xlim([msSamples[0], msSamples[-1]])
ax[1,0].set_ylabel('sdv')
ax[2,0].plot(msSamples, eegT.T)
ax[2,0].set_xlim([msSamples[0], msSamples[-1]])
ax[2,0].set_ylabel('Normal')
ax[2,0].axhline(vu, color='black', lw=2)
ax[2,0].axhline(vl, color='black', lw=2)
ax[3,0].scatter(msSamples,nSigChan, s=1)
ax[3,0].set_xlim([msSamples[0], msSamples[-1]])
ax[3,0].set_ylabel('Num Sig.')
ax[0,1].plot(eegMean)
ax[0,1].set_xlim([0, nC])
ax[0,1].set_ylabel('mean')
ax[0,1].set_title('channels')
ax[1,1].plot(eegStd)
ax[1,1].set_xlim([0, nC])
ax[1,1].set_ylabel('sdv')
ax[2,1].plot(eegT)
ax[2,1].set_xlim([0, nC])
ax[2,1].set_ylabel('Normal')
ax[2,1].axhline(vu, color='black', lw=2)
ax[2,1].axhline(vl, color='black', lw=2)
ax[3,1].scatter(range(nE), nSigTime, s=1)
ax[3,1].set_xlim([0, nC])
ax[3,1].set_ylabel('Num Sig.')
#fig.savefig('EEG105_wave.png')    
plt.show()


fig = plt.figure()
plt.scatter(msSamples,nSigChan, s=1)
plt.xlim([msSamples[0], msSamples[-1]])
plt.ylim([-2, np.max(nSigChan)+2])
#plt.hlines(cutNv, msSamples[0], msSamples[-1], color='red')  
plt.vlines(MicroMiddle-100, 0, nSigChan[MicroMiddle], color='red')  
plt.ylabel('Num Signif.')
#fig.savefig('EEG105_SigTimes.png')    
plt.show()

for tt in MicroMiddle: 
    fig = plt.figure()
    pltscatter = plt.scatter(x2, y2, c=eegT[:, tt], s=100, cmap=plt.get_cmap('seismic') )
    plt.colorbar(pltscatter,  shrink=0.9, pad= 0.01) 
    plt.xlim([-1.6, 1.6])
    plt.ylim([-1.3, 1.2])
    plt.axis('off')
#    fig.savefig('Topo105_eegT_%dms.png'%(tt-100))    
    plt.show()


for tt in MicroMiddle: 
    fig = plt.figure()
    pltscatter = plt.scatter(x2, y2, c=-eee[:, tt], s=200, cmap=plt.get_cmap('hot') )
    plt.xlim([-1.6, 1.6])
    plt.ylim([-1.3, 1.2])
    plt.axis('off') 
    cbar = plt.colorbar(pltscatter,  shrink=0.9, pad= 0.01) 
    cbar.set_ticks([0,1])
    cbar.set_ticklabels(['Yes', 'No']) 
#    plt.title('%dms'%(tt-100))
#    fig.savefig('Topo105_sig_eeg_%dms.png'%(tt-100))    
    plt.show()

##-----------------------------------------------------------------------------
#nC = 257  
#nE = nC - 1 

lfmfilename = '/Users/jesong1126/Work/Data/Ryan/223_RJ/IHM_32yr/223_IHM_LFM/223_IHM_LFM_output/223_IHM_LFM_ori.lfm' #2400
#lfmfilename = '/Users/jesong1126/Work/Data/Ryan/223_RJ/CAHM_32yr/223_32yr_LFM/223_32yr_LFM_output/223_32yr_LFM_trip.lfm' #1818, 5454
#lfmfilename = '/Users/jesong1126/Work/Data/Ryan/223_RJ/IHM_32yr/223_IHM_LFM/223_IHM_LFM_output/223_IHM_LFM_trip.lfm' #7200 
#lfmfilename = '/Users/jesong1126/Work/Data/Ryan/223_RJ/CAHM_40yr/223_40yr_LFM/223_40yr_LFM_output/223_40yr_LFM_trip.lfm' #1612, 4836
fd = open(lfmfilename , 'r')
data = np.fromfile(file=fd, dtype=np.dtype('float'))  
fd.close()
nD = data.shape[0] / nE
K = data.reshape((nE, nD))  

nV = nD 
nV = nD / 3 #1818 #2447 #1612 #2400  

fig = plt.figure() 
plt.plot(K[:,:3])
plt.show()    

#K = data.reshape((nD, nE)) 
#K = K.T 

#fig = plt.figure() 
#plt.imshow(K, cmap='seismic', aspect="auto")
#plt.show()    

#VoxCoord = np.loadtxt('/Users/jesong1126/Documents/MATLAB/SWCM/data/TalCoord2447.txt', delimiter='\t', dtype=np.dtype('i4') )
from3To1 = np.kron(np.eye(nV), np.ones((1, 3)))

alpha = 0.0001
Imn = Imatrix.MN(alpha, K) 
sdeneegMean = np.dot(Imn, eegMean)   
sdeneegT = np.dot(Imn, eegT)   

fig, ax = plt.subplots(2,2) 
im = ax[0,0].imshow(sdeneegMean, cmap='seismic', aspect="auto")
#ax[0,0].colorbar(im)
ax[0,1].plot(msSamples, sdeneegMean.T)
im = ax[0,0].imshow(sdeneegT, cmap='seismic', aspect="auto")
#ax[0,0].colorbar(im)
ax[0,1].plot(msSamples, sdeneegT.T)
plt.show()    

#Imn = Imatrix.LORETA(alpha, K) 
temp = np.kron(np.eye(nV), np.ones((1, 3)))

#rms = np.zeros((nV,nSamples,nTrials))
sden = np.zeros((nD, nSamples, nTrials))
for ii in range(nTrials):
    sinv = np.dot(Imn, s[:,:, ii])
    sden[:,:, ii] = sinv
#    sinv2 = np.multiply(sinv, sinv) 
#    rms[:,:,ii] = np.sqrt(np.dot(from3To1, sinv2) / 3) 

       
sdenMean = np.mean(sden, axis=2)
sdenSigma = np.std(sden, axis=2)
sdenT =  sdenMean / sdenSigma # * np.sqrt(44) 

[np.min(sdenMean), np.max(sdenMean)]
[np.min(sdenSigma), np.max(sdenSigma)]
[np.min(sdenT), np.max(sdenT)]

ut = 2.58 * np.std(sdenT)
sss = (abs(sdenT) > ut) 
 
fig = plt.figure()
im = plt.imshow(sdenMean, cmap='seismic', aspect="auto")
plt.colorbar(im)
#fig.savefig('Source105_mean_image.png')    
plt.show()    

fig = plt.figure()
im = plt.imshow(sdenSigma, cmap='hot', aspect="auto")
plt.colorbar(im)
plt.show()    
#fig.savefig('Source105_std_image.png')    

fig = plt.figure()
im = plt.imshow(sdenT, cmap='seismic', aspect="auto")
plt.colorbar(im)
plt.show()    
#fig.savefig('Source105_sn_image.png')    

fig = plt.figure()
im = plt.imshow(sss, cmap='hot', aspect="auto", extent=[-100,1500, 7341,0]) 
cbar = plt.colorbar(im) 
cbar.set_ticks([0,1])
cbar.set_ticklabels(['No', 'Yes']) 
#fig.savefig('Source105_sig_image.png')    
plt.show()

sdenTvector = sdenT.reshape((nD*1600,1))
mu = np.mean(sdenTvector)
sigma = np.std(sdenTvector)
vu = mu + 2.58 * sigma 
vl = mu - 2.58 * sigma 

fig = plt.figure()
num_bins = 50
n, bins, patches = plt.hist(sdenTvector, num_bins, normed=1, facecolor='green', alpha=0.5)
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, 'r--')
plt.vlines([vl, vu], 0, 0.1, color='black', lw=2)
plt.xlim([-2, 2])
plt.xlabel('T')
plt.ylabel('Probability')
plt.title(r'Normal: $\mu=%5.2f$, $\sigma=%3.2f$'%(mu, sigma)) 
#fig.savefig('Source105_NormalHist.png')    
plt.show()

##-----------------------------------------------------------------------------
rmsMean = np.mean(rms, axis=2)
rmsSigma = np.std(rms, axis=2)
rmsT =  rmsMean / rmsSigma # * np.sqrt(44) 

[np.min(rmsMean), np.max(rmsMean)]
[np.min(rmsSigma), np.max(rmsSigma)]
[np.min(rmsT), np.max(rmsT)]

rmsTvector = rmsT.reshape((2447*1600,1))
mu = np.mean(rmsTvector)
sigma = np.std(rmsTvector)
vu = mu + 2.58 * sigma 
vl = mu - 2.58 * sigma 

rrr = (rmsT > vu)  | (rmsT < vl)
 
fig = plt.figure()
im = plt.imshow(rmsMean, cmap='hot', aspect="auto")
plt.colorbar(im)
#fig.savefig('Source105_mean_image.png')    
plt.show()    

fig = plt.figure()
im = plt.imshow(rmsSigma, cmap='hot', aspect="auto")
plt.colorbar(im)
plt.show()    
#fig.savefig('Source105_std_image.png')    

fig = plt.figure()
im = plt.imshow(rmsT, cmap='seismic', aspect="auto")
plt.colorbar(im)
plt.show()    
#fig.savefig('Source105_sn_image.png')    

fig = plt.figure()
im = plt.imshow(rrr, cmap='hot', aspect="auto", extent=[-100,1500, 2447,0]) 
cbar = plt.colorbar(im) 
cbar.set_ticks([0,1])
cbar.set_ticklabels(['No', 'Yes']) 
#fig.savefig('Source105_sig_image.png')    
plt.show()


fig = plt.figure()
num_bins = 50
n, bins, patches = plt.hist(rmsTvector, num_bins, normed=1, facecolor='green', alpha=0.5)
y = mlab.normpdf(bins, mu, sigma)
plt.plot(bins, y, 'r--')
plt.vlines([vl, vu], 0, 0.1, color='black', lw=2)
plt.xlim([0, np.max(rmsT)])
plt.xlabel('T')
plt.ylabel('Probability')
plt.title(r'Normal: $\mu=%5.2f$, $\sigma=%3.2f$'%(mu, sigma)) 
#fig.savefig('Source105_NormalHist.png')    
plt.show()

##-----------------------------------------------------------------------------

RHemi = np.where(VoxCoord[:,1] >= 0)[0]
LHemi = np.where(VoxCoord[:,1] <0)[0] 
np.sum(VoxCoord[:,1] >= 0) 

nSigR = np.sum(rrr[RHemi,:], axis=0)
nSigL = np.sum(rrr[LHemi,:], axis=0)  

nSig = np.sum(rrr, axis=0)
nSigT = np.sum(rrr, axis=1)

rangeNV = np.arange(nV) 
plt.scatter(rangeNV[RHemi], nSigT[RHemi], c='red')
plt.scatter(rangeNV[LHemi], nSigT[LHemi], c='blue')
plt.show()

fig, ax = plt.subplots(4,2,figsize=(20,10))
ax[0,0].plot(msSamples, rmsMean.T)
ax[0,0].set_xlim([msSamples[0], msSamples[-1]])
ax[0,0].set_ylabel('mean')
ax[0,0].set_title('times')
ax[1,0].plot(msSamples, rmsSigma.T)
ax[1,0].set_xlim([msSamples[0], msSamples[-1]])
ax[1,0].set_ylabel('sdv')
ax[2,0].plot(msSamples, rmsT.T)
ax[2,0].set_xlim([msSamples[0], msSamples[-1]])
ax[2,0].set_ylabel('Normal')
ax[2,0].axhline(vu, color='black', lw=2)
ax[2,0].axhline(vl, color='black', lw=2)
ax[3,0].plot(msSamples, nSigL,  c='blue', label= 'left')
ax[3,0].plot(msSamples, nSigR,  c='red', label= 'right') 
ax[3,0].set_ylim([-5, max(max(nSigL), max(nSigR))+5])
ax[3,0].set_ylim([-5, max(max(nSigL), max(nSigR))+5])
ax[3,0].set_xlim([msSamples[0], msSamples[-1]])
ax[3,0].set_ylabel('Num Sig.')
ax[0,1].plot(rmsMean)
ax[0,1].set_xlim([0, nV])
ax[0,1].set_ylabel('mean')
ax[0,1].set_title('dipoles')
ax[1,1].plot(rmsSigma)
ax[1,1].set_xlim([0, nV])
ax[1,1].set_ylabel('sdv')
ax[2,1].plot(rmsT)
ax[2,1].set_xlim([0, nV])
ax[2,1].set_ylabel('Normal')
ax[2,1].axhline(vu, color='black', lw=2)
ax[2,1].axhline(vl, color='black', lw=2)
ax[3,1].scatter(rangeNV[LHemi], nSigT[LHemi], c='blue', label= 'left') #, s=1)
ax[3,1].scatter(rangeNV[RHemi], nSigT[RHemi], c='red', label= 'right')
ax[3,1].set_ylabel('Num Sig.')
ax[3,1].set_xlim([0, nV])
#fig.savefig('Rms105_wave.png')    
plt.show()

fig = plt.figure()
plt.plot(msSamples, rmsSigma.T)
plt.xlim([msSamples[0], msSamples[-1]])
plt.ylabel('sdv')
plt.show()

fig = plt.figure()
plt.plot(msSamples, rmsT.T)
plt.xlim([msSamples[0], msSamples[-1]])
plt.ylabel('Normal')
plt.hlines(vu, msSamples[0], msSamples[-1], color='black', lw=2)
plt.hlines(vl, msSamples[0], msSamples[-1], color='black', lw=2)
plt.show()

#
#sigma800 = rmsSigma[:, 900:1100] 
#
#fig = plt.figure()
#plt.plot(sigma800.T)
##plt.xlim([msSamples[0], msSamples[-1]])
#plt.ylabel('sdv [800, 1000] ')
#plt.show()
#
#sigvox800 = np.where(sigma800 > 10)[0]
#voxId800 = np.unique(sigvox800)/3 
##voxId800= np.array(voxId800, 'i4')
#LhemiVoxId800 = np.intersect1d(LHemi ,voxId800)
#RhemiVoxId800 = np.intersect1d(RHemi ,voxId800)
#
#fig, ax = plt.subplots(2,2)#,figsize=(40,30))
#ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c='0.7', s=10)  
#ax[0,0].scatter(VoxCoord[voxId800,1], VoxCoord[voxId800,2], c='red', s=100 )
#ax[0,0].set_axis_off()
#ax[0,0].set_aspect('equal')    
#ax[0,0].set_xlim([-70, 70])
#ax[0,0].set_ylim([-106, 72])
#     
#ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c='0.7', s=10)  
#ax[0,1].scatter(VoxCoord[voxId800,1], VoxCoord[voxId800,3], c='red', s=100 )
#ax[0,1].set_aspect('equal')
#ax[0,1].set_axis_off()
#ax[0,1].set_xlim([-70, 70])
#ax[0,1].set_ylim([-45,82])
#    
#ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3], c='0.7', s=10)  
#ax[1,0].scatter(-VoxCoord[LhemiVoxId800,2], VoxCoord[LhemiVoxId800,3], c='red', s=100 )
#ax[1,0].set_aspect('equal')
#ax[1,0].set_axis_off()    
#ax[1,0].set_xlim([-70, 106])
#ax[1,0].set_ylim([-45, 82])
#
#ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3], c='0.7', s=10)  
#ax[1,1].scatter(VoxCoord[RhemiVoxId800,2], VoxCoord[RhemiVoxId800,3], c='red', s=100 )
#ax[1,1].set_aspect('equal')
#ax[1,1].set_axis_off()
#ax[1,1].set_xlim([-106, 70])
#ax[1,1].set_ylim([-45, 82]) 
##fig.savefig('Source105_Dist_scatter_%dms.png'%(MicroMiddle[i]-100))    
##    plt.close() 
##fig.savefig('Source105_Sigma_800_1000_Voxels.png')    
#plt.show()    


fig = plt.figure()
plt.scatter(msSamples, nSigL,  c='blue')
plt.scatter(msSamples, nSigR,  c='red') 
plt.ylim([-5, max(max(nSigL), max(nSigR))+5])
plt.xlim([msSamples[0], msSamples[-1]])
plt.vlines(MicroMiddle-100, 0, nSigL[MicroMiddle], color='black')  
plt.ylabel('Num Sig.')
#fig.savefig('Source105_SigTimes.png')    
plt.show()


# MicroMiddle
for i in range(4):
    i = 2
    fig, ax = plt.subplots(2,2) 
    nSigVoxGrpI= rmsT[:, MicroMiddle[i]]  
    ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,0].set_axis_off()
    ax[0,0].set_aspect('equal')    
    ax[0,0].set_xlim([-70, 70])
    ax[0,0].set_ylim([-106, 72])
    
    ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,1].set_aspect('equal')
    ax[0,1].set_axis_off()
    ax[0,1].set_xlim([-70, 70])
    ax[0,1].set_ylim([-45,82])
    
    ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c=nSigVoxGrpI[LHemi], s=100, cmap='hot')    
    ax[1,0].set_aspect('equal')
    ax[1,0].set_axis_off()    
    ax[1,0].set_xlim([-70, 106])
    ax[1,0].set_ylim([-45, 82])

    ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c=nSigVoxGrpI[RHemi], s=100, cmap='hot')   
    ax[1,1].set_aspect('equal')
    ax[1,1].set_axis_off()
    ax[1,1].set_xlim([-106, 70])
    ax[1,1].set_ylim([-45, 82])
#    fig.savefig('Source105_Dist_scatter_%dms.png'%(MicroMiddle[i]-100))    
#    plt.close() 
#fig.savefig('Source105_Dist_scatter.png')    
plt.show()    

fig, ax = plt.subplots() 
ax.scatter(VoxCoord[:,1], VoxCoord[:,2], c='gray', s=100, cmap='hot')  
ax.set_aspect('equal')
ax.set_xlim([-70,70])
ax.set_ylim([-106, 70])
ax.set_xticks([-65, 65])
ax.set_xticklabels(['left','right'])
ax.set_yticks([-98, 62])
ax.set_yticklabels(['back','front'], rotation=90)
#fig.savefig('axial.png')    
plt.show()    

fig, ax = plt.subplots() 
ax.scatter(VoxCoord[:,1], VoxCoord[:,3], c='gray', s=100, cmap='hot')  
ax.set_aspect('equal')
ax.set_xlim([-70, 70])
ax.set_ylim([-45, 82])
ax.set_xticks([-65, 65])
ax.set_xticklabels(['left','right'])
ax.set_yticks([-35, 78])
ax.set_yticklabels(['bottom','top'], rotation=90)
#fig.savefig('coronal.png')    
plt.show()    

fig, ax = plt.subplots() 
ax.scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c='gray', s=100, cmap='hot')    
ax.set_aspect('equal')
ax.set_xlim([-70, 106])
ax.set_ylim([-45, 82])
ax.set_xticks([-65, 100])
ax.set_xticklabels(['front','back'])
ax.set_yticks([-35, 76])
ax.set_yticklabels(['bottom','top'], rotation=90)
#fig.savefig('sagitalL.png')    
plt.show()    

fig, ax = plt.subplots() 
ax.scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c='gray', s=100, cmap='hot')   
ax.set_aspect('equal')
ax.set_xlim([-106, 70])
ax.set_ylim([-45, 82])
ax.set_xticks([-100, 65])
ax.set_xticklabels(['back', 'front'])
ax.set_yticks([-35, 76])
ax.set_yticklabels(['bottom','top'], rotation=90)
#fig.savefig('sagitalR.png')    
plt.show()    


#a = np.where(np.diff(vlinesS)>1)[0] + 1
#MicroStart = np.concatenate((np.array(vlinesS[0]).reshape((1)), vlinesS[a]))
#MicroEnd = np.concatenate((vlinesS[a-1], np.array(vlinesS[-1]).reshape((1))))
#MicroMiddle = (MicroStart + MicroEnd)/2
#len(MicroMiddle) #9 array([ 285,  575,  617,  656,  748,  752, 1362, 1434, 1540])
MicroMiddle = [284, 538, 610, 1500] 
LatIdx = []
for i in range(len(MicroMiddle)): 
    SigVox = 1. * rrr[:, MicroMiddle[i]]
    LatIdx.append([MicroMiddle[i]-100, np.sum(SigVox[LHemi])/np.sum(VoxCoord[:,1]<0), 
                   np.sum(SigVox[RHemi])/np.sum((VoxCoord[:,1] >= 0))])

for i in range(len(MicroMiddle)):
    fig, ax = plt.subplots(2,2)  
#    nSigVoxGrpI = 1. * rmsT[:, MicroMiddle[i]]  
    nSigVoxGrpI = rrr[:, MicroMiddle[i]]  
    ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,0].set_axis_off()
    ax[0,0].set_aspect('equal')    
    ax[0,0].set_xlim([-70, 70])
    ax[0,0].set_ylim([-106, 72])
    
    ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,1].set_aspect('equal')
    ax[0,1].set_axis_off()
    ax[0,1].set_xlim([-70, 70])
    ax[0,1].set_ylim([-45,82])
    
    ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c=nSigVoxGrpI[LHemi], s=100, cmap='hot')    
    ax[1,0].set_aspect('equal')
    ax[1,0].set_axis_off()    
    ax[1,0].set_xlim([-70, 106])
    ax[1,0].set_ylim([-45, 82])
    ax[1,0].set_title('Left %3.2f ' %(LatIdx[i][1]))

    ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c=nSigVoxGrpI[RHemi], s=100, cmap='hot')   
    ax[1,1].set_aspect('equal')
    ax[1,1].set_axis_off()
    ax[1,1].set_xlim([-106, 70])
    ax[1,1].set_ylim([-45, 82])
    ax[1,1].set_title('Right %3.2f' %(LatIdx[i][2]))
#    fig.savefig('LatIdx105_%dms.png' %(MicroMiddle[i]-100))    

plt.show()    


for i in range(40):
    fig, ax = plt.subplots(2,2) #,figsize=(20,10))
    nSigVoxGrpI = rrr[:, i*40]  
    ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,0].set_axis_off()
    ax[0,0].set_aspect('equal')    
    ax[0,0].set_xlim([-70, 70])
    ax[0,0].set_ylim([-106, 72])
    ax[0,0].set_title('%d ms ' %(i*40-100))
    
    ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,1].set_aspect('equal')
    ax[0,1].set_axis_off()
    ax[0,1].set_xlim([-70, 70])
    ax[0,1].set_ylim([-45,82])
    
    ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c=nSigVoxGrpI[LHemi], s=100, cmap='hot')    
    ax[1,0].set_aspect('equal')
    ax[1,0].set_axis_off()    
    ax[1,0].set_xlim([-70, 106])
    ax[1,0].set_ylim([-45, 82])

    ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c=nSigVoxGrpI[RHemi], s=100, cmap='hot')   
    ax[1,1].set_aspect('equal')
    ax[1,1].set_axis_off()
    ax[1,1].set_xlim([-106, 70])
    ax[1,1].set_ylim([-45, 82])
    fig.savefig('/Users/jesong1126/Work/VGT/tex/Figure11142014/RmsLIMovie/%04d.png' %(i*40))         
    plt.close()
  
#ffmpeg -framerate 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p RmsLI105.mp4


for i in range(40):
    fig, ax = plt.subplots(2,2) 
    nSigVoxGrpI = rmsT[:, i*40]  
    ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,0].set_axis_off()
    ax[0,0].set_aspect('equal')    
    ax[0,0].set_xlim([-70, 70])
    ax[0,0].set_ylim([-106, 72])
    ax[0,0].set_title('%d ms ' %(i*40-100))
    
    ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,1].set_aspect('equal')
    ax[0,1].set_axis_off()
    ax[0,1].set_xlim([-70, 70])
    ax[0,1].set_ylim([-45,82])
    
    ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c=nSigVoxGrpI[LHemi], s=100, cmap='hot')    
    ax[1,0].set_aspect('equal')
    ax[1,0].set_axis_off()    
    ax[1,0].set_xlim([-70, 106])
    ax[1,0].set_ylim([-45, 82])

    ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c=nSigVoxGrpI[RHemi], s=100, cmap='hot')   
    ax[1,1].set_aspect('equal')
    ax[1,1].set_axis_off()
    ax[1,1].set_xlim([-106, 70])
    ax[1,1].set_ylim([-45, 82])
    fig.savefig('/Users/jesong1126/Work/VGT/tex/Figure11142014/RmsDistMovie/%04d.png' %(i*40))    
    plt.close()
    
#ffmpeg -framerate 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p RmsDist105.mp4




##-----------------------------------------------------------------------------

RHemi = np.where(VoxCoord[:,1] >= 0)[0]
LHemi = np.where(VoxCoord[:,1] <0)[0] 
np.sum(VoxCoord[:,1] >= 0) 

sdenMean = np.mean(sden, axis=2)
sdenSigma = np.std(sden, axis=2)
sdenT =  sdenMean / sdenSigma # * np.sqrt(44) 

[np.min(sdenMean), np.max(sdenMean)]
[np.min(sdenSigma), np.max(sdenSigma)]
[np.min(sdenT), np.max(sdenT)]

sdenTvector = sdenT.reshape((nD*nSamples,1))
mu = np.mean(sdenTvector)
sigma = np.std(sdenTvector)
vu = mu + 2.58 * sigma 
vl = mu - 2.58 * sigma 

sss = (sdenT > vu)  | (sdenT < vl)

rr = np.kron(3*RHemi, np.ones(3)) + np.kron(np.ones(RHemi.shape[0]),  np.arange(3)) 
rr = np.array(rr, dtype='i4')
nSigR = np.sum(sss[rr,:], axis=0)

ll = np.kron(3*LHemi, np.ones(3)) + np.kron(np.ones(LHemi.shape[0]),  np.arange(3)) 
ll = np.array(ll, dtype='i4')
nSigL = np.sum(sss[ll,:], axis=0)  

nSig = np.sum(sss, axis=0)
nSigT = np.sum(sss, axis=1)
nSigT3 = np.dot(from3To1, nSigT) 

rangeND = np.arange(nD) 
plt.scatter(rangeND[rr], nSigT[rr], c='red')
plt.scatter(rangeND[ll], nSigT[ll], c='blue')
plt.show()

sdenTx = sdenT[range(0, nD, 3),:] 
sdenTy = sdenT[range(1, nD, 3),:] 
sdenTz = sdenT[range(2, nD, 3),:] 
fig, ax = plt.subplots(2,2,figsize=(20,10))
#ax[0,0].plot(msSamples, sdenMean.T)
ax[0,0].plot(msSamples, sdenTx.T, c= 'blue')
ax[0,0].plot(msSamples, sdenTy.T, c= 'red')
ax[0,0].plot(msSamples, sdenTz.T, c= 'green')
ax[0,0].set_xlim([msSamples[0], msSamples[-1]])
ax[0,0].set_ylabel('mean')
ax[0,0].set_title('times')
ax[1,0].plot(msSamples, sdenSigma.T)
ax[1,0].set_xlim([msSamples[0], msSamples[-1]])
ax[1,0].set_ylabel('sdv')
ax[2,0].plot(msSamples, sdenT.T)
ax[2,0].set_xlim([msSamples[0], msSamples[-1]])
ax[2,0].set_ylabel('Normal')
ax[2,0].axhline(vu, color='black', lw=2)
ax[2,0].axhline(vl, color='black', lw=2)
ax[3,0].plot(msSamples, nSigL,  c='blue', label= 'left')
ax[3,0].plot(msSamples, nSigR,  c='red', label= 'right') 
ax[3,0].set_ylim([-5, max(max(nSigL), max(nSigR))+5])
ax[3,0].set_ylim([-5, max(max(nSigL), max(nSigR))+5])
ax[3,0].set_xlim([msSamples[0], msSamples[-1]])
ax[3,0].set_ylabel('Num Sig.')
#ax[0,1].plot(sdenMean)
#ax[0,1].set_xlim([0, nD])
#ax[0,1].set_ylabel('mean')
#ax[0,1].set_title('dipoles')
#ax[1,1].plot(sdenSigma)
#ax[1,1].set_xlim([0, nD])
#ax[1,1].set_ylabel('sdv')
#ax[2,1].plot(sdenT)
#ax[2,1].set_xlim([0, nD])
#ax[2,1].set_ylabel('Normal')
#ax[2,1].axhline(vu, color='black', lw=2)
#ax[2,1].axhline(vl, color='black', lw=2)
#ax[3,1].scatter(rangeND[rr], nSigT[rr], c='red', label= 'right')
#ax[3,1].scatter(rangeND[ll], nSigT[ll], c='blue', label= 'left') #, s=1)
#ax[3,1].set_ylabel('Num Sig.')
#ax[3,1].set_xlim([0, nD])
#fig.savefig('Source105_wave.png')    
plt.show()

fig = plt.figure()
plt.plot(msSamples, sdenSigma1.T)
plt.xlim([msSamples[0], msSamples[-1]])
plt.ylabel('sdv')
plt.show()

sigma800 = sdenSigma1[:, 900:1100] 


fig = plt.figure()
plt.plot(sigma800.T)
#plt.xlim([msSamples[0], msSamples[-1]])
plt.ylabel('sdv [800, 1000] ')
plt.show()

sigvox800 = np.where(sigma800 > 10)[0]
voxId800 = np.unique(sigvox800)/3 
#voxId800= np.array(voxId800, 'i4')
LhemiVoxId800 = np.intersect1d(LHemi ,voxId800)
RhemiVoxId800 = np.intersect1d(RHemi ,voxId800)

fig, ax = plt.subplots(2,2)#,figsize=(40,30))
ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c='0.7', s=10)  
ax[0,0].scatter(VoxCoord[voxId800,1], VoxCoord[voxId800,2], c='red', s=100 )
ax[0,0].set_axis_off()
ax[0,0].set_aspect('equal')    
ax[0,0].set_xlim([-70, 70])
ax[0,0].set_ylim([-106, 72])
     
ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c='0.7', s=10)  
ax[0,1].scatter(VoxCoord[voxId800,1], VoxCoord[voxId800,3], c='red', s=100 )
ax[0,1].set_aspect('equal')
ax[0,1].set_axis_off()
ax[0,1].set_xlim([-70, 70])
ax[0,1].set_ylim([-45,82])
    
ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3], c='0.7', s=10)  
ax[1,0].scatter(-VoxCoord[LhemiVoxId800,2], VoxCoord[LhemiVoxId800,3], c='red', s=100 )
ax[1,0].set_aspect('equal')
ax[1,0].set_axis_off()    
ax[1,0].set_xlim([-70, 106])
ax[1,0].set_ylim([-45, 82])

ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3], c='0.7', s=10)  
ax[1,1].scatter(VoxCoord[RhemiVoxId800,2], VoxCoord[RhemiVoxId800,3], c='red', s=100 )
ax[1,1].set_aspect('equal')
ax[1,1].set_axis_off()
ax[1,1].set_xlim([-106, 70])
ax[1,1].set_ylim([-45, 82]) 
#fig.savefig('Source105_Dist_scatter_%dms.png'%(MicroMiddle[i]-100))    
#    plt.close() 
#fig.savefig('Source105_Sigma_800_1000_Voxels.png')    
plt.show()    


fig = plt.figure()
plt.scatter(msSamples, nSigL,  c='blue')
plt.scatter(msSamples, nSigR,  c='red') 
plt.ylim([-5, max(max(nSigL), max(nSigR))+5])
plt.xlim([msSamples[0], msSamples[-1]])
plt.vlines(MicroMiddle-100, 0, nSigL[MicroMiddle], color='black')  
plt.ylabel('Num Sig.')
#fig.savefig('Source105_SigTimes.png')    
plt.show()


# MicroMiddle
for i in range(4):
    fig, ax = plt.subplots(2,2)#,figsize=(40,30))
    SigVoxGrpI = abs(sdenT1[:, MicroMiddle[i]]) 
    nSigVoxGrpI= np.dot(from3To1, SigVoxGrpI)
    ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,0].set_axis_off()
    ax[0,0].set_aspect('equal')    
    ax[0,0].set_xlim([-70, 70])
    ax[0,0].set_ylim([-106, 72])
    
    ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,1].set_aspect('equal')
    ax[0,1].set_axis_off()
    ax[0,1].set_xlim([-70, 70])
    ax[0,1].set_ylim([-45,82])
    
    ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c=nSigVoxGrpI[LHemi], s=100, cmap='hot')    
    ax[1,0].set_aspect('equal')
    ax[1,0].set_axis_off()    
    ax[1,0].set_xlim([-70, 106])
    ax[1,0].set_ylim([-45, 82])

    ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c=nSigVoxGrpI[RHemi], s=100, cmap='hot')   
    ax[1,1].set_aspect('equal')
    ax[1,1].set_axis_off()
    ax[1,1].set_xlim([-106, 70])
    ax[1,1].set_ylim([-45, 82])
    fig.savefig('Source105_Dist_scatter_%dms.png'%(MicroMiddle[i]-100))    
#    plt.close() 
#fig.savefig('Source105_Dist_scatter.png')    
plt.show()    

fig, ax = plt.subplots()#,figsize=(40,30))
nSigVoxGrpI= np.dot(from3To1, SigVoxGrpI)
ax.scatter(VoxCoord[:,1], VoxCoord[:,2], c='gray', s=100, cmap='hot')  
ax.set_aspect('equal')
ax.set_xlim([-70,70])
ax.set_ylim([-106, 70])
ax.set_xticks([-65, 65])
ax.set_xticklabels(['left','right'])
ax.set_yticks([-98, 62])
ax.set_yticklabels(['back','front'], rotation=90)
#fig.savefig('axial.png')    
plt.show()    

fig, ax = plt.subplots()#,figsize=(40,30))
ax.scatter(VoxCoord[:,1], VoxCoord[:,3], c='gray', s=100, cmap='hot')  
#plt.axis('equal')
ax.set_aspect('equal')
ax.set_xlim([-70, 70])
ax.set_ylim([-45, 82])
ax.set_xticks([-65, 65])
ax.set_xticklabels(['left','right'])
ax.set_yticks([-35, 78])
ax.set_yticklabels(['bottom','top'], rotation=90)
#fig.savefig('coronal.png')    
plt.show()    


fig, ax = plt.subplots()#,figsize=(40,30))
ax.scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c='gray', s=100, cmap='hot')    
ax.set_aspect('equal')
ax.set_xlim([-70, 106])
ax.set_ylim([-45, 82])
ax.set_xticks([-65, 100])
ax.set_xticklabels(['front','back'])
ax.set_yticks([-35, 76])
ax.set_yticklabels(['bottom','top'], rotation=90)
#fig.savefig('sagitalL.png')    
plt.show()    

fig, ax = plt.subplots()#,figsize=(40,30))
ax.scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c='gray', s=100, cmap='hot')   
ax.set_aspect('equal')
ax.set_xlim([-106, 70])
ax.set_ylim([-45, 82])
ax.set_xticks([-100, 65])
ax.set_xticklabels(['back', 'front'])
ax.set_yticks([-35, 76])
ax.set_yticklabels(['bottom','top'], rotation=90)
#fig.savefig('sagitalR.png')    
plt.show()    


#a = np.where(np.diff(vlinesS)>1)[0] + 1
#MicroStart = np.concatenate((np.array(vlinesS[0]).reshape((1)), vlinesS[a]))
#MicroEnd = np.concatenate((vlinesS[a-1], np.array(vlinesS[-1]).reshape((1))))
#MicroMiddle = (MicroStart + MicroEnd)/2
#len(MicroMiddle) #9 array([ 285,  575,  617,  656,  748,  752, 1362, 1434, 1540])

LatIdx = []
for i in range(len(MicroMiddle)): 
    SigDip = 1 * sss[:, MicroMiddle[i]]
    SigVox = np.array(np.dot(from3To1, SigDip), dtype='f4')
    LatIdx.append([MicroMiddle[i]-100, np.sum(SigVox[RHemi])/np.sum((VoxCoord[:,1] >= 0)), 
                   np.sum(SigVox[LHemi])/np.sum(VoxCoord[:,1]<0)])


for i in range(len(MicroMiddle)):
    fig, ax = plt.subplots(2,2) #,figsize=(20,10))
    SigVoxGrpI = sss[:, MicroMiddle[i]]  
    nSigVoxGrpI= np.dot(from3To1, SigVoxGrpI)     
    ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,0].set_axis_off()
    ax[0,0].set_aspect('equal')    
    ax[0,0].set_xlim([-70, 70])
    ax[0,0].set_ylim([-106, 72])
    
    ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,1].set_aspect('equal')
    ax[0,1].set_axis_off()
    ax[0,1].set_xlim([-70, 70])
    ax[0,1].set_ylim([-45,82])
    
    ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c=nSigVoxGrpI[LHemi], s=100, cmap='hot')    
    ax[1,0].set_aspect('equal')
    ax[1,0].set_axis_off()    
    ax[1,0].set_xlim([-70, 106])
    ax[1,0].set_ylim([-45, 82])
    ax[1,0].set_title('Left %3.2f ' %( LatIdx[i][2] ))

    ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c=nSigVoxGrpI[RHemi], s=100, cmap='hot')   
    ax[1,1].set_aspect('equal')
    ax[1,1].set_axis_off()
    ax[1,1].set_xlim([-106, 70])
    ax[1,1].set_ylim([-45, 82])
    ax[1,1].set_title('Right %3.2f' %(  LatIdx[i][1]))
    fig.savefig('LatIdx105_%dms.png' %(MicroMiddle[i]-100))    
    
plt.show()    



for i in range(40):
    fig, ax = plt.subplots(2,2) #,figsize=(20,10))
    SigVoxGrpI = sss[:, i*40]  
    nSigVoxGrpI= np.dot(from3To1, SigVoxGrpI)     
    ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,0].set_axis_off()
    ax[0,0].set_aspect('equal')    
    ax[0,0].set_xlim([-70, 70])
    ax[0,0].set_ylim([-106, 72])
    ax[0,0].set_title('%d ms ' %(i*40-100))
    
    ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,1].set_aspect('equal')
    ax[0,1].set_axis_off()
    ax[0,1].set_xlim([-70, 70])
    ax[0,1].set_ylim([-45,82])
    
    ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c=nSigVoxGrpI[LHemi], s=100, cmap='hot')    
    ax[1,0].set_aspect('equal')
    ax[1,0].set_axis_off()    
    ax[1,0].set_xlim([-70, 106])
    ax[1,0].set_ylim([-45, 82])
#    ax[1,0].set_title('Left %3.2f ' %( LatIdx[i][2] ))

    ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c=nSigVoxGrpI[RHemi], s=100, cmap='hot')   
    ax[1,1].set_aspect('equal')
    ax[1,1].set_axis_off()
    ax[1,1].set_xlim([-106, 70])
    ax[1,1].set_ylim([-45, 82])
#    ax[1,1].set_title('Right %3.2f' %(  LatIdx[i][1]))
    fig.savefig('/Users/jesong1126/Work/VGT/tex/Figure11142014/LatIdxMovie/%04d.png' %(i*40))    
    

#ffmpeg -framerate 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p LatIdx105.mp4



for i in range(40):
    fig, ax = plt.subplots(2,2) #,figsize=(20,10))
    SigVoxGrpI = abs(sdenT1[:, i*40] ) 
    nSigVoxGrpI= np.dot(from3To1, SigVoxGrpI)     
    ax[0,0].scatter(VoxCoord[:,1], VoxCoord[:,2], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,0].set_axis_off()
    ax[0,0].set_aspect('equal')    
    ax[0,0].set_xlim([-70, 70])
    ax[0,0].set_ylim([-106, 72])
    ax[0,0].set_title('%d ms ' %(i*40-100))
    
    ax[0,1].scatter(VoxCoord[:,1], VoxCoord[:,3], c=nSigVoxGrpI, s=100, cmap='hot')  
    ax[0,1].set_aspect('equal')
    ax[0,1].set_axis_off()
    ax[0,1].set_xlim([-70, 70])
    ax[0,1].set_ylim([-45,82])
    
    ax[1,0].scatter(-VoxCoord[LHemi,2], VoxCoord[LHemi,3],  c=nSigVoxGrpI[LHemi], s=100, cmap='hot')    
    ax[1,0].set_aspect('equal')
    ax[1,0].set_axis_off()    
    ax[1,0].set_xlim([-70, 106])
    ax[1,0].set_ylim([-45, 82])
#    ax[1,0].set_title('Left %3.2f ' %( LatIdx[i][2] ))

    ax[1,1].scatter(VoxCoord[RHemi,2], VoxCoord[RHemi,3],  c=nSigVoxGrpI[RHemi], s=100, cmap='hot')   
    ax[1,1].set_aspect('equal')
    ax[1,1].set_axis_off()
    ax[1,1].set_xlim([-106, 70])
    ax[1,1].set_ylim([-45, 82])
#    ax[1,1].set_title('Right %3.2f' %(  LatIdx[i][1]))
    fig.savefig('/Users/jesong1126/Work/VGT/tex/Figure11142014/SourceDistMovie/%04d.png' %(i*40))    
    
#ffmpeg -framerate 1 -pattern_type glob -i '*.png' -c:v libx264 -pix_fmt yuv420p SourceDist105.mp4
##-----------------------------------------------------------------------------

 
#
#import matplotlib.animation as animation
#
#def update_line(num, data, line):
#    line.set_data(data[...,:num])
#    return line,
#
#fig1 = plt.figure() 
#
#data = np.random.rand(2, 25)
#l, = plt.plot([], [], 'r-')
#plt.xlim(0, 1)
#plt.ylim(0, 1)
#plt.xlabel('x')
#plt.title('test')
#line_ani = animation.FuncAnimation(fig1, update_line, 25, fargs=(data, l),
#    interval=50, blit=True)
##line_ani.save('lines.mp4')
#
#fig2 = plt.figure()
#
#x = np.arange(-9, 10)
#y = np.arange(-9, 10).reshape(-1, 1)
#base = np.hypot(x, y)
#ims = []
#for add in np.arange(15):
#    ims.append((plt.pcolor(x, y, base + add, norm=plt.Normalize(0, 30)),))
#
#im_ani = animation.ArtistAnimation(fig2, ims, interval=50, repeat_delay=3000,
#    blit=True)
#im_ani.save('im.mp4', metadata={'artist':'Guido'})
#
#plt.show() 
#
#  
 
##-----------------------------------------------------------------------------
  