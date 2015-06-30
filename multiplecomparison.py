# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 11:20:57 2015
@author: jesong1126
"""
 
import numpy as np   
import matplotlib.pyplot as plt
#from swcm import Imatrix, Lcurve #from matplotlib import gridspec
from mff import read_mff_header, read_mff_data #,getEpochInfos, mff_getSummaryInfo
from braink import read_lfm                 
import os 
import glob 
from permustat import permustat_fun, perm_plot 


filePath = '/Users/jesong1126/Python27/PermuStat/data/SEP_107_0046_seg_486.mff'

hdr = read_mff_header.read_mff_header(filePath)
nC = hdr['nChans']
nSamples = hdr['nSamples']
srate = hdr['Fs']
nSamplesPre = hdr['nSamplesPre']
nTrials = hdr['nTrials']
summaryInfo = hdr['orig'] 
trialsName = summaryInfo['epochLabels']   
categoryName = list(set(trialsName))
nCategory = len(categoryName)
baseline = (nSamplesPre * 1000 / srate)  
msSamples = np.arange(0, nSamples,1) * 1000/srate  - baseline 
xlimMin = msSamples[0]; xlimMax = msSamples[-1];

s = read_mff_data.read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr)    

whichCat = np.zeros((nTrials, nCategory))
for j in range(nCategory):
    for i in range(nTrials):
        if (trialsName[i] == categoryName[j]):
            whichCat[i,j] = 1
whichCat = np.array(whichCat, dtype='i')
 
 
catName1 = categoryName[0]; catName2 = categoryName[1] ;
print('Two sample T-test: %s vs %s' %(catName1, catName2)) 
X1 = s[:,:, (whichCat[:,0]>0)]; X2 = s[:,:, (whichCat[:,1]>0)];        

nX1 = X1.shape[2]
SigLevel = 0.05     

XMean1 = np.mean(X1, axis=2); XSigma1 = np.std(X1, axis=2);
statOut1 = permustat_fun.one_sample_t(X1, SigLevel/2) 

Gfp1 = np.std(XMean1,axis=0) 
fig = plt.figure()
plt.plot(msSamples, Gfp1)
plt.axvline(43, color='black', lw=2)

fig = plt.figure()
plt.imshow(statOut1['Sig']) 

tsid = np.where(msSamples == 43)[0]
statOut1['Sig'][:, tsid]

up = XMean1[:, tsid] + 3* XSigma1[:, tsid] / np.sqrt(nX1)
bt = XMean1[:, tsid] - 3* XSigma1[:, tsid] / np.sqrt(nX1)

fig = plt.figure()
ax1 = fig.add_subplot(211)
ax1.plot(np.arange(nC), XMean1[:, tsid])#, usevlines=True, maxlags=50, normed=True, lw=2)
ax1.vlines(np.arange(nC), bt ,up)
ax1.grid(True)
ax1.axhline(0, color='black', lw=2)
ax1.set_xlim([0, nC])
ax2 = fig.add_subplot(212, sharex=ax1)
ax2.set_xlim([0, nC])
ax2.vlines(np.arange(nC), [0], statOut1['Sig'][:, tsid])
plt.show() 


from scipy.stats import t 

chI = 10; chJ = 70; Xi = X1[chI, tsid, :]; Xj = X1[chJ, tsid, :]; 

MitI = np.mean(Xi); MitJ = np.mean(Xj); nume = MitI - MitJ
SigmaI = np.var(Xi); SigmaJ = np.var(Xj); denume = np.sqrt((SigmaI /Xi.shape[1]) + (SigmaJ /Xj.shape[1]))
T = np.divide(nume, denume) 
temp = 1-t.cdf(T, 2*nX1-2)
Pval = 1-temp 
(Pval < SigLevel/2)





#    return statOut #[XT, PvalST, SigST, nSig, nSigT]
     

#chid = 10; statOut1['Sig'][chid,:]; 
#
#up = XMean1[chid,:] + 3* XSigma1[chid,:]/np.sqrt(nX1)
#bt = XMean1[chid,:] - 3* XSigma1[chid,:]/np.sqrt(nX1)
#
#fig = plt.figure()
#ax1 = fig.add_subplot(211)
#ax1.plot(msSamples, XMean1[chid,:])#, usevlines=True, maxlags=50, normed=True, lw=2)
#ax1.vlines(msSamples, bt ,up)
#ax1.grid(True)
#ax1.axhline(0, color='black', lw=2) 
#ax1.axhline(np.mean(XMean1[chid,:] ), color='red', lw=2) 
##ax1.set_xlim([msSamples[0], msSamples[1]])
#ax2 = fig.add_subplot(212, sharex=ax1)
#ax2.vlines(msSamples, [0], statOut1['Sig'][chid,:])
##ax2.set_xlim([msSamples[0], msSamples[1]])
#plt.show()




#
#XMean2 = np.mean(X2, axis=2); XSigma2 = np.std(X2, axis=2);   
#
#statOut1 = permustat_fun.one_sample_t(X1, SigLevel/2) 
#statOut2 = permustat_fun.one_sample_t(X2, SigLevel/2) 
#
#
#
#if self.rbPara.isChecked():
#    statOut = permustat_fun.two_sample_t(X1, X2, SigLevel) 
#    perm_plot.twoT(XMean1, XSigma1, statOut1['T'], statOut1['nSigS'], catName1, 
#                   XMean2, XSigma2, statOut2['T'], statOut2['nSigS'], catName2, statOut['T'], statOut['nSigS'], msSamples)
#else:
#    statOut = permustat_fun.two_sample_permu(X1, X2, SigLevel/2, NumSim)
#    logPvalST = np.log(statOut['Pval'])
#    perm_plot.permT(statOut['T'], logPvalST, msSamples) 
#
#config.Stat = statOut['T']                  

 
# 
# 
#nE = nC - 1
#HMdir = '/Users/jesong1126/Python27/GeoPy/data/Head_Model_Data07' 
#lfm_fname = glob.glob(HMdir+'/*.lfm')
#if os.path.exists(HMdir+'/fdmForwardMatrixOriented'):          
#    K  = read_lfm.forward(HMdir+'/fdmForwardMatrixOriented', nE)  
#    print("fdmForwardMatrixOriented is loaded. ")            
#elif len(lfm_fname) > 0:  # elif os.path.exists(HMdir+'/Leadfield.lfm'): 
#    K  = read_lfm.lfm(lfm_fname[0], nE)                          
#    print("Leadfield.lfm is loaded. ")      
#else: 
#    print("Either fdmForwardMatrixOriented or *.lfm are not loaded. ")            
#        
#
#print K.shape, s.shape
#nSamples = s.shape[1] 
#if len(s.shape)>2 :
#    nTrials = s.shape[2] 
#else: 
#    nTrials = 1            
#    
#nV = K.shape[1]


##----------------------------------------------------------------------------- 
#phi = s[:,0]
#np.seterr(divide = 'ignore')
# 
#Logalpha = np.arange(-6, 7, 1)
#nLcurve = len(Logalpha) 
#Lcurve = [] # Lcurve Residual L2 
#for i in Logalpha: 
#    alphai = np.power(10, i)
#    Imat = Imatrix.MN(alphai, K)
#    beta = np.dot(Imat, phi) 
#    phihat = np.dot(K, beta) 
#    phidiff = phi-phihat 
##    Resid2 = np.sqrt(phi.shape[0]) * np.std(phi_diff)     
##    L2 =  np.sqrt(beta.shape[0]) * np.std(beta)
##    Lcurve = np.concatenate((Lcurve, [i, Resid2, L2]), axis=0)
#    Resid2 = np.sqrt(sum(phidiff*phidiff)) 
#    L2 = np.sqrt(sum(beta*beta))
#    Lcurve = np.concatenate((Lcurve, [i, Resid2, L2]), axis=0)
#
#Lcurve =  Lcurve.reshape(nLcurve,3)    
#Lcurve0 = np.array([Lcurve[:,0], Lcurve[:,1]/np.max(Lcurve[:,1]), Lcurve[:,2]/np.max(Lcurve[:,2])])    
#Lcurve = Lcurve0.T
#  
#Resid2D1 = np.diff(Lcurve[:,1]); Resid2D2 = np.diff(Resid2D1);
#L2D1 = np.diff(Lcurve[:,2]); L2D2 = np.diff(L2D1);
#
#pho = np.array([Lcurve[:,1], np.concatenate((Resid2D1,[0])), np.concatenate((Resid2D2,[0,0]))])
#pho = pho.T
#eta = np.array([Lcurve[:,2], np.concatenate((L2D1,[0])), np.concatenate((L2D2,[0,0]))])
#eta = eta.T
# 
#nume =(np.multiply(pho[:,1], eta[:,2]) - np.multiply(pho[:,2], eta[:,1])) 
#denume = np.power((np.multiply(pho[:,1], pho[:,1]) -  np.multiply(eta[:,1],eta[:,1])), 3/2) 
#
##(denume == 0) 
#kappa = np.divide(nume,denume)
##kappa = np.array(np.divide(nume,denume))
#
#figLC, axLC = plt.subplots(2,2)
#axLC[0,0].plot(Lcurve[:,0], Lcurve[:,1]); axLC[0,0].set_xlim((-6, 6)); 
#axLC[0,0].set_xlabel('alpha'); axLC[0,0].set_ylabel('Residual')
#axLC[0,1].plot(Lcurve[:,0], Lcurve[:,2]); axLC[0,1].set_xlim((-6, 6)); 
#axLC[0,1].set_xlabel('alpha'); axLC[0,1].set_ylabel('|J|')
#axLC[1,0].scatter(Lcurve[:,1], Lcurve[:,2]);
#axLC[1,0].set_xlabel('Residual'); axLC[1,0].set_ylabel('|J|')
#axLC[1,1].scatter(Lcurve[:,0], kappa);
#axLC[1,1].set_xlabel('alpha'); axLC[1,1].set_ylabel('kappa')

# 
#fig, ax = plt.subplots()
#ax.scatter(Lcurve[:,1],Lcurve[:,2])
#plt.xlabel('Residual'); plt.ylabel('|J|')
#for i, txt in enumerate(Logalpha):
#    ax.annotate(txt, (Lcurve[i,1],Lcurve[i,2]), xytext = (5,1), textcoords = 'offset points')
#
#maxKappa = np.where(kappa == np.nanmax(kappa))[0]
#logalphaS = Lcurve[maxKappa[0],0]
# 
#logalpha = logalphaS, Lcurve
#  
###-----------------------------------------------------------------------------  
#
#poweralphaS = -2
#poweralphaS = int(poweralphaS)
#alpha = pow(10, poweralphaS)  
#alphaL = 1.0e-4
#print("alpha = %3.2e and alphaL = %3.2e"  %(alpha, alphaL))
#
#Imat = Imatrix.MN(alpha, K)
##Imat = Imatrix.sMN(alpha, K)
#
#if nTrials > 1:
#    sden = np.zeros((nV, nSamples, nTrials)) 
#    for i in range(nTrials):             
#        sden[:,:,i] = np.dot(Imat, s[:,:,i]) 
#else:
#    sden = np.dot(Imat, s) 
#    
#print("Source is estimated")                    
#
#MS = 0         
#msI = np.where(msSamples == MS)[0]
#
#abssden = np.abs(sden[:,msI]) 
#abssden.sort(axis=0) 
#cutoff = abssden[int(np.rint(( 95/100.0 ) * nV))]
#abssden = np.abs(sden[:,msI]) 
#abssdenTh = (abssden > cutoff) * abssden 
#sdenTh = (abssden > cutoff) * sden[:,msI]
# 
#np.sum(abssden > cutoff)
#
#fig, ax = plt.subplots(4,1) 
#ax[0].plot(np.arange(nV), sden[:,msI])
#ax[1].plot(np.arange(nV), abssden)
#ax[2].plot(np.arange(nV), abssdenTh)
#ax[3].plot(np.arange(nV), sdenTh)
#
#EEG = s
#sdenI = sden 
#AllForward = np.dot(K, sdenI[:,msI])        
#AllResidual = EEG[:,msI] - AllForward   
#PartForward = np.dot(K, sdenTh)
#PartResidual = EEG[:, msI] - PartForward 
#
#sEEGI = EEG[:,msI] / np.std(EEG[:,msI]) 
#sAllForward = AllForward / np.std(AllForward) 
#sAllResidual = sEEGI- sAllForward 
#sPartForward = PartForward / np.std(PartForward) 
#sPartResidual = sEEGI- sPartForward  
#
#fig, ax = plt.subplots(5,2) 
#ax[0,0].plot(np.arange(K.shape[0]), EEG[:,msI])
#ax[1,0].plot(np.arange(K.shape[0]), AllForward)
#ax[2,0].plot(np.arange(K.shape[0]), AllResidual)
#ax[3,0].plot(np.arange(K.shape[0]), PartForward)
#ax[4,0].plot(np.arange(K.shape[0]), PartResidual)
#
#ax[0,1].plot(np.arange(K.shape[0]), sEEGI)
#ax[1,1].plot(np.arange(K.shape[0]), sAllForward)
#ax[2,1].plot(np.arange(K.shape[0]), sAllResidual)
#ax[3,1].plot(np.arange(K.shape[0]), sPartForward)
#ax[4,1].plot(np.arange(K.shape[0]), sPartResidual)
#
#fig, ax = plt.subplots(5,2) 
#ax[0,1].scatter(np.arange(K.shape[1]), sdenI[:,msI])
#ax[1,1].scatter(np.arange(K.shape[1]), abssden)
#ax[2,1].scatter(np.arange(K.shape[1]), sdenTh)
#
#        figSden, ax = plt.subplots(4,1) 
#        ax[0].scatter(np.arange(K.shape[1]), sdenI[:,msI])
#        ax[0].set_xlim([0, K.shape[1]])
#        ax[0,0].set_title('source')
#        ax[1,0].scatter(np.arange(K.shape[1]), abssden)
#        ax[1,0].set_xlim([0, K.shape[1]])
#        ax[1,0].set_title('abs(source) ')
#        ax[2,0].scatter(np.arange(K.shape[1]), sdenTh)
#        ax[2,0].set_xlim([0, K.shape[1]])
#        ax[2,0].set_title('Threshold(source) ')
#        ax[3,0].scatter(np.arange(K.shape[1]), abssdenTh)
#        ax[3,0].set_xlim([0, K.shape[1]])
#        ax[3,0].set_title('abs(Threshold(source))')
#
#
### 
#    Logalpha = np.arange(-6, 7, 1)
#    nLcurve = len(Logalpha) 
#    Lcurve = [] # Lcurve Residual L2 
#    for i in Logalpha: 
#        alphai = 10 ** i
#        Imat = Imatrix.MN(alphai, K)
#        beta = np.dot(Imat, phi) 
#        phihat = np.dot(K, beta) 
#        Resid2 = np.sqrt(sum((phi-phihat)*(phi-phihat))) 
#        L2 = np.sqrt(sum(beta*beta))
#        Lcurve = np.concatenate((Lcurve, [i, Resid2, L2]), axis=0)
#
#    Lcurve =  Lcurve.reshape(nLcurve,3)    
#    Lcurve0 = np.array([Lcurve[:,0], Lcurve[:,1]/max(Lcurve[:,1]), Lcurve[:,2]/max(Lcurve[:,2])])    
#    Lcurve = Lcurve0.T
#  
#    Resid2D1 = np.diff(Lcurve[:,1]); Resid2D2 = np.diff(Resid2D1);
#    L2D1 = np.diff(Lcurve[:,2]); L2D2 = np.diff(L2D1);
#
#    pho = np.array([Lcurve[:,1], np.concatenate((Resid2D1,[0])), np.concatenate((Resid2D2,[0,0]))])
#    pho = pho.T
#    eta = np.array([Lcurve[:,2], np.concatenate((L2D1,[0])), np.concatenate((L2D2,[0,0]))])
#    eta = eta.T
# 
#    nume =(np.multiply(pho[:,1], eta[:,2]) - np.multiply(pho[:,2], eta[:,1])) 
#    denume = np.power((np.multiply(pho[:,1], pho[:,1]) -  np.multiply(eta[:,1],eta[:,1])), 3/2)
#    kappa = np.array(np.divide(nume,denume))
#
#    plt.subplot(221)
#    plt.plot(Lcurve[:,0], Lcurve[:,1]); plt.xlim((-6, 6)); 
#    plt.xlabel('alpha'); plt.ylabel('Residual')
#    plt.subplot(222)
#    plt.plot(Lcurve[:,0], Lcurve[:,2]); plt.xlim((-6, 6)); 
#    plt.xlabel('alpha'); plt.ylabel('|J|')
#    plt.subplot(223)
#    plt.scatter(Lcurve[:,1], Lcurve[:,2]);
#    plt.xlabel('Residual'); plt.ylabel('|J|')
#    plt.subplot(224)
#    plt.scatter(Lcurve[:,0], kappa);
#    plt.xlabel('alpha'); plt.ylabel('kappa')
#
#    fig, ax = plt.subplots()
#    ax.scatter(Lcurve[:,1],Lcurve[:,2])
#    plt.xlabel('Residual'); plt.ylabel('|J|')
#    for i, txt in enumerate(Logalpha):
#        ax.annotate(txt, (Lcurve[i,1],Lcurve[i,2]), xytext = (5,1), textcoords = 'offset points')
#
#    maxKappa = (kappa == max(kappa)).nonzero()
#    logalphaS = Lcurve[maxKappa[0],0]
#     
#    logalpha = logalphaS, Lcurve
#      
##plt.hist(sden[:, msI])
##np.mean(sden[:, msI])
#np.std( sden[:, msI])

#from pandas import DataFrame 
#import pandas as pd 
#pd.options.display.float_format = '{:,.3f}'.format
#
#gsn257 = pd.read_table('/Users/jesong1126/Python27/GeoPy/nscolor_hgsn/GSN257ToPy2.dat')  
#frame = DataFrame(gsn257 , columns=['ChLabel','X3','Y3','Z3','X2','Y2'])
#x2 = frame.values[:,4]; y2 = frame.values[:,5]  

#    X1 = sden[:,:, (whichCat[:, config.TwoTCatVal1]>0)]        
#    X2 = sden[:,:, (whichCat[:, config.TwoTCatVal2]>0)]        


 
 
