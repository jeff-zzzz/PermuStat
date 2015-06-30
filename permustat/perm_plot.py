# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 12:33:27 2015
@author: jesong1126
"""

import matplotlib.pyplot as plt 
plt.ion()

def oneT(XMean, XSigma, XTobs, nSig, catName, msSamples):

    figOneT1, ax = plt.subplots(4,1) 
    ax[0].plot(msSamples, XMean.T)
    ax[0].set_xlim([msSamples[0], msSamples[-1]])
    ax[0].set_ylabel('mean')
    ax[0].set_title(catName)
    ax[1].plot(msSamples, XSigma.T)
    ax[1].set_xlim([msSamples[0], msSamples[-1]])
    ax[1].set_ylabel('sdv')
    ax[2].plot(msSamples, XTobs.T)
    ax[2].set_xlim([msSamples[0], msSamples[-1]])
    ax[2].set_ylabel('Normal')  
    ax[3].plot(msSamples, nSig)
    ax[3].set_xlim([msSamples[0], msSamples[-1]])
    ax[3].set_ylabel('Num Sig.Spaces')
    ax[3].set_xlabel('ms')
    
    figOneT2 = plt.figure()#, figsize=(20,10))
    ax = figOneT2.add_subplot(111)
    ax.imshow(XTobs, extent=[msSamples[0], msSamples[-1], XTobs.shape[0],0]) 
    ax.set_aspect('auto') 
    ax.set_title('T')
    
    
def twoT(XMean1, XSigma1, XT1, nSig1, catName1, XMean2, XSigma2, XT2, nSig2, catName2, XTobs, nSig, msSamples):    
    figTwoT1, ax = plt.subplots(4,2) 
    ax[0,0].plot(msSamples, XMean1.T)
    ax[0,0].set_xlim([msSamples[0], msSamples[-1]])
    ax[0,0].set_ylabel('mean')
    ax[0,0].set_title(catName1)
    ax[1,0].plot(msSamples, XSigma1.T)
    ax[1,0].set_xlim([msSamples[0], msSamples[-1]])
    ax[1,0].set_ylabel('sdv')
    ax[2,0].plot(msSamples, XT1.T)
    ax[2,0].set_xlim([msSamples[0], msSamples[-1]])
    ax[2,0].set_ylabel('Normal') # ax[2,0].axhline(vu, color='black', lw=2)
    ax[3,0].plot(msSamples, nSig1)
    ax[3,0].set_xlim([msSamples[0], msSamples[-1]])
    ax[3,0].set_ylabel('Num Sig.Spaces')
    ax[3,0].set_xlabel('ms')
    ax[0,1].plot(msSamples, XMean2.T)
    ax[0,1].set_xlim([msSamples[0], msSamples[-1]])
    ax[0,1].set_ylabel('mean')            
    ax[0,1].set_title(catName2)
    ax[1,1].plot(msSamples, XSigma2.T)
    ax[1,1].set_xlim([msSamples[0], msSamples[-1]])
    ax[1,1].set_ylabel('sdv')
    ax[2,1].plot(msSamples, XT2.T)
    ax[2,1].set_xlim([msSamples[0], msSamples[-1]])
    ax[2,1].set_ylabel('Normal')  
    ax[3,1].plot(msSamples, nSig2)
    ax[3,1].set_xlim([msSamples[0], msSamples[-1]])
    ax[3,1].set_ylabel('Num Sig.Spaces')
    ax[3,1].set_xlabel('ms')

    figTwoT2, ax = plt.subplots(2,1)  
    ax[0].plot(msSamples, XTobs.T)
    ax[0].set_ylabel('%s vs %s' %(catName1, catName2))            
    ax[0].set_xlim([msSamples[0], msSamples[-1]])
    ax[0].set_ylabel('Normal')  
    ax[1].plot(msSamples, nSig)
    ax[1].set_xlim([msSamples[0], msSamples[-1]])
    ax[1].set_ylabel('Num Sig.Spaces')
    ax[1].set_xlabel('ms') 

    figTwoT3 = plt.figure() 
    ax = figTwoT3.add_subplot(111)
    ax.imshow(XTobs, extent=[msSamples[0], msSamples[-1], XTobs.shape[0],0]) 
    ax.set_aspect('auto') 
    ax.set_title('T')


def pairT(XMean1, XSigma1, XT1, nSig1, catName1, XMean2, XSigma2, XT2, nSig2, catName2, XTobs, nSig, msSamples):    
    figPairT1, ax = plt.subplots(4,2) 
    ax[0,0].plot(msSamples, XMean1.T)
    ax[0,0].set_xlim([msSamples[0], msSamples[-1]])
    ax[0,0].set_ylabel('mean')
    ax[0,0].set_title(catName1)
    ax[1,0].plot(msSamples, XSigma1.T)
    ax[1,0].set_xlim([msSamples[0], msSamples[-1]])
    ax[1,0].set_ylabel('sdv')
    ax[2,0].plot(msSamples, XT1.T)
    ax[2,0].set_xlim([msSamples[0], msSamples[-1]])
    ax[2,0].set_ylabel('Normal') # ax[2,0].axhline(vu, color='black', lw=2)
    ax[3,0].plot(msSamples, nSig1)
    ax[3,0].set_xlim([msSamples[0], msSamples[-1]])
    ax[3,0].set_ylabel('Num Sig.Spaces')
    ax[3,0].set_xlabel('ms')
    ax[0,1].plot(msSamples, XMean2.T)
    ax[0,1].set_xlim([msSamples[0], msSamples[-1]])
    ax[0,1].set_ylabel('mean')            
    ax[0,1].set_title(catName2)
    ax[1,1].plot(msSamples, XSigma2.T)
    ax[1,1].set_xlim([msSamples[0], msSamples[-1]])
    ax[1,1].set_ylabel('sdv')
    ax[2,1].plot(msSamples, XT2.T)
    ax[2,1].set_xlim([msSamples[0], msSamples[-1]])
    ax[2,1].set_ylabel('Normal')  
    ax[3,1].plot(msSamples, nSig2)
    ax[3,1].set_xlim([msSamples[0], msSamples[-1]])
    ax[3,1].set_ylabel('Num Sig.Spaces')
    ax[3,1].set_xlabel('ms')

    figPairT2, ax = plt.subplots(2,1)  
    ax[0].plot(msSamples, XTobs.T)
    ax[0].set_ylabel('%s vs %s' %(catName1, catName2))            
    ax[0].set_xlim([msSamples[0], msSamples[-1]])
    ax[0].set_ylabel('Normal')  
    ax[1].plot(msSamples, nSig)
    ax[1].set_xlim([msSamples[0], msSamples[-1]])
    ax[1].set_ylabel('Num Sig.Spaces')
    ax[1].set_xlabel('ms') 

    figPairT3 = plt.figure()
    ax = figPairT3.add_subplot(111)
    ax.imshow(XTobs, extent=[msSamples[0], msSamples[-1], XTobs.shape[0],0]) 
    ax.set_aspect('auto') 
    ax.set_title('T')


def anova(XF, PvalF, SigF, msSamples):    
    fig, ax = plt.subplots(1,2, figsize=(12,6))
    ax[0].imshow(XF, cmap='hot', extent=[msSamples[0], msSamples[-1], XF.shape[0],0])
    ax[0].set_aspect('auto') 
    ax[0].set_title('F')
    ax[1].imshow(PvalF, cmap='hot_r', extent=[msSamples[0], msSamples[-1], XF.shape[0],0])
    ax[1].set_aspect('auto') 
    ax[1].set_title('Pval')

def permT(XTobs, PvalST,msSamples):    
#    import numpy as np 
#    logPvalST = -np.log(PvalST)
    fig, ax = plt.subplots(1,2, figsize=(12,6))
    ax[0].imshow(XTobs, extent=[msSamples[0], msSamples[-1], XTobs.shape[0],0])
    ax[0].set_aspect('auto') 
    ax[0].set_title('observed T')
    ax[1].imshow(PvalST, cmap='hot_r', extent=[msSamples[0], msSamples[-1], XTobs.shape[0],0])
#    ax[1].imshow(logPvalST, cmap= 'hot')
    ax[1].set_aspect('auto') 
    ax[1].set_title('Pval')



         
    
    
    