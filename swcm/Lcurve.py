import numpy as np
from swcm import Imatrix
from matplotlib import pyplot as plt
plt.ion()
import glob
from braink import read_bk

def MN(phi, K):
    
    Logalpha = np.arange(-6, 7, 1)
    nLcurve = len(Logalpha) 
    Lcurve = [] # Lcurve Residual L2 
    for i in Logalpha: 
        alphai = 10 ** i
        Imat = Imatrix.MN(alphai, K)
        beta = np.dot(Imat, phi) 
        phihat = np.dot(K, beta) 
        phidiff = phi-phihat 
        Resid2 = np.sqrt(sum(phidiff*phidiff)) 
        L2 = np.sqrt(sum(beta*beta))
        Lcurve = np.concatenate((Lcurve, [i, Resid2, L2]), axis=0)

    Lcurve =  Lcurve.reshape(nLcurve,3)    
    Lcurve0 = np.array([Lcurve[:,0], Lcurve[:,1]/max(Lcurve[:,1]), Lcurve[:,2]/max(Lcurve[:,2])])    
    Lcurve = Lcurve0.T
  
    Resid2D1 = np.diff(Lcurve[:,1]); Resid2D2 = np.diff(Resid2D1);
    L2D1 = np.diff(Lcurve[:,2]); L2D2 = np.diff(L2D1);

    pho = np.array([Lcurve[:,1], np.concatenate((Resid2D1,[0])), np.concatenate((Resid2D2,[0,0]))])
    pho = pho.T
    eta = np.array([Lcurve[:,2], np.concatenate((L2D1,[0])), np.concatenate((L2D2,[0,0]))])
    eta = eta.T
 
    nume =(np.multiply(pho[:,1], eta[:,2]) - np.multiply(pho[:,2], eta[:,1])) 
    denume = np.power((np.multiply(pho[:,1], pho[:,1]) - np.multiply(eta[:,1], eta[:,1])), 3/2) 
    kappa = np.divide(nume,denume) 
    
    figLC, axLC = plt.subplots(2,2)
    axLC[0,0].plot(Lcurve[:,0], Lcurve[:,1]); 
    axLC[0,0].set_xlim((Logalpha[0],Logalpha[-1])); 
    axLC[0,0].set_xlabel('alpha'); 
    axLC[0,0].set_ylabel('Residual')
    axLC[0,1].plot(Lcurve[:,0], Lcurve[:,2]); 
    axLC[0,1].set_xlim((Logalpha[0],Logalpha[-1])); 
    axLC[0,1].set_xlabel('alpha'); 
    axLC[0,1].set_ylabel('|J|')
    axLC[1,0].scatter(Lcurve[:,1], Lcurve[:,2]);
    axLC[1,0].set_xlabel('Residual'); 
    axLC[1,0].set_ylabel('|J|')
    axLC[1,1].scatter(Lcurve[:,0], kappa);
    axLC[1,1].set_xlabel('alpha'); 
    axLC[1,1].set_ylabel('kappa')

    fig, ax = plt.subplots()
    ax.scatter(Lcurve[:,1],Lcurve[:,2])
    plt.xlabel('Residual'); plt.ylabel('|J|')
    for i, txt in enumerate(Logalpha):
        ax.annotate(txt, (Lcurve[i,1],Lcurve[i,2]), xytext = (5,1), textcoords = 'offset points')

    maxKappa = np.where(abs(kappa) == np.nanmax(abs(kappa)))[0]
    logalphaS = Lcurve[maxKappa[0],0]
     
    logalpha = logalphaS, Lcurve
     
    return logalpha  
##------------------------------------------------------------------------------

def CSL(phi, K, bkdir):
    ## alphaS 
    Logalpha = np.arange(-6, 7, 1)
    nLcurve = len(Logalpha) 
    LcurveS = [] # Lcurve Residual L2 
    for i in Logalpha: 
        alphai = 10 ** i
        Imat = Imatrix.MN(alphai, K)
        beta = np.dot(Imat, phi) 
        phihat = np.dot(K, beta) 
        phidiff = phi-phihat
        Resid2 = np.sqrt(sum(phidiff*phidiff)) 
        L2 = np.sqrt(sum(beta*beta))
        LcurveS = np.concatenate((LcurveS, [i, Resid2, L2]), axis=0)

    LcurveS =  LcurveS.reshape(nLcurve,3)    
    LcurveS0 = np.array([LcurveS[:,0], LcurveS[:,1]/max(LcurveS[:,1]), LcurveS[:,2]/max(LcurveS[:,2])])    
    LcurveS = LcurveS0.T
  
    Resid2D1 = np.diff(LcurveS[:,1]); Resid2D2 = np.diff(Resid2D1);
    L2D1 = np.diff(LcurveS[:,2]); L2D2 = np.diff(L2D1);

    pho = np.array([LcurveS[:,1], np.concatenate((Resid2D1,[0])), np.concatenate((Resid2D2,[0,0]))])
    pho = pho.T
    eta = np.array([LcurveS[:,2], np.concatenate((L2D1,[0])), np.concatenate((L2D2,[0,0]))])
    eta = eta.T
 
    nume =(np.multiply(pho[:,1], eta[:,2]) - np.multiply(pho[:,2], eta[:,1])) 
    denume = np.power((np.multiply(pho[:,1], pho[:,1]) -  np.multiply(eta[:,1],eta[:,1])), 3/2)
    kappa = np.array(np.divide(nume,denume))
    maxKappa = np.where(abs(kappa) == np.nanmax(abs(kappa)))[0]
    logalphaS = LcurveS[maxKappa[0],0] 

    figLC, axLC = plt.subplots(2,2)
    axLC[0,0].plot(LcurveS[:,0], LcurveS[:,1]); 
    axLC[0,0].set_xlim((Logalpha[0],Logalpha[-1])); 
    axLC[0,0].set_xlabel('alpha'); 
    axLC[0,0].set_ylabel('Residual')
    axLC[0,1].plot(LcurveS[:,0], LcurveS[:,2]); 
    axLC[0,1].set_xlim((Logalpha[0],Logalpha[-1])); 
    axLC[0,1].set_xlabel('alpha'); 
    axLC[0,1].set_ylabel('|J|')
    axLC[1,0].scatter(LcurveS[:,1], LcurveS[:,2]);
    axLC[1,0].set_xlabel('Residual'); 
    axLC[1,0].set_ylabel('|J|')
    axLC[1,1].scatter(LcurveS[:,0], kappa);
    axLC[1,1].set_xlabel('alpha'); 
    axLC[1,1].set_ylabel('kappa')


    fig, ax = plt.subplots()
    ax.scatter(LcurveS[:,1],LcurveS[:,2])
    plt.xlabel('Residual'); plt.ylabel('|J|')
    for i, txt in enumerate(Logalpha):
        ax.annotate(txt, (LcurveS[i,1],LcurveS[i,2]), xytext = (5,1), textcoords = 'offset points')

    ## alphaL 
    Ne, Nd = K.shape   
    for dsafilename1 in glob.iglob(bkdir+"/*.Adjacents_1_*"):         
        A_local1 = read_bk.dsa(dsafilename1)
    for dsafilename2 in glob.iglob(bkdir+"/*.Adjacents_2_*"): 
        A_local2 = read_bk.dsa(dsafilename2)
    
    nD1 = A_local1.shape[0]-1
    nD2 = A_local2.shape[0]-1 
    Nei1 = A_local1[:-1, :-1]
    Nei2 = A_local2[:-1, :-1]  

    r0 = np.concatenate((Nei1, np.zeros((nD1,nD2))), axis=1)  
    r1 = np.concatenate((np.zeros((nD2,nD1)), Nei2), axis=1) 
    Neighbor = np.concatenate((r0, r1), axis=0)

    B = np.identity(Nd)
    for i in range(Nd):
        Ni = np.sum(Neighbor[i,:])-1   
        if Ni > 0:
            B[i,:] = -Neighbor[i,:]/Ni 
            B[i,i] = 1 
        else:
            B[i,:] = -1/(Nd-1)
            B[i,i] = 1    

    alphaS = 10 ** logalphaS 
    Logalpha = np.arange(-6, 7, 1)
    nLcurve = len(Logalpha) 
    LcurveL = [] # np.zeros((nLcurve, 4)) 
    for i in Logalpha: 
        alphaLi = 10 ** i
        Imat = Imatrix.CSL(alphaS, alphaLi, K, bkdir)   
        beta = np.dot(Imat, phi) 
        phihat = np.dot(K, beta) 
        phidiff = phi-phihat
        Resid2 = np.sqrt(sum(phidiff*phidiff)) 
        L2 = np.sqrt(sum(beta * beta))
        Bbeta = np.dot(B, beta)
        BL2 = np.sqrt(sum(Bbeta * Bbeta))
        LcurveL = np.concatenate((LcurveL, [i, Resid2, L2, BL2]), axis=0) 
    
    LcurveL =  LcurveL.reshape(nLcurve, 4)     
    LcurveL0 = np.array([LcurveL[:,0], LcurveL[:,1]/max(LcurveL[:,1]), LcurveL[:,2]/max(LcurveL[:,2]), LcurveL[:,3]/max(LcurveL[:,3])])    
    LcurveL = LcurveL0.T
  
    Resid2D1 = np.diff(LcurveL[:,1]); Resid2D2 = np.diff(Resid2D1);
    L2D1 = np.diff(LcurveL[:,3]); L2D2 = np.diff(L2D1);

    pho = np.array([LcurveL[:,1], np.concatenate((Resid2D1,[0])), np.concatenate((Resid2D2,[0,0]))])
    pho = pho.T
    eta = np.array([LcurveL[:,3], np.concatenate((L2D1,[0])), np.concatenate((L2D2,[0,0]))])
    eta = eta.T
 
    nume =(np.multiply(pho[:,1], eta[:,2]) - np.multiply(pho[:,2], eta[:,1])) 
    denume = np.power((np.multiply(pho[:,1], pho[:,1]) -  np.multiply(eta[:,1],eta[:,1])), 3/2)
    kappa = np.array(np.divide(nume,denume))
    maxKappa = np.where(abs(kappa) == np.nanmax(abs(kappa)))[0]
    logalphaL = LcurveL[maxKappa[0],0] 

    figLCL, axLCL = plt.subplots(2,2)
    axLCL[0,0].plot(LcurveL[:,0], LcurveL[:,1]); 
    axLCL[0,0].set_xlim((Logalpha[0],Logalpha[-1])); 
    axLCL[0,0].set_xlabel('alpha'); 
    axLCL[0,0].set_ylabel('Residual')
    axLCL[0,1].plot(LcurveL[:,0], LcurveL[:,3]); 
    axLCL[0,1].set_xlim((Logalpha[0],Logalpha[-1])); 
    axLCL[0,1].set_xlabel('alpha'); 
    axLCL[0,1].set_ylabel('|BJ|')
    axLCL[1,0].scatter(LcurveL[:,1], LcurveL[:,3]);
    axLCL[1,0].set_xlabel('Residual'); 
    axLCL[1,0].set_ylabel('|BJ|')
    axLCL[1,1].scatter(LcurveL[:,0], kappa);
    axLCL[1,1].set_xlabel('alpha'); 
    axLCL[1,1].set_ylabel('kappa')

    fig, ax = plt.subplots()
    ax.scatter(LcurveL[:,1],LcurveL[:,3])
    plt.xlabel('Residual'); plt.ylabel('|BJ|')
    for i, txt in enumerate(Logalpha):
        ax.annotate(txt, (LcurveL[i,1],LcurveL[i,3]), xytext = (5,1), textcoords = 'offset points')

    logalpha = logalphaS, LcurveS, logalphaL, LcurveL
    
    return logalpha  