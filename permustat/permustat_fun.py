# -*- coding: utf-8 -*-
"""
Created on Mon Mar 23 14:03:57 2015
@author: jesong1126
"""

def one_sample_t(X, SigLevel):
    from scipy.stats import t 
    import numpy as np

    nume = np.mean(X, axis=2)    
    Xvar = np.var(X, axis=2)
    denume = np.sqrt(Xvar / X.shape[2])
    XT = np.divide(nume, denume) 
    PvalST = np.zeros(XT.shape)
    for ii in range(XT.shape[0]): 
        for jj in range(XT.shape[1]):
            temp = 1-t.cdf(XT[ii,jj], X.shape[2] - 1)
            if temp > 0.5 : 
                PvalST[ii, jj] = 1 - temp
            else: 
                PvalST[ii, jj] = temp
                
    SigST = 1.0 * (PvalST < SigLevel/2)
    nSigS  = np.sum(SigST, axis=0)     
    nSigT  = np.sum(SigST, axis=1) 
    statOut = {'T':XT,'Pval':PvalST,'Sig':SigST,'nSigS':nSigS,'nSigT':nSigT} 
    return  statOut #[XT, PvalST, SigST, nSig, nSigT]
    

def two_sample_t(X0, X1, SigLevel):
    from scipy.stats import t 
    import numpy as np
    from permustat import perm_t
    
    XT = perm_t.get_T(X0, X1)
    PvalST = np.zeros(XT.shape)
    for ii in range(XT.shape[0]): 
        for jj in range(XT.shape[1]):
            temp = 1-t.cdf(XT[ii,jj], X0.shape[2] + X1.shape[2]-2)
            if temp > 0.5 : 
                PvalST[ii, jj] = 1 - temp
            else: 
                PvalST[ii, jj] = temp
                
    SigST = 1.0 * (PvalST < SigLevel/2)
    nSigS  = np.sum(SigST, axis=0) 
    nSigT  = np.sum(SigST, axis=1) 
    statOut = {'T':XT,'Pval':PvalST,'Sig':SigST,'nSigS':nSigS,'nSigT':nSigT} 

    return statOut #[XT, PvalST, SigST, nSig, nSigT]
     

def pairt(X0, X1, SigLevel):     
    from scipy.stats import t    
    import numpy as np

    Diff = X0 - X1
    nume = np.mean(Diff, axis=2) 
    denume = np.std(Diff, axis=2) / np.sqrt(Diff.shape[2])
    pairT = np.divide(nume, denume) 
    XT = pairT 
    PvalST = np.zeros(XT.shape)
    for ii in range(XT.shape[0]): 
        for jj in range(XT.shape[1]):
            temp = 1-t.cdf(XT[ii,jj], Diff.shape[2]-1)
            if temp > 0.5 : 
                PvalST[ii, jj] = 1 - temp
            else: 
                PvalST[ii, jj] = temp
                
    SigST = 1.0 * (PvalST < SigLevel/2)
    nSigS  = np.sum(SigST, axis=0) 
    nSigT  = np.sum(SigST, axis=1) 
    statOut = {'T':XT,'Pval':PvalST,'Sig':SigST,'nSigS':nSigS,'nSigT':nSigT} 

    return statOut #[XT, PvalST, SigST, nSig, nSigT]
 
    
def anova_f(s, whichCat, nCategory, SigLevel):
    import numpy as np
    from scipy.stats import f 

    Xmean = np.mean(s, axis=2)
    
    df_b = nCategory - 1 
    df_w = s.shape[2] - nCategory

    SS_b = np.zeros((s.shape[0], s.shape[1])) 
    SS_w = np.zeros((s.shape[0], s.shape[1]))
    for kk in range(nCategory):
        Xk = s[:,:, whichCat[:, kk]>0]  
        Xkmean = np.mean(Xk, axis=2)
        dtemp = Xkmean-Xmean
        SS_b = SS_b + Xk.shape[2] * (dtemp * dtemp)
        for i in range(Xk.shape[2]):
            temp = (Xk[:,:,i] - Xkmean)    
            SS_w = SS_w + (temp * temp) 
    SS_t = SS_b + SS_w; 
    
    MS_b = SS_b / df_b
    MS_w = SS_w / df_w  
    XF = np.divide(MS_b, MS_w)
    PvalF = np.ones((XF.shape[0],XF.shape[1]))
    for ii in range(XF.shape[0]): 
        for jj in range(XF.shape[1]):
            PvalF[ii,jj] = 1-f.cdf(XF[ii,jj], df_b, df_w)
        
    SigF = 1.0 * (PvalF < SigLevel)
    nSigS  = np.sum(SigF, axis=0)  
    nSigT  = np.sum(SigF, axis=1) 
    statOut = {'F':XF,'Pval':PvalF,'Sig':SigF,'df_b':df_b,'df_w':df_w,'SS_b':SS_b,'SS_w':SS_w,'SS_t':SS_t,'MS_b':MS_b,'MS_w':MS_w,'nSigS':nSigS,'nSigT':nSigT}  
    #'T':XT,'Pval':PvalST,'Sig':SigST} 

    return statOut #[XF, PvalF, SigF, nSig, nSigT, df_b, df_w, SS_b, SS_w, MS_b, MS_w]
 

def two_sample_permu(X0, X1, SigLevel, nSim):
    import numpy as np
    from permustat import perm_t
    XTobs = perm_t.get_T(X0, X1)
    XTm = perm_t.get_Ts(X0, X1, nSim)

    PvalST = np.zeros((XTm.shape[0], XTm.shape[1]))
    for ii in range(XTm.shape[0]):
        for tt in range(XTm.shape[1]):
            a, b = np.histogram(XTm[ii, tt,:], bins=20)
            a = np.cumsum(a)/float(nSim) 
            if XTobs[ii, tt] < 0: 
                cutID = np.where(b < XTobs[ii, tt])[0]
                if len(cutID)>0 :    
                    PvalST[ii, tt] = a[cutID[-1]]
                else:
                    PvalST[ii, tt] = 0.0    
            else :  
                cutID = np.where(b >= XTobs[ii, tt])[0]
                if len(cutID)>0 and (cutID[0]<len(a)):     
                    PvalST[ii, tt] = 1.0 - a[cutID[0]]
                else:
                    PvalST[ii, tt] = 0.0

    SigST = 1.0 * (PvalST < (SigLevel/2))
    nSigS  = np.sum(SigST, axis=0)  
    nSigT  = np.sum(SigST, axis=1)  
    statOut = {'Obs':XTobs,'Pval':PvalST,'Sig':SigST,'nSigS':nSigS,'nSigT':nSigT} #,'XTm':XTm}
    
    return statOut #[XTobs, PvalST, SigST, XTm, nSig, nSigT]


#    from scipy.stats import norm
#    Xstd = np.std(X, axis=2)
##        XT = Xmean / Xstd  
#    XT = np.sqrt(nX) * Xmean / Xstd  
##        ut = t.isf(SigLevel/2, nX-1, 0, 1)   
#    ut = norm.isf(SigLevel/2, 0, 1)   
#    
#    Tvector = XT.reshape((nSpaces*nSamples,1))
#    mu = np.mean(Tvector)
#    sig = np.std(Tvector)
#    vu = mu + ut * sig
#    vl = mu - ut * sig
#     
#    uu = (Tvector > vu ) 
#    ll = (Tvector < vl ) 
#    eee = uu + ll 
#    eee = eee.reshape((nSpaces, nSamples))
#    SigST = 1. * eee
    

    