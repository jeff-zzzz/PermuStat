# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 16:39:11 2015
@author: jesong1126
"""

 #Kmeans.py
 #This provides the kmeans clustering EEG data given predetermined clustering size. 
 #Input  
 #  Data: EEG data, the size of matrix is the number of channels (row) by
 #  the number of samples (column). 
 #  k   : predetermined number of clusters. If Data is the ERP type, k=5 is
 #  recommanded. But if Data is the Spikes type (times course is not
 #  important), k=3 is recommended
 #Output  
 #  Tmaps   : Template Maps. The matrix is the size of the number of channels by number of clusters (k) 
 #  KmeanID : Each sample is associated with KmeanID. It is the vector of
 #  length as number of samples. 
 #  GEV : global explained variance by Tmaps

def Change(Data, k, nStable): 
#    import random   
    import numpy as np 
    
    # Read Data 
    nC, nT = Data.shape 
    GFP = np.std(Data, axis= 0) 

    # Initial Tmaps are defined as the first k samples. 
    CutId = np.append(0, np.sort(np.random.choice(nT,k-1)))
    CutId = np.append(CutId, nT)
    Tmaps = np.zeros((nC, k))
    for i in range(k):
        ithID = range(CutId[i], CutId[i+1]) 
        if ithID != [] :
            Tmaps[:,i] = np.mean(Data[:, ithID], axis=1) 
            
    # Assign each sample to one of Tmaps by finding the maximal of Correlation.
    KmeanId = np.zeros((nT,1)) 
    for i in range(k):
        KmeanId[range(CutId[i], CutId[i+1])] = i  
    KmeanId = np.array(KmeanId, dtype='i4') 

    # Correlation between Tmaps and Samples 
    Cuv = np.zeros((nT,1))  
    for j in range(nT):
        u = np.reshape(Tmaps[:,KmeanId[j]],(nC, -1) )
        v = np.reshape(Data[:,j], (nC,-1))
        u_modulus = np.sqrt((u*u).sum()) 
        v_modulus = np.sqrt((v*v).sum())
        cos_angle = np.dot(u.T,v) / u_modulus / v_modulus  
        Cuv[j] = cos_angle 
    
    # Compute the GEV for this Tmaps. 
    GEV = np.sum(GFP * Cuv * GFP * Cuv)/np.sum(GFP*GFP)/nT  
    
    # Among Kmeans clusting output of size k, find the optimal set of Tmaps
#    nStable = 100 
    CutIds = np.zeros((nStable, k+1))    
    GEVs = np.zeros((nStable,1))
    
    for m in range(nStable):      
        CutIds[m, :] = CutId
        GEVs[m,:] =  GEV 
        deathId = 1 + np.random.choice(k-1, 1) 
        CutId_ = np.delete(CutId, deathId)
        birthId0 = np.delete(np.arange(nT), CutId)
        birthId = birthId0[np.random.choice(birthId0.shape[0], 1)]
        CutIdProp = np.sort(np.concatenate((CutId_, birthId), axis=0)) 
        TmapsProp = np.zeros((nC, k))
        for i in range(k):
            ithID = range(CutIdProp[i], CutIdProp[i+1])
            if ithID != [] :
                TmapsProp[:,i] = np.mean(Data[:, ithID], axis=1) 

        KmeanIdProp = np.zeros((nT,1)) 
        for i in range(k):
            KmeanIdProp[range(CutIdProp[i], CutIdProp[i+1])] = i  
        KmeanIdProp = np.array(KmeanIdProp, dtype='i4') 

        CuvProp = np.zeros((nT,1))  
        for j in range(nT):
            u = np.reshape(TmapsProp[:,KmeanIdProp[j]],(nC, -1) )
            v = np.reshape(Data[:,j], (nC,-1))
            u_modulus = np.sqrt((u*u).sum()) 
            v_modulus = np.sqrt((v*v).sum())
            cos_angle = np.dot(u.T,v) / u_modulus / v_modulus  
            CuvProp[j] = cos_angle 

        GEVProp = np.sum(GFP * CuvProp * GFP * CuvProp)/np.sum(GFP*GFP)/nT  

        if (GEVProp > GEV):
            GEV = GEVProp   
            KmeanId = KmeanIdProp 
            Tmaps = TmapsProp 
            CutId =CutIdProp
         
#    colorkeys = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','firebrick','darkgoldenrod','purple','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00') 
#    KmeanColorId = [None] * nT 
#    for i in range(nT):
#        KmeanColorId[i] = colorkeys[KmeanId[i]] 
    
    ChangeOut = {'Tmaps':Tmaps, 'KmeanId':KmeanId, 'GEVs':GEVs, 'CutId':CutId , 'CutIds':CutIds} #,'KmeanColorId':KmeanColorId} 
    return ChangeOut 
    

