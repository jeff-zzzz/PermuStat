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

def Kmeans(Data, k, nStable, Thred): 
    import numpy as np 
    import scipy.stats 
    
    # Read Data 
    nC, nT = Data.shape 
    GFP = np.std(Data, axis= 0) 

    # Initial Tmaps are defined as the first k samples. 
    TmapId = np.sort(np.random.choice(nT,k))
    Tmaps = Data[:, TmapId]  

    # Correlation between Tmaps and Samples 
    Cuv = np.zeros((k, nT))  
    for j in range(nT):
        for i in range(k):
            u = Tmaps[:,i] 
            v = Data[:,j]
            u_modulus = np.sqrt((u*u).sum()) 
            v_modulus = np.sqrt((v*v).sum())
            cos_angle= np.dot(u,v) / u_modulus / v_modulus  
            Cuv[i, j] = cos_angle      
    
    DISSk = 1-Cuv 

    # Assign each sample to one of Tmaps by finding the maximal of Correlation.
    KmeanId = np.zeros(nT) 
    for j in range(nT):
        KmeanId[j] = np.where((DISSk[:,j])== min(DISSk[:,j]))[0] + 1 
 
    DISS = np.zeros(nT)
    for i in range(nT):
        DISS[i] = DISSk[KmeanId[i]-1,i] 
    
    # Compute the average within the same cluster.
    C_UTmaps = np.zeros(nT)
    for j in range(nT):
        C_UTmaps[j] = Cuv[KmeanId[j]-1, j]  
        
    # Compute the GEV for this Tmaps. 
    GEV = np.sum(GFP * C_UTmaps * GFP * C_UTmaps)/np.sum(GFP*GFP)/nT  
    
    # Among Kmeans clusting output of size k, find the optimal set of Tmaps
    TmapIds = np.zeros((nStable, k))
    
    GEVs = np.zeros(nStable)
    for m in range(nStable):      
        TmapIds[m, :] = TmapId
        GEVs[m] =  GEV 
        deathId = np.random.choice(k, 1) 
        TmapId_ = np.delete(TmapId, deathId)
        birthId0 = np.delete(np.arange(nT), TmapId)
        birthId = birthId0[np.random.choice(birthId0.shape[0], 1)]
        TmapIdProp = np.sort(np.concatenate((TmapId_, birthId), axis=0)) 
        TmapsProp = Data[:, TmapIdProp]  
    
        CuvProp = np.zeros((k, nT)) 
        for j in range(nT):
            for i in range(k):
                u = TmapsProp[:,i] 
                v = Data[:,j]
                u_modulus = np.sqrt((u*u).sum()) 
                v_modulus = np.sqrt((v*v).sum())
                cos_angle= np.dot(u,v) / u_modulus / v_modulus #-> cosine of the angle
                CuvProp[i, j] = cos_angle 
        
        DISSkProp = 1-CuvProp
                             
        KmeanIdProp = np.zeros(nT) 
        for j in range(nT):
            KmeanIdProp[j] = int(np.where((DISSkProp[:,j])== min(DISSkProp[:,j]))[0]) + 1    

        C_UTmapsProp = np.zeros(nT)
        for j in range(nT):
            C_UTmapsProp[j] = CuvProp[KmeanIdProp[j]-1, j]  

        DISSProp = np.zeros(nT)
        for i in range(nT):
            DISSProp[i] = DISSkProp[KmeanIdProp[i]-1,i]  
        
        GEVProp = np.sum(GFP * C_UTmapsProp * GFP * C_UTmapsProp)/np.sum(GFP*GFP)/nT  

        if (GEVProp > GEV):
            GEV = GEVProp   
            KmeanId = KmeanIdProp 
            C_UTmaps = C_UTmapsProp
            Tmaps = TmapsProp 
            TmapId =TmapIdProp
            DISS = DISSProp

    KmeanId = np.array(KmeanId, dtype='i4')
        
    SD = np.zeros(k)
    for i in range(k):
        Ith = np.where(KmeanId == (i+1))[0]
        SD[i] = np.std(DISS[Ith])

    NormIsf = scipy.stats.norm.isf((1.0-Thred/100.0)/2.0) 
    CInt = NormIsf * SD #; AngInt = np.arccos(1-CInt)

    KmeanId1 = -np.zeros(nT)
    for i in range(k):
        Ith = np.where(KmeanId == (i+1))[0]  
        KmeanId1[Ith[DISS[Ith] < CInt[i]]] = i+1         
    KmeanId1 = np.array(KmeanId1, dtype='i4')       
    
    # Generate Tmaps 
    Tmaps1 = np.zeros((nC, k))
    for i in range(k): 
        Ith = np.where(KmeanId1 == (i+1))[0] 
        Tmaps1[:,i] = np.mean(Data[:,Ith], axis=1)
    
    KmeanId2 = -np.zeros(nT) 
    for i in range(k):
        Ith = np.where(KmeanId1 == (i+1))[0]
        temp = np.where(np.diff(Ith) > 1)[0]
        temp = np.insert(temp,0,0)
        temp = np.append(temp, Ith.shape[0])
        if len(temp) > 2 :
            templen = np.diff(temp) 
            longest = np.where(templen == max(templen))[0][0]
            KmeanId2[Ith[np.arange(temp[longest]+1, temp[longest+1])]] = i+1
        else:         
            KmeanId2[Ith] = i+1
    KmeanId2 = np.array(KmeanId2, dtype='i4')
    
    # Generate Tmaps 
    Tmaps2 = np.zeros((nC, k))
    for i in range(k): 
        Ith = np.where(KmeanId2 == (i+1))[0] 
        Tmaps2[:,i] = np.mean(Data[:,Ith], axis=1)
    
    startTP = np.zeros(k) # np.array ([np.arange(k+1), np.zeros(k+1)]).reshape(2, k+1)
    for i in range(k):
       startTP[i] = np.where(KmeanId2==(i+1))[0][0]
    rankTmaps = np.argsort(startTP)
    
    KmeanId3 = -np.zeros(nT)
    for i in range(k):
        KmeanId3[KmeanId2==(rankTmaps[i]+1)] = i+1
    KmeanId3 = np.array(KmeanId3, dtype='i4')
    
    Tmaps3 = np.zeros((nC, k))
    for i in range(k): 
        Ith = np.where(KmeanId3 == (i+1))[0] 
        Tmaps3[:,i] = np.mean(Data[:,Ith], axis=1)
    

    KmeansOut = {'Tmaps':Tmaps,'KmeanId':KmeanId, 'Tmaps1':Tmaps1,'KmeanId1':KmeanId1, 'Tmaps2':Tmaps2,'KmeanId2':KmeanId2, 'Tmaps3':Tmaps3,'KmeanId3':KmeanId3,'DISS':DISS,'TmapId':TmapId,'GEVs':GEVs}#,'C_UTmaps':C_UTmaps}#,'GEV':GEV}#,'TmapIds':TmapIds}
    
    return KmeansOut 
    
#    colorkeys = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','firebrick','darkgoldenrod','purple','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00') 
#    KmeanColorId = [None] * nT 
#    for i in range(nT):
#        KmeanColorId[i] = colorkeys[KmeanId[i]] 
#,'KmeanColorId':KmeanColorId} 

