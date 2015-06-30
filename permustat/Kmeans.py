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

def Kmeans(Data, k, nStable): 
    import random   
    import numpy as np 
    
    # Read Data 
    nC, nT = Data.shape 
    GFP = np.std(Data, axis= 0) 

    # Initial Tmaps are defined as the first k samples. 
    TmapId = np.sort(random.sample(range(nT), k))  
    Tmaps = Data[:, TmapId]  

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
            
    # Assign each sample to one of Tmaps by finding the maximal of Correlation.
    KmeanId = np.zeros((nT,), 'int32') 
    for j in range(nT):
        KmeanId[j] = int(np.where((Cuv[:,j])== max(Cuv[:,j]))[0])  
        
    # Compute the average within the same cluster.
    C_UTmaps = np.zeros((nT,))
    for j in range(nT):
        C_UTmaps[j] = Cuv[KmeanId[j], j]  
        
    # Compute the GEV for this Tmaps. 
    GEV = sum(GFP * C_UTmaps * GFP * C_UTmaps)/sum(GFP*GFP)  
    for i in range(k) :
        ithID = np.where(KmeanId == i)[0] 
        Tmaps[:,i] = np.mean(Data[:, ithID], axis=1) 
    
    # Among Kmeans clusting output of size k, find the optimal set of Tmaps
#    nStable = 100 
    GEVs = np.zeros((nStable,))
    for m in range(nStable):        
        TmapIdProp = np.sort(random.sample(range(nT), k))  
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
                             
        KmeanIdProp = np.zeros((nT,), 'int32') 
        for j in range(nT):
            KmeanIdProp[j] = int(np.where((CuvProp[:,j])== max(CuvProp[:,j]))[0])              

        C_UTmapsProp = np.zeros((nT,))
        for j in range(nT):
            C_UTmapsProp[j] = CuvProp[KmeanIdProp[j], j]  

        GEVProp = sum(GFP * C_UTmapsProp * GFP * C_UTmapsProp)/sum(GFP*GFP)  
        for i in range(k) :
            ithID = np.where(KmeanIdProp == i)[0] 
            TmapsProp[:,i] = np.mean(Data[:, ithID], axis=1) 

        if (GEVProp > GEV):
            GEV = GEVProp  
            KmeanId = KmeanIdProp 
            C_UTmaps = C_UTmapsProp
            Tmaps = TmapsProp 
         
        GEVs[m] =  GEV 

    colorkeys = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','firebrick','darkgoldenrod','purple','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00') 
    KmeanColorId = [None] * nT 
    for i in range(nT):
        KmeanColorId[i] = colorkeys[KmeanId[i]]
    
    
    KmeansOut = {'Tmaps':Tmaps, 'KmeanId':KmeanId, 'C_UTmaps':C_UTmaps,'GEV':GEV,'GEVs':GEVs,'KmeanColorId':KmeanColorId} 
    return KmeansOut 
    

