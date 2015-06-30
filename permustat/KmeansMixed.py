

def KmeansMixed(Data, nClusters, nSim): 

    import random 
    import numpy as np

    nC, nT = Data.shape  
    GFP = np.std(Data, axis = 0) 
    k = nClusters - 1 # number of inner partition edges  
    
    PartitionId = np.concatenate([[0],np.sort(random.sample(range(nT), k)),[nT]], axis=0)  
    KmeanId = np.array([]) 
    for i in range(nClusters):
        KmeanId = np.concatenate((KmeanId, i*np.ones((PartitionId[i+1]-PartitionId[i]))))
    
    KmeanId = np.array(KmeanId, 'i4')
    
    Tmaps = np.zeros((nC, nClusters))
    for i in range(nClusters):
        DataI = Data[:, range(PartitionId[i], PartitionId[i+1])]
        DataI = np.reshape(DataI, (nC,-1)) 
        Tmaps[:,i] = np.mean(DataI, axis=1) 
        
    # Correlation between Tmaps and Samples 
    Cuv = np.array([]) 
    for i in range(nClusters) :
        DataI = Data[:, range(PartitionId[i], PartitionId[i+1])]
        DataI = np.reshape(DataI, (nC,-1))  
        CuvI = np.zeros((DataI.shape[1]))
        u = Tmaps[:,i] 
        u_modulus = np.sqrt((u*u).sum()) 
        for j in range(DataI.shape[1]): 
            v = DataI[:,j]
            v_modulus = np.sqrt((v*v).sum())
            cos_angle= np.dot(u,v) / u_modulus / v_modulus  
            CuvI[j] = cos_angle
        Cuv = np.concatenate((Cuv, CuvI))   
       
    GEV = sum(GFP * Cuv * GFP * Cuv)/sum(GFP * GFP) 
    
    GEVs = np.zeros((nSim))         
    for ii in range(nSim): 
                
        PartitionIdProp = np.concatenate([[0],np.sort(random.sample(range(nT), k)),[nT]], axis=0)  

        KmeanIdProp = np.array([]) 
        for i in range(nClusters):
            KmeanIdProp = np.concatenate((KmeanIdProp, i*np.ones((PartitionIdProp[i+1]-PartitionIdProp[i])))) 
                
        KmeanIdProp = np.array(KmeanIdProp, 'i4')
        
        TmapsProp = np.zeros((nC, nClusters))
        for i in range(nClusters):
            DataI = Data[:, range(PartitionIdProp[i], PartitionIdProp[i+1])] 
            DataI = np.reshape(DataI, (nC,-1))          
            TmapsProp[:,i] = np.mean(DataI, axis=1)
        
        # Correlation between Tmaps and Samples 
        CuvProp = np.array([])  
        for i in range(k+1) :
            DataI = Data[:, range(PartitionIdProp[i], PartitionIdProp[i+1])]
            DataI = np.reshape(DataI, (nC,-1))
            CuvI = np.zeros((DataI.shape[1]))
            u = Tmaps[:,i] 
            u_modulus = np.sqrt((u*u).sum()) 
            for j in range(DataI.shape[1]): 
                v = DataI[:,j]
                v_modulus = np.sqrt((v*v).sum())
                cos_angle= np.dot(u,v) / u_modulus / v_modulus  
                CuvI[j] = cos_angle
            CuvProp = np.concatenate((CuvProp, CuvI))   
       
        GEVProp = sum(GFP * CuvProp * GFP * CuvProp)/sum(GFP * GFP)  

        if GEVProp > GEV: 
            GEV = GEVProp
            PartitionId = PartitionIdProp
            KmeanId = KmeanIdProp
            Tmaps = TmapsProp

        GEVs[ii] = GEV
    
        
    colorkeys = ('b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','firebrick','darkgoldenrod','purple','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00') 
    KmeanColorId = [None] * nT 
    for i in range(nT):
        KmeanColorId[i] = colorkeys[KmeanId[i]]
    
            
    KmeansRJOut = {'Tmaps':Tmaps, 'KmeanId':KmeanId, 'Cuv':Cuv,'GEV':GEV, 'GEVs':GEVs,'nClusters':nClusters, 
    'PartitionId':PartitionId, 'KmeanColorId':KmeanColorId} 
    return KmeansRJOut 

