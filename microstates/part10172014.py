k0 = 5 
nSim = 1000 

for sIdx in range(8): 
    
    Data = s_blc[:,:,sIdx] 
    blc_gfp = s_blc_GFP[:,sIdx] 

    plt.plot(msSamples,Data.T) 
    plt.xlim([msSamples[0], msSamples[-1]])
    plt.savefig('VGT%sButterfly.png' %SubjId[sIdx])
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
    plt.savefig('VGT%sKmeanInitial.png' %SubjId[sIdx])
    plt.close()
        
    StableLegnth = 50 # ms 
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
    #plt.title(('Stable %d ms ' % (StableLegnth)))
    for ii in range(nMS):
        tempii = range(BeginMS[ii], EndMS[ii])     
        plt.vlines(msSamples[tempii], [0], blc_gfp[tempii], color=colorkeys[ii], alpha= 0.2)    
        plt.text(msSamples[tempii[round(len(tempii)/2)]], 1, ('%d' % ii), fontsize=15)    
        
    plt.savefig('VGT%sGfpMicroState.png' %SubjId[sIdx])
    plt.close() 

            
  
    fig,axes = plt.subplots(1, nMS, figsize=(3*nMS,3))#, sharex=True)
    for ii in range(nMS):
        axes[ii].scatter(x2, y2, c=Microstate[:,ii], s=30, cmap=plt.get_cmap('seismic'), alpha= .5)
        axes[ii].set_alpha(0.75)
        axes[ii].set_title(('MicroState %d' % ii))
        axes[ii].set_xticks([])
        axes[ii].set_yticks([])

    plt.savefig('VGT%sMicroStateTMaps.png' %SubjId[sIdx])
    plt.close() 

#plt.show()

  
    