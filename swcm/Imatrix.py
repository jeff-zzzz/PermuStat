import numpy as np
from braink import read_bk
import glob

#K = read_bk.lfm('/Users/jsong/Python2.7/108/108_HeadModel_MRI/fdmForwardMatrixOriented', 256)

##------------------------------------------------------------------------------
def MN(alpha, K):
    
    Ne, Nd = K.shape 
    H = np.eye((Ne))- np.ones((Ne, Ne))/Ne
    K = np.dot(H, K) 
    svdK = np.linalg.svd(K, compute_uv=False)    
    
    KtK = np.dot(K.T, K)
    M = KtK + alpha * svdK[0] * np.identity(Nd)  
    invM = np.linalg.pinv(M) 
    Imatrix = np.dot(invM, K.T) 

    #KKt = np.dot(K, K.T) 
    #diagave = np.mean(np.diag(KKt)) 
    #M = KKt + alpha * diagave * H  
    #invM = np.linalg.pinv(M) 
    #Imatrix = np.dot(K.T, invM) 

    return Imatrix 

##------------------------------------------------------------------------------
def sMN(alpha, K):
    
    Ne, Nd = K.shape 
    H = np.eye((Ne))- np.ones((Ne, Ne))/Ne
    K = np.dot(H, K) 

    KKt = np.dot(K, K.T) 
    diagave = np.mean(np.diag(KKt)) 
    M = KKt + alpha * diagave * H  
    invM = np.linalg.pinv(M) 
    ImatMN = np.dot(K.T, invM) 
    Cjhat = np.dot(ImatMN, K)
    Imatrix = np.dot(np.diag(1/np.sqrt(np.diag(Cjhat))), ImatMN)

    return Imatrix 

##------------------------------------------------------------------------------
def LORETA(alpha, K, bkdir, MaxNeighbor):    
  
    Ne, Nd = K.shape 
    H = np.eye((Ne))- np.ones((Ne, Ne))/Ne
    K = np.dot(H, K) 
    
    for bkdfilename1 in glob.iglob(bkdir+"/*.Dipoles_1_*"):         
        D1 = read_bk.bkd(bkdfilename1)
    for bkdfilename2 in glob.iglob(bkdir+"/*.Dipoles_2_*"): 
        D2 = read_bk.bkd(bkdfilename2)

    Nd = D1['ndipoles'] + D2['ndipoles']

    Mat1 = np.concatenate((np.matrix(D1['dipole_location_index']['x']),np.matrix(D1['dipole_location_index']['y']), np.matrix(D1['dipole_location_index']['z'])), axis=0)
    Mat2 = np.concatenate((np.matrix(D2['dipole_location_index']['x']),np.matrix(D2['dipole_location_index']['y']), np.matrix(D2['dipole_location_index']['z'])), axis=0)
    VoxCoord  = np.array(np.concatenate((Mat1, Mat2), axis=1))

    Vox2 = VoxCoord*VoxCoord 
    #vcnorms = np.matrix(np.sum(VoxCoord*VoxCoord, axis=0))
    #vcnorms = np.matrix(Vox2.sum(axis=0)) 
    vcnorms = np.matrix(np.sum(Vox2, axis=0)) 
    tridist = vcnorms - 2 * np.dot( VoxCoord.T, VoxCoord) + vcnorms.T
    tridist= np.sqrt(tridist)

    Neighbor = np.array(tridist < MaxNeighbor, dtype=float)  

    B = np.identity(Nd)
    for i in range(Nd):   
        Ni = np.sum(Neighbor[i,:])-1   
        if Ni > 0:
            B[i,:] = -Neighbor[i,:]/Ni 
            B[i,i] = 1 
  
    BigK = np.concatenate((K, alpha * B), axis=0)
    KtK = np.dot(BigK.T, BigK)
    
    T = np.dot(np.linalg.inv(KtK),  BigK.T) 
    Imatrix = T[:, :Ne]

    return Imatrix 

##------------------------------------------------------------------------------

def CSL(alphaS, alphaL, K, bkdir): 
    
    Ne, Nd = K.shape 
    H = np.eye((Ne))- np.ones((Ne, Ne))/Ne
    K = np.dot(H, K) 

    for dsafilename1 in glob.iglob(bkdir+"/*.Adjacents_1_*"):         
        A_local1 = read_bk.dsa(dsafilename1)

    for dsafilename2 in glob.iglob(bkdir+"/*.Adjacents_2_*"): 
        A_local2 = read_bk.dsa(dsafilename2)

    #dsafilename1 = bkdir+prefix+'Adjacents_1_'+scenario+'.dsa'
    #dsafilename2 = bkdir+prefix+'Adjacents_2_'+scenario+'.dsa'
    
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

    svdK = np.linalg.svd(K, compute_uv=False)
    svdB = np.linalg.svd(B, compute_uv=False)
    condNum = svdK[0]/svdB[0]
    
    #BigK = np.concatenate((K, alpha * B), axis=0)
    #KtK = np.dot(BigK.T, BigK)

    I = np.identity(Nd)
    BigK = np.concatenate((K, alphaS*svdK[0]*I,  alphaL*condNum*B), axis=0)
    KtK = np.dot(BigK.T, BigK)

    T = np.dot(np.linalg.inv(KtK),  BigK.T) 
    Imatrix = T[:, :Ne]
    
    return Imatrix 


##------------------------------------------------------------------------------
def tCSL(alphaS, alphaL, alphaT, K, bkdir):#, prefix, scenario):
    
    Ne, Nd = K.shape 
    H = np.eye((Ne))- np.ones((Ne, Ne))/Ne
    K = np.dot(H, K) 
    
    for dsafilename1 in glob.iglob(bkdir+"/*.Adjacents_1_*"): 
        A_local1 = read_bk.dsa(dsafilename1)

    for dsafilename2 in glob.iglob(bkdir+"/*.Adjacents_2_*"):  
        A_local2 = read_bk.dsa(dsafilename2)

    #dsafilename1 = bkdir+prefix+'Adjacents_1_'+scenario+'.dsa'
    #A_local1 = read_dsa.read_dsa(dsafilename1)
    #dsafilename2 = bkdir+prefix+'Adjacents_2_'+scenario+'.dsa'
    #A_local2 = read_dsa.read_dsa(dsafilename2)
    
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
            
    Nd3 = 3 * Nd
    K3 = np.concatenate((np.zeros((Ne,Nd)),  K,  np.zeros((Ne,Nd))), axis=1) 
    B3 = np.kron(np.identity(3), B)
    I3 = np.identity(Nd3)
    Bt = np.concatenate((-np.identity(Nd),  2*np.identity(Nd), -np.identity(Nd)), axis=1)  
    BigK = np.concatenate((K3, alphaS * I3,  alphaL * B3,  alphaT * Bt), axis=0)
    KtK = np.dot(BigK.T, BigK)  
    T = np.dot(np.linalg.inv(KtK),  BigK.T) 
    Imatrix = T[Nd:(2*Nd), :Ne]

    return Imatrix
    
##------------------------------------------------------------------------------




##------------------------------------------------------------------------------
