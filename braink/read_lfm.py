import numpy as np 

def forward(lfmfilename, nE):
    nC = nE + 1 
    fd = open(lfmfilename, 'rb')
    Kori = np.fromfile(file=fd, dtype=np.dtype('>d')) #big endian 
    fd.close()
 
    nV =  Kori.shape[0] / nE
    Kori = Kori.reshape(nE, nV)
 
    H = np.eye((nC))- np.ones((nC, nC))/nC
    Kori2 = np.concatenate((Kori, np.zeros((1, nV)))) 
    K = np.dot(H, Kori2) 
    return K 
 

def lfm(lfmfilename, nE):
    nC = nE + 1 
    fd = open(lfmfilename, 'r')
    Kori = np.fromfile(file=fd, dtype=np.dtype('d')) #big endian ? 
    fd.close()
 
    nV =  Kori.shape[0] / nE
    Kori = Kori.reshape(nE, nV)
 
    H = np.eye((nC))- np.ones((nC, nC))/nC
    Kori2 = np.concatenate((Kori, np.zeros((1, nV)))) 
    K = np.dot(H, Kori2) 
    return K 

 
#nE = 256
#nV = 2219
#nC = nE + 1

#lfmfilename = "/Users/jsong/Python2.7/108/fdmForwardMatrix/fdmForwardMatrixOriented" 
#        nE = 256
#        nV = 2219
#        nC = nE + 1
##
##        lfmfilename = askopenfilename(initialdir='.', title="Select a LFM. ")
#
#        fd = open(lfmfilename, 'rb')
#        Kori = np.fromfile(file=fd, dtype=np.dtype('>d')).reshape(nE, nV)
#        fd.close()
#
#        H = np.eye((nC))- np.ones((nC, nC))/nC
#        Kori2 = np.concatenate((Kori, np.zeros((1, nV)))) 
#        K = np.dot(H, Kori2) 
#        config.K = K
#        print("LFM is loaded.")
#

