import numpy as np 

def read_dsa(dsafilename):
    fd = open(dsafilename, 'rb')
    A_dsa = np.fromfile(file=fd, dtype=np.dtype('i4'))
    fd.close()

    nV = int(np.sqrt(A_dsa.shape[0]))
    A_dsa = A_dsa.reshape(nV, nV)

    return A_dsa 


#    dsafilename = "/Users/jsong/Python2.7/108/fdmForwardMatrix/108_01_2000.Adjacents_1_MRI.dsa" 
#    fd = open(dsafilename, 'rb')
#    A_dsa = np.fromfile(file=fd, dtype=np.dtype('i4'))
#    fd.close()
#
#    nV = int(np.sqrt(A_dsa.shape[0]))
#    A_dsa = A_dsa.reshape(nV, nV)
 
