# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 14:46:38 2015

@author: jesong1126
"""
import numpy as np   
import matplotlib.pyplot as plt 
#from braink import read_lfm  

nE = 128    
nC = nE + 1 

fname ='/Users/jesong1126/Work/Support/Atlas_LF_Compare/AMale_40_128_trip.lfm'
fd = open(fname, 'r')
Kvec = np.fromfile(file=fd, dtype=np.dtype('d')) #big endian ? 
fd.close()
nV =  Kvec.shape[0] / nE
K = Kvec.reshape(nE, nV)
 
fig = plt.figure()
plt.plot(K[:,0])
plt.show()



 
