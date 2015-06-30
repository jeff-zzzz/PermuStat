# -*- coding: utf-8 -*-
"""
Created on Mon May  4 11:20:07 2015
@author: jesong1126
"""
def sden2rms(sden): 
    import numpy as np 
    nV = sden.shape[0]/3
    sden2 = sden * sden  
    rms = np.zeros(nV)
    for i in range(nV):
        rms[i] = np.sqrt(np.mean(sden2[3*i+np.arange(3)]))
    
    return rms 


#    sden = sdenI 
#    nD = sden.shape[0]
#    eyeD = np.eye(nD/3)  
#    eye3 = np.kron(eyeD, np.ones((1,3)))  
#    rms = np.sqrt(eye3 * sden2 / 3)  


