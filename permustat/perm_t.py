# -*- coding: utf-8 -*-
"""
Created on Thu Mar 19 12:33:27 2015
@author: jesong1126
"""
import numpy as np

def get_T(X0, X1):
    Mit0 = np.mean(X0, axis=2)
    Mit1 = np.mean(X1, axis=2)
    nume = Mit0 - Mit1
    Sigma0 = np.var(X0, axis=2)
    Sigma1 = np.var(X1, axis=2)
#    Sigma0 = np.std(X0, axis=2)
#    Sigma1 = np.std(X1, axis=2)
#    denume = np.sqrt((Sigma0 * Sigma0)/X0.shape[2] + (Sigma1 * Sigma1)/X1.shape[2])
    denume = np.sqrt((Sigma0 /X0.shape[2]) + (Sigma1 /X1.shape[2]))
    T = np.divide(nume, denume) 
    return T
    

def get_Ts(X0, X1, nSim):
    nX0 = X0.shape[2]
    Pooled = np.concatenate((X0, X1), axis=2) 
    [nV, nSamples, nPooled]= Pooled.shape
    Ts = np.zeros((nV, nSamples, nSim))
    for i in range(nSim):
        Perm = np.random.permutation(nPooled)
        X0m = Pooled[:,:, Perm[:nX0]]
        X1m = Pooled[:,:, Perm[nX0:]]
        Ts[:,:,i] = get_T(X0m, X1m)
    return Ts
    