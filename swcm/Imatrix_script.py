import sys 
sys.path.append('/Users/jsong/Python3/gs3')
import numpy as np 
from braink import read_bk 
from swcm import Lcurve
from mff import *  
from matplotlib import pyplot as plt
#import Imatrix 
#import math 
#import glob

bkdir = '/Users/jesong1126/Work/Data/SEP_MEP_108/108_HeadModel'
nE = 256
K  = read_bk.lfm(bkdir+'/fdmForwardMatrixOriented', nE)  

mfffilepath ='/Users/jesong1126/Work/Data/SEP_MEP_108/SEP/' 
filePath = mfffilepath+'SEP_108_0691.blc.ave.aref.mff'

hdr = read_mff_header.read_mff_header(filePath)

hdr['nTrials']
dataEpoch = read_mff_data.read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr)
dataEpoch.shape
ts = (np.arange(0,  hdr['nSamples']) - hdr['nSamplesPre']) * (1000 / hdr['Fs'])
nSamples = dataEpoch.shape[1] 
nC, nV = K.shape

s = dataEpoch[:,:,0] 
gfp = np.std(s, axis=0)
maxGfp = np.where(gfp == max(gfp))[0]
plt.plot(gfp)

alpha = 0.01
ImatMN = Imatrix.MN(alpha, K)
H = np.eye((nC))- np.ones((nC,nC))/nC
s = np.dot(H, s) 

#mfffilename = '/Users/jsong/Python2.7/108/SEP/SEP108_RT.bin' 
#nC = nE+1
#fd = open(mfffilename , 'rb')
#data = np.fromfile(file=fd, dtype=np.dtype('f4'))  
#fd.close()
#nSamples = data.shape[0]/nC
#sT = np.array(data.reshape(nSamples, nC))
#s = sT.T
#s2 = np.multiply(s, s)
#gfp = s2.sum(axis=0)

maxi = (gfp == max(gfp)).nonzero()
phi = s[:, maxi[0]] 

plt.plot(phi)  
plt.xlim((0, nC))
plt.xlabel('channel ')
plt.ylabel('mV')

LcMN = Lcurve.LcurveMN(phi, K) 
LcCSL = Lcurve.LcurveCSL(phi, K, bkdir) 
 
 
##------------------------------------------------------------------------------ 
 
#import sys 
#sys.path.append('/Users/jsong/Python27')

from __future__ import print_function

from IPython.parallel import Client
from matplotlib import pyplot as plt
import numpy
#import numpy as np
import random
from timeit import default_timer as clock
 
##------------------------------------------------------------------------------ 
nV = 2229
nTime = 150 
nSim = 15

filename = '/Users/jesong1126/Python27/Data/SEP107_RT1_369_sden.bin' 
with open(filename , 'r') as fd:
    data = numpy.fromfile(file=fd, dtype=np.dtype('d'))  
    fd.close()

nSamples = data.shape[0]/nV
sdenT = numpy.array(data.reshape(nSamples, nV))
sden = sdenT.T
nX0 = nSamples/nTime
X0 = numpy.zeros((nV, nTime, nX0))
for i in numpy.arange(nX0):
    X0[:,:,i] = sden[:,nTime*i+np.arange(nTime)]

filename = '/Users/jesong1126/Python27/Data/SEP107_LT1_364_sden.bin' 
with open(filename , 'r') as fd: 
    data = numpy.fromfile(file=fd, dtype=np.dtype('d'))  
    fd.close()

nSamples = data.shape[0]/nV
sdenT = numpy.array(data.reshape(nSamples, nV))
sden = sdenT.T
nX1 = nSamples/nTime
X1 = numpy.zeros((nV, nTime, nX1))
for i in numpy.arange(nX1):
    X1[:,:,i] = sden[:,nTime*i+np.arange(nTime)]
    
del data 
del sden 
del sdenT 
del nSamples 

## 
Pooled = numpy.concatenate((X0, X1), axis=2) 
nPooled = nX0 + nX1 

def get_Ts(Pooled, nX0, nSim):
    [nV, nTime, nPooled]= Pooled.shape
    Ts = zeros((nV, nTime, nSim))
    items = range(nPooled)
    for i in range(nSim):
        random.shuffle(items)
        X0m = Pooled[:,:, items[:nX0]]
        X1m = Pooled[:,:, items[nX0:]]
        Mit0m = numpy.mean(X0m, axis=2)
        Mit1m = numpy.mean(X1m, axis=2)
        Sigma0m = numpy.std(X0m, axis=2)
        Sigma1m = numpy.std(X1m, axis=2)
        nume = Mit0m - Mit1m
        denume = numpy.sqrt((Sigma0m * Sigma0m)/nX0 + (Sigma1m * Sigma1m)/nX1)
        Ts[:,:,i] = numpy.divide(nume, denume) 
    return Ts

##
# Connect to the IPython cluster
c = Client()

# the number of engines
n = len(c)
dview = c[:]
dview.block=True
print("%i engines "%n)

with dview.sync_imports():
    import numpy
    import random 

dview.apply(get_Ts, Pooled, nX0, nSim)

