# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 10:45:29 2015
@author: jesong1126
"""

import numpy as np
nV = 1196 
nD = 3588 
nSamples = 1000
#nSamples = 1600 
#scalars = np.arange(1196.0)
#scalars /= 2.0 #1195.0 
#scalars = abs(np.random.normal(0.0, 1., (nV,1)))
#x = np.ones(1196.)
#x /= np.sqrt(3.0)
#x = np.random.normal(0.0, 1., (nV,3))
#x.reshape((nV,3))

#sdenMean = np.mean(sden, axis=2)
#sdenSigma = np.std(sden, axis=2)
#sdenT =  sdenMean / sdenSigma # * np.sqrt(44) 
#sdenMean3 = sdenMean
#sdenSigma3 = sdenSigma
#sdenT3 = sdenT

#H3 = np.kron(np.eye(nV), [1, 1, 1])
#sdenT = np.dot(H3, sdenT3)

MM = [228, 281, 466, 524, 606] #105
i=1
nSamples = 1600

f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sdenMean3.bin","r")
sdenMean3 = np.fromfile(file=f, dtype="f8").reshape((nD,nSamples))
f.close() 

f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sssTh3.bin","r")
sssTh3 = np.fromfile(file=f, dtype="f8").reshape((nD,nSamples))
f.close() 

f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sdenT3.bin","r")
sdenT3 = np.fromfile(file=f, dtype="f8").reshape((nD,nSamples))
f.close() 
 
f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sssTh.bin","r")
sssTh = np.fromfile(file=f, dtype="f8").reshape((nV,nSamples))
f.close() 

f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sdenT.bin","r")
sdenT = np.fromfile(file=f, dtype="f8").reshape((nV,nSamples))
f.close() 


#MicroMiddle = [169, 522, 571, 666, 691, 709] # 130  
#i=1
 
#f=open("/Users/jesong1126/Work/Data/VGT/VGT_130/sssTh3.bin","r")
#sssTh3 = np.fromfile(file=f, dtype="f8").reshape((nD,nSamples))
#f.close() 
#
#f=open("/Users/jesong1126/Work/Data/VGT/VGT_130/sdenT3.bin","r")
#sdenT3 = np.fromfile(file=f, dtype="f8").reshape((nD,nSamples))
#f.close() 
# 
#f=open("/Users/jesong1126/Work/Data/VGT/VGT_130/sssTh.bin","r")
#sssTh = np.fromfile(file=f, dtype="f8").reshape((nV,nSamples))
#f.close() 
#
#f=open("/Users/jesong1126/Work/Data/VGT/VGT_130/sdenT.bin","r")
#sdenT = np.fromfile(file=f, dtype="f8").reshape((nV,nSamples))
#f.close() 
#

PyEGI = egiData.PythonAdapter()
#c = PyEGI.getCondition("Triple")
PyEGI.addTriplesCondition("Triple")
c = PyEGI.getCondition("Triple")

MM = [228, 281, 466, 524, 606] #105
i=2
scalars = abs(sssTh[:, MM[i]])
x = sdenT3[:, MM[i]]

for i in range(nV):
  c.setDipoleScalar(i, scalars[i])
  c.setDipoleDirection(i, x[3*i], x[3*i+1], x[3*i+2], False)

#PyEGI.addOrientedCondition("O:test")


