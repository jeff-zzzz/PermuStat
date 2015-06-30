# -*- coding: utf-8 -*-
"""
Created on Fri Jan 30 10:45:29 2015
@author: jesong1126
"""
#import os 
#os.chdir("/Users/jesong1126/Work/Data/EAV/2941_HM") 

import numpy as np
import matplotlib.pyplot as plt 

## 
from braink import read_bk
lbkd = read_bk.bkd('/Users/jesong1126/Python27/GeoPy/data/108_MIE/108_3.Dipoles_left_MRI.bkd')
rbkd = read_bk.bkd('/Users/jesong1126/Python27/GeoPy/data/108_MIE/108_3.Dipoles_right_MRI.bkd')
lbkd['ndipoles'] + rbkd['ndipoles']

lcoord = np.array((lbkd['dipole_location_index']['x'], lbkd['dipole_location_index']['y'], lbkd['dipole_location_index']['z'])).T
rcoord = np.array((rbkd['dipole_location_index']['x'], rbkd['dipole_location_index']['y'], rbkd['dipole_location_index']['z'])).T
VoxCoord = np.concatenate((lcoord, rcoord), axis=0)
nVL = lbkd['ndipoles'] 
nVR = rbkd['ndipoles'] 
mmin = np.min(VoxCoord, axis=0)
mmax = np.max(VoxCoord, axis=0)

fig, ax = plt.subplots(2,2) 
ax[0,0].scatter(VoxCoord[:,0], VoxCoord[:,1], s=100, cmap='hot')   #c=nSigVoxGrpI, 
ax[0,0].set_axis_off()
ax[0,0].set_aspect('equal')    
ax[0,0].set_xlim([mmin[0]-10, mmax[0]+10])
ax[0,0].set_ylim([mmin[1]-10, mmax[1]+10])

ax[0,1].scatter(VoxCoord[:,0], VoxCoord[:,2], s=100, cmap='hot')  
ax[0,1].set_aspect('equal')
ax[0,1].set_axis_off()
ax[0,1].set_xlim([mmin[0]-10, mmax[0]+10])
ax[0,1].set_ylim([mmin[2]-10, mmax[2]+10])

ax[1,0].scatter(-lcoord[:,1], lcoord[:,2],   s=100, cmap='hot')    
ax[1,0].set_aspect('equal')
ax[1,0].set_axis_off()    
ax[1,0].set_xlim([-mmax[1]-10, -mmin[1]+10 ])
ax[1,0].set_ylim([mmin[2]-10, mmax[2]+10])

ax[1,1].scatter(rcoord[:,1], rcoord[:,2],   s=100, cmap='hot')   
ax[1,1].set_aspect('equal')
ax[1,1].set_axis_off()
ax[1,1].set_xlim([mmin[1]-10, mmax[1]+10])
ax[1,1].set_ylim([mmin[2]-10, mmax[2]+10])
  
  
f=open("/Users/jesong1126/Python27/GeoPy/data/108_HM/Leadfield.lfm","r")
K0 = np.fromfile(file=f, dtype="f8") #.reshape((nD,nSamples))
f.close()

nE = 256 
nV = K0.shape[0] / nE 

K = K0.reshape((nE,nV))

fig, ax = plt.subplots(1,1)
ax.plot(K[:,0])
plt.show()


i=0
VoxCoord[i,:]

fig, ax = plt.subplots(2,2) 
ax[0,0].scatter(VoxCoord[:,0], VoxCoord[:,1], s=10, color='red')   #c=nSigVoxGrpI, 
ax[0,0].scatter(VoxCoord[i,0], VoxCoord[i,1], s=100, color='blue')   #c=nSigVoxGrpI, 
ax[0,0].set_axis_off()
ax[0,0].set_aspect('equal')    
ax[0,0].set_xlim([mmin[0]-10, mmax[0]+10])
ax[0,0].set_ylim([mmin[1]-10, mmax[1]+10])

ax[0,1].scatter(VoxCoord[:,0], VoxCoord[:,2], s=10, color='red')  
ax[0,1].scatter(VoxCoord[i,0], VoxCoord[i,2], s=100, color='blue')  
ax[0,1].set_aspect('equal')
ax[0,1].set_axis_off()
ax[0,1].set_xlim([mmin[0]-10, mmax[0]+10])
ax[0,1].set_ylim([mmin[2]-10, mmax[2]+10])

ax[1,0].scatter(-lcoord[:,1], lcoord[:,2],   s=10, color='red')    
ax[1,0].scatter(-lcoord[i,1], lcoord[i,2],   s=100, color='blue')    
ax[1,0].set_aspect('equal')
ax[1,0].set_axis_off()    
ax[1,0].set_xlim([-mmax[1]-10, -mmin[1]+10 ])
ax[1,0].set_ylim([mmin[2]-10, mmax[2]+10])

ax[1,1].scatter(rcoord[:,1], rcoord[:,2],   s=10, color='red')   
#ax[1,1].scatter(rcoord[i,1], rcoord[i,2],   s=100, color='blue')   
ax[1,1].set_aspect('equal')
ax[1,1].set_axis_off()
ax[1,1].set_xlim([mmin[1]-10, mmax[1]+10])
ax[1,1].set_ylim([mmin[2]-10, mmax[2]+10])
  
  
#c = PyEGI.getCondition("Oriented") 
#eegData = numpy.arange(256.0)
#c.clearSensorData()
#for i in range(eegData.shape[0]):
#  c.addSensorData(eegData[i])

PyEGI = egiData.PythonAdapter()
PyEGI.addOrientedCondition("Oriented")
c = PyEGI.getCondition("Oriented")

j = 968
scalar = range(nV)
scalar[j] = 1
x = np.ones(3*nV)

for i in range(nV):
  c.setDipoleScalar(i, scalar[i])
  c.setDipoleDirection(i, x[3*i], x[3*i+1], x[3*i+2], False)

for i in range(nE):
  c.addSensorData(K[i,j]) 

#>>> for i in range(nV):
#...  c.setDipoleScalar(i, scalar[i])
#...  c.setDipoleDirection(i,x[3*i], x[3*i+1], x[3*i+2], False) 
#... 
#Traceback (most recent call last):
#  File "<string>", line 2, in <module>
#Boost.Python.ArgumentError: Python argument types in
#    Condition.setDipoleScalar(Condition, int, numpy.ndarray)
#did not match C++ signature:
#    setDipoleScalar(Condition {lvalue}, int, float)
#        
#nV = 1196 
#nD = 3588 
#nSamples = 1000
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

#MM = [228, 281, 466, 524, 606] #105
#i=1
#nSamples = 1600
#
#f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sdenMean3.bin","r")
#sdenMean3 = np.fromfile(file=f, dtype="f8").reshape((nD,nSamples))
#f.close() 
#
#f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sssTh3.bin","r")
#sssTh3 = np.fromfile(file=f, dtype="f8").reshape((nD,nSamples))
#f.close() 
#
#f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sdenT3.bin","r")
#sdenT3 = np.fromfile(file=f, dtype="f8").reshape((nD,nSamples))
#f.close() 
# 
#f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sssTh.bin","r")
#sssTh = np.fromfile(file=f, dtype="f8").reshape((nV,nSamples))
#f.close() 
#
#f=open("/Users/jesong1126/Work/Data/VGT/VGT105_ForErik/sdenT.bin","r")
#sdenT = np.fromfile(file=f, dtype="f8").reshape((nV,nSamples))
#f.close() 


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
##f.close() 
##
#
#PyEGI = egiData.PythonAdapter()
##c = PyEGI.getCondition("Triple")
#PyEGI.addTriplesCondition("Triple")
#c = PyEGI.getCondition("Triple")
#
#MM = [228, 281, 466, 524, 606] #105
#i=2
#scalars = abs(sssTh[:, MM[i]])
#x = sdenT3[:, MM[i]]
#
#for i in range(nV):
#  c.setDipoleScalar(i, scalars[i])
#  c.setDipoleDirection(i, x[3*i], x[3*i+1], x[3*i+2], False)
#
##PyEGI.addOrientedCondition("O:test")
