# -*- coding: utf-8 -*-
"""
Created on Wed Jun  3 13:24:49 2015 
@author: jesong1126
"""
 
import numpy as np   
import matplotlib.pyplot as plt
from braink import read_lfm, read_bk

nE = 256
#K  = read_lfm.lfm('/Users/jesong1126/Python27/GeoPy/data/108_HM/Leadfield.lfm', nE)  # OK!!! left=1 and right=2                                 
K  = read_lfm.lfm('/Users/jesong1126/Work/EAV_MIE_BrainK/r842_ihm_test1/r842_108_LFM/r842_108_output/r842_108_ori.lfm', nE)        
#K  = read_lfm.lfm('/Users/jesong1126/Work/Data/EAV/2941_HM/Leadfield.lfm', nE)        

#/Users/jesong1126/Work/EAV_MIE_BrainK/r842_ihm_test1/r842_108_ori.lfm/
 
nD = K.shape[1] # 2383

import os
os.system("open /Applications/EAV/Mimir.app")
from RabbitMQ import Connection 

i=83; DipI= np.zeros(nD); DipI[i] = 1 
plt.plot(K[:,i]); plt.show()


c = Connection.Connection('localhost')
c.connect()
c.sendOrientedData(DipI, "Ori", "4:00pm")
c.sendEEGData(K[:, i], "EEG", "")
c.disconnect()       


Dip1 = read_bk.bkd('/Users/jesong1126/Work/EAV_MIE_BrainK/r842_ihm_test1/r842_ihm_test1.Dipoles_1_MRI.bkd')
Dip2 = read_bk.bkd('/Users/jesong1126/Work/EAV_MIE_BrainK/r842_ihm_test1/r842_ihm_test1.Dipoles_2_MRI.bkd')

d1 = Dip1['dipole_location_index']; d2 = Dip2['dipole_location_index']; 
Dip_x = np.concatenate((d1['x'], d2['x']), axis=0) 
Dip_y = np.concatenate((d1['y'], d2['y']), axis=0) 
Dip_z = np.concatenate((d1['z'], d2['z']), axis=0) 

Dip_xyz = np.concatenate((Dip_x, Dip_y, Dip_z), axis=0)
Dip_xyz = Dip_xyz.reshape(3,nD)
Dip_xyz = Dip_xyz.T 

Dip_xyz[i,:]
 

print(np.min(d1['x']), np.max(d1['x']), np.min(d2['x']), np.max(d2['x']))
 
 
 