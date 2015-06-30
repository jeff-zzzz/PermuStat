# -*- coding: utf-8 -*-
"""
Created on Tue May 19 11:04:41 2015

@author: jesong1126
"""

import os 
from braink import read_lfm                 
import glob 

        
HMdir = '/Users/jesong1126/Python27/GeoPy/data/108_HM'
nE = 256
 
lfm_fname = glob.glob(HMdir+'/*.lfm')
if os.path.exists(HMdir+'/fdmForwardMatrixOriented'):          
    K  = read_lfm.forward(HMdir+'/fdmForwardMatrixOriented', nE)  
    print("fdmForwardMatrixOriented is loaded. ")            
elif len(lfm_fname) > 0:  # elif os.path.exists(HMdir+'/Leadfield.lfm'): 
    K  = read_lfm.lfm(lfm_fname[0], nE)                          
    print("Leadfield.lfm is loaded. ")      
else: 
    print("Either fdmForwardMatrixOriented or *.lfm are not loaded. ")            

nV = K.shape[1]

os.system("open /Applications/EAV/EAV.app")
#os.system("open /Users/jesong1126/Desktop/EAV.app")
from RabbitMQ import Connection  
SampleTime = 0 

c = Connection.Connection('localhost')

c.connect()
c.sendEEGData(K[:,SampleTime], "StatEEG", "")    
c.disconnect()

