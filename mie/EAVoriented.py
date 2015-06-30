# -*- coding: utf-8 -*-
"""
Created on Wed May  6 11:12:28 2015

@author: jesong1126
"""



import numpy as np 
import os
os.system("open /Applications/EAV/EAV.app")
from RabbitMQ import Connection 

from braink import read_lfm
#from swcm import Imatrix      
nC = 257 
nE = nC -1  
K  = read_lfm.lfm('/Users/jesong1126/Python27/data_GeoPy/108_HM/Leadfield.lfm', nE)                                  
nV = K.shape[1]
nV = 2383 
Dip1 = np.zeros(nV)
Dip1[0] = 1 

c = Connection.Connection('localhost')
c.connect()
#c.sendEEGData(Tmaps[:, i], "Tmap", "")
c.sendOrientedData(Dip1, "Oriented", "4:00pm")
#c.sendTriplesData(bbb, "Triples", "4:00pm")
#c.sendOrientedData(abs(sdenTmaps[:, i]), "sdenTmap", "")
#c.sendTriplesData(sdenTmaps[:,i] , "Triples", "4:00pm")
c.disconnect()       
 




