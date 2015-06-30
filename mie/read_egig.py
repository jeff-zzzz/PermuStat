# -*- coding: utf-8 -*-
"""
Created on Thu Mar 26 13:42:24 2015
@author: jesong1126
"""
# parsalate egig

import numpy as np
import glob 
from xml.dom.minidom import parse
import xml.etree.ElementTree as ET

filePath = '/Users/jesong1126/Python27/GeoPy/data/2941_HM'
ff = glob.glob(filePath+'/geometry.egig')[0]
tree = ET.parse(ff)     
root = tree.getroot()

for child in root: # THIS GETS ME ALL THE PATIENT ATTRIBUTES
    print child.tag 
    
    
    
a =  infoObj.getElementsByTagName('DIPOLE')
tree = ET.parse('Untitled.xml')
root = tree.getroot()
for child in root:
    print child.tag 
    child.find( "visits" )
    for x in child.iter("visit"):
        print x.tag, x.text
        

len(a)   

    infofile = filePath+'/info.xml'
    infotime = parse(infofile)    
    beginTime = infotime.getElementsByTagName('recordTime')[0].firstChild.data.encode()
    bTime = ns2pyTime(beginTime)

   
binfiles = []  
for ff in glob.glob(filePath+'/geometry.egig'):
    binfiles.append(ff[len(filePath)+1:])

signalFile = []
infoFile = []
if len(binfiles) > 1:
    for p in range(1, len(binfiles)): # p is from 1 to * to find PIBfile
        binFilename = binfiles[p]
        # All this to strip the number (binNumStr) from the signal file in order to apply
        # it to the info file. 
        prefix = 'signal' 
        prefixLen = len(prefix)
        extension = '.bin'
        extensionLen = len(extension)
        binNumDigits = len(binFilename) - (prefixLen + extensionLen)
        binNumStr = binFilename[np.asarray(range(prefixLen,  prefixLen + binNumDigits))]

        infoObjFile = filePath+'/info'+binNumStr+'.xml'
        infoObj = parse(infoObjFile) 
        # 'PNSData' 'Spectral' 'sourceData' 'JTF' 'TValues' 'Filter'
        if infoObj.getElementsByTagName('PNSData') != None :
            infoFile = 'info'+binNumStr+'.xml'
            signalFile = 'signal'+binNumStr+'.bin'

pibInfo = {'pibSignalFile':signalFile, 'pibInfoFile': infoFile}







