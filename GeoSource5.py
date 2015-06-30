# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'GeoSource.ui'
#
# Created: Mon Feb 16 12:30:09 2015
#      by: PyQt5 UI code generator 5.3.1
#
# WARNING! All changes made in this file will be lost!

#import sys
#from PyQt5 import QtCore, QtWidgets 
#from PyQt5.QtCore import QDir, Qt
#from PyQt5.QtGui import QFont, QPalette
#from PyQt5.QtWidgets import (QApplication, QCheckBox, QColorDialog, QDialog,
#        QErrorMessage, QFileDialog, QFontDialog, QFrame, QGridLayout,
#        QInputDialog, QLabel, QLineEdit, QMessageBox, QPushButton)

#from PyQt5.QtCore import pyqtSlot #, QtGui, 
#from PyQt5.QtWidgets import QApplication, QWidget #, QtGui,                
#import config 

#import matplotlib
#matplotlib.use("Qt5Agg")

import sip
sip.setapi('QString', 2)

import sys
from PyQt4 import QtCore, QtGui

try:
    _fromUtf8 = QtCore.QString.fromUtf8
except AttributeError:
    def _fromUtf8(s):
        return s

try:
    _encoding = QtGui.QApplication.UnicodeUTF8
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig, _encoding)
except AttributeError:
    def _translate(context, text, disambig):
        return QtGui.QApplication.translate(context, text, disambig)

    
class Ui_FormGS(object):
    def setupUi(self, FormGS):
        FormGS.setObjectName("FormGS")
        FormGS.resize(790, 373)
        self.pbEEG = QtWidgets.QPushButton(FormGS)
        self.pbEEG.setGeometry(QtCore.QRect(20, 20, 180, 60))
        self.pbEEG.setObjectName("pbEEG")        
        self.pbHM = QtWidgets.QPushButton(FormGS)
        self.pbHM.setGeometry(QtCore.QRect(210, 20, 180, 60))
        self.pbHM.setObjectName("pbHM")
        self.pbEst = QtWidgets.QPushButton(FormGS)
        self.pbEst.setGeometry(QtCore.QRect(400, 20, 180, 60))
        self.pbEst.setObjectName("pbEst")
        self.pbStat = QtWidgets.QPushButton(FormGS)
        self.pbStat.setGeometry(QtCore.QRect(590, 20, 180, 60))
        self.pbStat.setObjectName("pbStat")
        self.pbVis = QtWidgets.QPushButton(FormGS)
        self.pbVis.setGeometry(QtCore.QRect(590, 88, 180, 150))
        self.pbVis.setObjectName("pbVis")
        self.pbQuit = QtWidgets.QPushButton(FormGS)
        self.pbQuit.setGeometry(QtCore.QRect(590, 245, 180, 110))
        self.pbQuit.setObjectName("pbQuit")
        self.pbTopo = QtWidgets.QPushButton(FormGS)
        self.pbTopo.setGeometry(QtCore.QRect(400, 245, 180, 110))
        self.pbTopo.setObjectName("pbTopo")
        self.frame = QtWidgets.QFrame(FormGS)
        self.frame.setGeometry(QtCore.QRect(25, 90, 170, 140))
        self.frame.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.frame.setFrameShadow(QtWidgets.QFrame.Raised)
        self.frame.setObjectName("frame")
        self.cbTrial = QtWidgets.QComboBox(self.frame)
        self.cbTrial.setGeometry(QtCore.QRect(50, 30, 110, 26))
        self.cbTrial.setObjectName("cbTrial")
        self.labelTrial = QtWidgets.QLabel(self.frame)
        self.labelTrial.setGeometry(QtCore.QRect(10, 35, 40, 16))
        self.labelTrial.setObjectName("labelTrial")
        self.labelThreshold = QtWidgets.QLabel(self.frame)
        self.labelThreshold.setGeometry(QtCore.QRect(10, 80, 80, 16))
        self.labelThreshold.setObjectName("labelThreshold")
        self.sliderThreshold = QtWidgets.QSlider(self.frame)
        self.sliderThreshold.setGeometry(QtCore.QRect(15, 105, 141, 22))
        self.sliderThreshold.setMaximum(100)
        self.sliderThreshold.setOrientation(QtCore.Qt.Horizontal)
        self.sliderThreshold.setObjectName("sliderThreshold")
        self.listThreshold = QtWidgets.QSpinBox(self.frame)
        self.listThreshold.setGeometry(QtCore.QRect(100, 77, 57, 26))
        self.listThreshold.setObjectName("listThreshold")
        self.txtEpoch = QtWidgets.QFrame(FormGS)
        self.txtEpoch.setGeometry(QtCore.QRect(215, 90, 170, 140))
        self.txtEpoch.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.txtEpoch.setFrameShadow(QtWidgets.QFrame.Raised)
        self.txtEpoch.setObjectName("txtEpoch")
        self.labelStart = QtWidgets.QLabel(self.txtEpoch)
        self.labelStart.setGeometry(QtCore.QRect(20, 35, 40, 16))
        self.labelStart.setObjectName("labelStart")
        self.labelEnd = QtWidgets.QLabel(self.txtEpoch)
        self.labelEnd.setGeometry(QtCore.QRect(20, 85, 40, 16))
        self.labelEnd.setObjectName("labelEnd")
        self.editStart = QtWidgets.QLineEdit(self.txtEpoch)
        self.editStart.setGeometry(QtCore.QRect(80, 30, 70, 23))
        self.editStart.setObjectName("editStart")
        self.editEnd = QtWidgets.QLineEdit(self.txtEpoch)
        self.editEnd.setGeometry(QtCore.QRect(80, 80, 70, 23))
        self.editEnd.setObjectName("editEnd")
#        self.editEnd.setReadOnly(True) 
        self.txtMovie = QtWidgets.QFrame(FormGS)
        self.txtMovie.setGeometry(QtCore.QRect(405, 90, 170, 140))
        self.txtMovie.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.txtMovie.setFrameShadow(QtWidgets.QFrame.Raised)
        self.txtMovie.setObjectName("txtMovie")
        self.rbMovie = QtWidgets.QRadioButton(self.txtMovie)
        self.rbMovie.setGeometry(QtCore.QRect(20, 20, 102, 20))
        self.rbMovie.setChecked(True)
        self.rbMovie.setObjectName("rbMovie")
        self.rbAveragEpoch = QtWidgets.QRadioButton(self.txtMovie)
        self.rbAveragEpoch.setGeometry(QtCore.QRect(20, 90, 120, 20))
        self.rbAveragEpoch.setObjectName("rbAveragEpoch")
        self.cbColor = QtWidgets.QCheckBox(self.txtMovie)
        self.cbColor.setGeometry(QtCore.QRect(40, 50, 90, 20))
        self.cbColor.setObjectName("cbColor")
        self.txtTopomap = QtWidgets.QFrame(FormGS)
        self.txtTopomap.setGeometry(QtCore.QRect(25, 250, 360, 100))
        self.txtTopomap.setFrameShape(QtWidgets.QFrame.StyledPanel)
        self.txtTopomap.setFrameShadow(QtWidgets.QFrame.Raised)
        self.txtTopomap.setObjectName("txtTopomap")
        self.rbEEG = QtWidgets.QRadioButton(self.txtTopomap)
        self.rbEEG.setGeometry(QtCore.QRect(10, 40, 50, 20))
        self.rbEEG.setChecked(True)
        self.rbEEG.setObjectName("rbEEG")
        self.rbAllForward = QtWidgets.QRadioButton(self.txtTopomap)
        self.rbAllForward.setGeometry(QtCore.QRect(100, 20, 100, 20))
        self.rbAllForward.setObjectName("rbAllForward")
        self.rbAllResidual = QtWidgets.QRadioButton(self.txtTopomap)
        self.rbAllResidual.setGeometry(QtCore.QRect(100, 60, 100, 20))
        self.rbAllResidual.setObjectName("rbAllResidual")
        self.rbPartForward = QtWidgets.QRadioButton(self.txtTopomap)
        self.rbPartForward.setGeometry(QtCore.QRect(225, 20, 102, 20))
        self.rbPartForward.setObjectName("rbPartForward")
        self.rbPartResidual = QtWidgets.QRadioButton(self.txtTopomap)
        self.rbPartResidual.setGeometry(QtCore.QRect(225, 60, 102, 20))
        self.rbPartResidual.setObjectName("rbPartResidual")

        self.retranslateUi(FormGS)
        self.pbEEG.clicked.connect(self.pbEEG_Callback)
        self.pbHM.clicked.connect(self.pbHM_Callback) 
        self.pbEst.clicked.connect(self.pbEst_Callback)
        self.pbStat.clicked.connect(self.pbStat_Callback)
        self.sliderThreshold.valueChanged['int'].connect(self.listThreshold.setValue)
        self.listThreshold.valueChanged['int'].connect(self.sliderThreshold.setValue) 
        self.listThreshold.editingFinished.connect(self.editTH_Callback) 
        
        self.editStart.editingFinished.connect(self.editStart_Callback)
        self.editEnd.editingFinished.connect(self.editEnd_Callback)
        
        
        self.pbVis.clicked.connect(self.pbVis_Callback)
        self.pbQuit.clicked.connect(FormGS.close)
#        import HeadModelOpen 
#        self.pbQuit.clicked.connect(self.FormHM.close)
#        self.pbQuit.clicked.connect(self.FormSourceRun.close)
                
        QtCore.QMetaObject.connectSlotsByName(FormGS)
        FormGS.setTabOrder(self.pbEEG, self.pbHM)
        FormGS.setTabOrder(self.pbHM, self.pbEst)
        FormGS.setTabOrder(self.pbEst, self.pbStat)
        FormGS.setTabOrder(self.pbStat, self.cbTrial)
        FormGS.setTabOrder(self.cbTrial, self.listThreshold)
        FormGS.setTabOrder(self.listThreshold, self.sliderThreshold)
        FormGS.setTabOrder(self.sliderThreshold, self.editStart)
        FormGS.setTabOrder(self.editStart, self.editEnd)
        FormGS.setTabOrder(self.editEnd, self.rbMovie)
        FormGS.setTabOrder(self.rbMovie, self.cbColor)
        FormGS.setTabOrder(self.cbColor, self.rbAveragEpoch)
        FormGS.setTabOrder(self.rbAveragEpoch, self.pbVis)
        FormGS.setTabOrder(self.pbVis, self.rbEEG)
        FormGS.setTabOrder(self.rbEEG, self.rbAllForward)
        FormGS.setTabOrder(self.rbAllForward, self.rbAllResidual)
        FormGS.setTabOrder(self.rbAllResidual, self.rbPartForward)
        FormGS.setTabOrder(self.rbPartForward, self.rbPartResidual)
        FormGS.setTabOrder(self.rbPartResidual, self.pbTopo)
        FormGS.setTabOrder(self.pbTopo, self.pbQuit)

    def retranslateUi(self, FormGS):
        _translate = QtCore.QCoreApplication.translate
        FormGS.setWindowTitle(_translate("FormGS", "GeoSource"))
        self.pbEEG.setText(_translate("FormGS", "EEG"))
        self.pbHM.setText(_translate("FormGS", "Head Model"))
        self.pbEst.setText(_translate("FormGS", "Estimate"))
        self.pbStat.setText(_translate("FormGS", "Stat"))
        self.pbVis.setText(_translate("FormGS", "Visualize"))
        self.pbQuit.setText(_translate("FormGS", "Quit"))
        self.pbTopo.setText(_translate("FormGS", "Topomap"))
        self.labelTrial.setText(_translate("FormGS", "Trial"))
        self.labelThreshold.setText(_translate("FormGS", "Threshold(%)"))
#        self.listThreshold.setText(_translate("FormGS", "0"))
        self.labelStart.setText(_translate("FormGS", "Start"))
        self.labelEnd.setText(_translate("FormGS", "End"))
        self.editStart.setText(_translate("FormGS", "0"))
        self.editEnd.setText(_translate("FormGS", "0"))
        self.rbMovie.setText(_translate("FormGS", "Movie"))
        self.rbAveragEpoch.setText(_translate("FormGS", "Average Epoch"))
        self.cbColor.setText(_translate("FormGS", "Color Lock"))
        self.rbEEG.setText(_translate("FormGS", "EEG"))
        self.rbAllForward.setText(_translate("FormGS", "All Forward"))
        self.rbAllResidual.setText(_translate("FormGS", "All Residual"))
        self.rbPartForward.setText(_translate("FormGS", "Part Forward"))
        self.rbPartResidual.setText(_translate("FormGS", "Part Residual"))


    def pbEEG_Callback(self):

        global s, nC, nSamples, nSamplesPre, nTrials, srate, trialsName 
        global categoryName, nCategory  

        import config 
        from mff import read_mff_header, read_mff_data 
        import matplotlib.pyplot as plt
        import numpy as np
        
        options = QtWidgets.QFileDialog.Options()
        mfffilename, _ = QtWidgets.QFileDialog.getOpenFileNames(self.pbEEG, "Select a mff file", "", "mff files (*.mff)", options= options) 
        filePath = mfffilename[0] 
        hdr = read_mff_header.read_mff_header(filePath)       
        
        nC = hdr['nChans'] ; config.nC = nC
        nSamples = hdr['nSamples'] ; config.nSamples = nSamples
        nSamplesPre = hdr['nSamplesPre'] ; config.nSamplesPre = nSamplesPre
        nTrials = hdr['nTrials'] ; config.nTrials = nTrials
        srate = hdr['Fs'] ; config.srate = srate
        summaryInfo = hdr['orig'] ; # config.nC = nC
        trialsName = summaryInfo['epochLabels'] ; config.trialsName = trialsName
        categoryName = list(set(trialsName)) ; config.categoryName = categoryName
        nCategory = len(categoryName) ; config.nCategory = nCategory

        for i in range(nTrials):
            self.cbTrial.addItem(trialsName[i])

        data = read_mff_data.read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr)    

        SegStatus = hdr['orig']['epochSegStatus']  
        GoodSeg = []; BadSeg = []; #epochLabels = []; 
        for i in range(len(SegStatus)):
            if SegStatus[i] == 'bad' :  
                BadSeg.append(i)
            else :  
                GoodSeg.append(i)
        
        nTrials = len(GoodSeg)

        print('nC:%d'  %nC)
        print('nSamples:%d'  %nSamples)
        print('nTrials:%d'  %nTrials)
        print('nCategory:%d' %nCategory)
        print('categoryName:%s' %categoryName)
#        print('categoryName:%s' %trialsName)

        ## average reference 
        H = np.identity(nC) - np.ones((nC, nC))/nC 
        s = np.zeros((nC, nSamples, nTrials))
        sGFP = np.zeros((nSamples, nTrials))
        if nTrials > 1:
            for i in range(nTrials):
                s[:,:,i] = np.dot(H, data[:,:,GoodSeg[i]])
                sGFP[:,i] = np.std(s[:,:,i], axis=0)                  
#                s2 = s[:,:,i] * s[:,:,i] 
#                gfp = s2.mean(axis=0)      
#                sGFP = np.sqrt(gfp)
        else :
            s = np.dot(H, data)  
            sGFP = np.std(s, axis=0)                  
#            s2 = s * s 
#            gfp = s2.mean(axis=0)      
#            sGFP = np.sqrt(gfp)

        baseline = (nSamplesPre * 1000 / srate) ; config.baseline = baseline
        msSamples = np.arange(0, nSamples,1) * 1000/srate  - baseline 
        config.msSamples = msSamples

#        global f1
        if nTrials == nCategory :
#            f1 = plt.figure 
            plt.plot(msSamples, sGFP) 
            plt.title('gfp') 
            plt.xlabel('ms') 
            plt.legend(trialsName, loc='best')
            plt.show() 
            print('EEG is ERP.')  
            #eventsName = trialsName 
        else : 
#            f1 = plt.figure  
            plt.plot(msSamples, sGFP) 
            plt.title('gfp') 
            plt.xlabel('ms') 
            plt.show() 
            print('EEG is Single trial. There are',nTrials,'trials. ') 
            #eventsName = categoryName
        
        config.nSamples = nSamples
        config.msSamples = msSamples
        config.nTrials = nTrials
        config.trialsName = trialsName
        config.s = s 
        config.sGFP = sGFP
        print("EEG data is loaded.")
        # self.pbEEG.configure(background = "red")          
        return s  

    def pbHM_Callback(self):
        
        options = QtWidgets.QFileDialog.DontResolveSymlinks | QtWidgets.QFileDialog.ShowDirsOnly
        HMdir = QtWidgets.QFileDialog.getExistingDirectory(self.pbHM,
                "Select a Head Model Directory.",
                self.pbHM.text(), options=options)    
                
        import os 
        from braink import read_lfm                 
        import config                
        
        nE = config.nE
        nC = nE + 1
        if os.path.exists(HMdir+'/fdmForwardMatrixOriented'):          
            K  = read_lfm.forward(HMdir+'/fdmForwardMatrixOriented', config.nE)  
            print K 
            print("fdmForwardMatrixOriented is loaded. ")            
        elif os.path.exists(HMdir+'/Leadfield.lfm'): 
            K  = read_lfm.lfm(HMdir+'/Leadfield.lfm', config.nE)                          
            print("Leadfield.lfm is loaded. ")      
            print K 
        else:
            print("LFM is not loaded. ")            
        
        nV = K.shape[1]
        config.K = K    
        config.nV = nV    
        config.nC = nC    
        
        os.system("open /Applications/EAV/EAV.app")

        return  

    def pbEst_Callback(self): 
        from SourceRun import Ui_FormSourceRun
        self.FormSourceRun = QtWidgets.QWidget()
        self.uiEst = Ui_FormSourceRun()
        self.uiEst.setupUi(self.FormSourceRun)
        self.FormSourceRun.show()        
        return #sden    
    
    def editTH_Callback(self):
        global thresholdcutoff ;
        thresholdcutoff = self.listThreshold.value()
        print('Lowest %d(percent) of dipoles are cut off.' %thresholdcutoff)
        return 
              
    def editStart_Callback(self):
        import config 
        startMS = self.editStart.text()
        startMS = int(startMS)
        config.startMS = startMS 
        return         

    def editEnd_Callback(self):
        import config 
        import matplotlib.pyplot as plt
        import numpy as np        
        endMS = self.editEnd.text()
        endMS = int(endMS)        
        config.endMS = endMS 
        print('From %d to %d ms' %(config.startMS, config.endMS)) 
#         f1 = plt.figure  
        plt.plot(config.msSamples, config.sGFP) 
        plt.title('sum(phi^2)') 
        plt.xlabel('ms')             
        plt.plot([config.startMS, config.startMS], [np.min(config.sGFP), np.max(config.sGFP)], 'k-', lw=2)
        plt.plot([config.endMS, config.endMS], [np.min(config.sGFP), np.max(config.sGFP)], 'k-', lw=2)
        if config.nTrials == config.nCategory :
            plt.legend(config.trialsName, loc='best')
#        else:            
        plt.show() 
        print('EEG is Single trial. There are %d trials.' %config.nTrials) 
        return         
    

    def pbStat_Callback(self):                      
        from StatRun import Ui_FormStat
        self.FormStat = QtWidgets.QWidget()
        self.uiStat = Ui_FormStat()
        self.uiStat.setupUi(self.FormStat)
        self.FormStat.show()        
        return #sden    


    def pbVis_Callback(self): 
        import matplotlib.pyplot as plt 
        import numpy as np 
        import config 
#        import os
#        os.system("open /Applications/EAV/EAV.app")

        sdenMean = config.sdenMean
        sdenSigma = config.sdenSigma
        sdenT = config.sdenT
        sMean = config.sMean
        sSigma = config.sSigma
        sT = config.sT
        nSigChan = config.nSigChan
        nSigVox = config.nSigVox

        msSamples = config.msSamples
        startMS = config.startMS 
        ii = np.where( msSamples == startMS )[0]

        fig, ax = plt.subplots(4,2,figsize=(20,10))
        ax[0,0].plot(msSamples, sMean.T)
        ax[0,0].set_xlim([msSamples[0], msSamples[-1]])
        ax[0,0].set_ylabel('mean')
        ax[0,0].set_title('EEG')
        ax[1,0].plot(msSamples, sSigma.T)
        ax[1,0].set_xlim([msSamples[0], msSamples[-1]])
        ax[1,0].set_ylabel('sdv')
        ax[2,0].plot(msSamples, sT.T)
        ax[2,0].set_xlim([msSamples[0], msSamples[-1]])
        ax[2,0].set_ylabel('Normal')
#        ax[2,0].axhline(vu, color='black', lw=2)
#        ax[2,0].axhline(vl, color='black', lw=2)
        ax[3,0].plot(msSamples, nSigChan )
        ax[3,0].axvline(msSamples[ii],  color='black', lw=2)
        #ax[3,0].set_ylim([-1, max(nSigChan)+1])
        ax[3,0].set_xlim([msSamples[0], msSamples[-1]])
        ax[3,0].set_ylabel('Num Sig.Channels')
        ax[3,0].set_xlabel('ms')
        #ax[3,0].legend()
        ax[0,1].plot(msSamples, sdenMean.T)
        ax[0,1].set_xlim([msSamples[0], msSamples[-1]])
        ax[0,1].set_ylabel('mean')
        ax[0,1].set_title('Source')
        ax[1,1].plot(msSamples, sdenSigma.T)
        ax[1,1].set_xlim([msSamples[0], msSamples[-1]])
        ax[1,1].set_ylabel('sdv')
        ax[2,1].plot(msSamples, sdenT.T)
        ax[2,1].set_xlim([msSamples[0], msSamples[-1]])
        ax[2,1].set_ylabel('Normal')
#        ax[2,1].axhline(sdenvu, color='black', lw=2)
#        ax[2,1].axhline(sdenvl, color='black', lw=2)
        ax[3,1].plot(msSamples, nSigVox,  c='black', label= 'l+r')
        ax[3,1].axvline(msSamples[ii],  color='black', lw=2)
        #ax[3,1].set_ylim([-5, max(nSigVox)+5])
        ax[3,1].set_xlim([msSamples[0], msSamples[-1]])
        ax[3,1].set_ylabel('Num Sig.Voxels')
        ax[3,1].set_xlabel('ms')
        #ax[3,1].legend()
        plt.show()
        
        
        sdenTSig = config.sdenTSig 
        sT = config.sT
#        sTSig = config.sTSig 
        
        from RabbitMQ import Connection
        c = Connection.Connection('localhost')
        c.connect()
        c.sendOrientedData(abs(sdenTSig[:, ii]) , "Oriented", "3:00pm")
        c.sendEEGData(sT[:,ii], "Oriented", "3:00pm")
        c.disconnect()       

#        print os.getcwd()
#        os.chdir("/Users/jesong1126/Python27/GeoPy") 
                   
    
if __name__ == "__main__":
    app = QtWidgets.QApplication(sys.argv)
    FormGS = QtWidgets.QWidget()
    ui = Ui_FormGS()
    ui.setupUi(FormGS)
    FormGS.show()
    sys.exit(app.exec_())


