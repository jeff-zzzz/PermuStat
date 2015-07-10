# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'Microstate.ui'
#
# Created: Thu Jul  2 09:33:47 2015
#      by: PyQt4 UI code generator 4.11.1
#
# WARNING! All changes made in this file will be lost!

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

class Ui_FormMicro(object):
    def setupUi(self, FormMicro):
        FormMicro.setObjectName(_fromUtf8("FormMicro"))
        FormMicro.resize(431, 618)
        self.pbClose = QtGui.QPushButton(FormMicro)
        self.pbClose.setGeometry(QtCore.QRect(10, 550, 401, 50))
        self.pbClose.setObjectName(_fromUtf8("pbClose"))
        self.toolBoxMicro = QtGui.QToolBox(FormMicro)
        self.toolBoxMicro.setGeometry(QtCore.QRect(20, 20, 391, 511))
        self.toolBoxMicro.setObjectName(_fromUtf8("toolBoxMicro"))
        self.pageData = QtGui.QWidget()
        self.pageData.setGeometry(QtCore.QRect(0, 0, 391, 375))
        self.pageData.setObjectName(_fromUtf8("pageData"))
        self.frameScalp = QtGui.QFrame(self.pageData)
        self.frameScalp.setGeometry(QtCore.QRect(30, 200, 320, 71))
        self.frameScalp.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frameScalp.setFrameShadow(QtGui.QFrame.Raised)
        self.frameScalp.setObjectName(_fromUtf8("frameScalp"))
        self.rbScalp = QtGui.QRadioButton(self.frameScalp)
        self.rbScalp.setGeometry(QtCore.QRect(20, 10, 70, 20))
        self.rbScalp.setChecked(True)
        self.rbScalp.setObjectName(_fromUtf8("rbScalp"))
        self.rbSource = QtGui.QRadioButton(self.frameScalp)
        self.rbSource.setGeometry(QtCore.QRect(20, 40, 70, 20))
        self.rbSource.setObjectName(_fromUtf8("rbSource"))
        self.frameLevel = QtGui.QFrame(self.pageData)
        self.frameLevel.setGeometry(QtCore.QRect(30, 290, 320, 71))
        self.frameLevel.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frameLevel.setFrameShadow(QtGui.QFrame.Raised)
        self.frameLevel.setObjectName(_fromUtf8("frameLevel"))
        self.rbIndividual = QtGui.QRadioButton(self.frameLevel)
        self.rbIndividual.setGeometry(QtCore.QRect(20, 10, 160, 20))
        self.rbIndividual.setChecked(True)
        self.rbIndividual.setObjectName(_fromUtf8("rbIndividual"))
        self.rbGroup = QtGui.QRadioButton(self.frameLevel)
        self.rbGroup.setGeometry(QtCore.QRect(20, 40, 160, 20))
        self.rbGroup.setObjectName(_fromUtf8("rbGroup"))
        self.frameData = QtGui.QFrame(self.pageData)
        self.frameData.setGeometry(QtCore.QRect(30, 0, 320, 181))
        self.frameData.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frameData.setFrameShadow(QtGui.QFrame.Raised)
        self.frameData.setObjectName(_fromUtf8("frameData"))
        self.pbMFF = QtGui.QPushButton(self.frameData)
        self.pbMFF.setGeometry(QtCore.QRect(10, 10, 301, 40))
        self.pbMFF.setObjectName(_fromUtf8("pbMFF"))
        self.pbHeadModel = QtGui.QPushButton(self.frameData)
        self.pbHeadModel.setGeometry(QtCore.QRect(10, 50, 301, 40))
        self.pbHeadModel.setObjectName(_fromUtf8("pbHeadModel"))
        self.labelMethod = QtGui.QLabel(self.frameData)
        self.labelMethod.setGeometry(QtCore.QRect(30, 100, 71, 16))
        self.labelMethod.setObjectName(_fromUtf8("labelMethod"))
        self.cbMethod = QtGui.QComboBox(self.frameData)
        self.cbMethod.setGeometry(QtCore.QRect(160, 95, 141, 26))
        self.cbMethod.setObjectName(_fromUtf8("cbMethod"))
        self.cbMethod.addItem(_fromUtf8(""))
        self.cbMethod.addItem(_fromUtf8(""))
        self.cbMethod.addItem(_fromUtf8(""))
        self.cbMethod.addItem(_fromUtf8(""))
        self.pbEstimate = QtGui.QPushButton(self.frameData)
        self.pbEstimate.setGeometry(QtCore.QRect(10, 130, 301, 40))
        self.pbEstimate.setObjectName(_fromUtf8("pbEstimate"))
        
        self.toolBoxMicro.addItem(self.pageData, _fromUtf8(""))
        self.pageMicro = QtGui.QWidget()
        self.pageMicro.setObjectName(_fromUtf8("pageMicro"))
        self.labelCategory = QtGui.QLabel(self.pageMicro)
        self.labelCategory.setGeometry(QtCore.QRect(60, 10, 81, 16))
        self.labelCategory.setObjectName(_fromUtf8("labelCategory"))
        self.cbCategory = QtGui.QComboBox(self.pageMicro)
        self.cbCategory.setGeometry(QtCore.QRect(215, 10, 110, 26))
        self.cbCategory.setEditable(False)
        self.cbCategory.setObjectName(_fromUtf8("cbCategory"))
        self.labelNumCluster = QtGui.QLabel(self.pageMicro)
        self.labelNumCluster.setGeometry(QtCore.QRect(60, 45, 110, 16))
        self.labelNumCluster.setObjectName(_fromUtf8("labelNumCluster"))
        self.editNumCluster = QtGui.QLineEdit(self.pageMicro)
        self.editNumCluster.setGeometry(QtCore.QRect(220, 45, 100, 21))
        self.editNumCluster.setCursorPosition(1)
        self.editNumCluster.setObjectName(_fromUtf8("editNumCluster"))
        
        self.labelSamplesMS = QtGui.QLabel(self.pageMicro)
        self.labelSamplesMS.setGeometry(QtCore.QRect(60, 80, 130, 16))
        self.labelSamplesMS.setObjectName(_fromUtf8("labelSamplesMS"))
        self.editSamplesMS = QtGui.QLineEdit(self.pageMicro)
        self.editSamplesMS.setGeometry(QtCore.QRect(220, 80, 100, 21))
        self.editSamplesMS.setCursorPosition(1)
        self.editSamplesMS.setObjectName(_fromUtf8("editSamplesMS"))
        
        self.labelSamplesMSto = QtGui.QLabel(self.pageMicro)
        self.labelSamplesMSto.setGeometry(QtCore.QRect(60, 115, 130, 16))
        self.labelSamplesMSto.setObjectName(_fromUtf8("labelSamplesMSto"))
        self.editSamplesMSto = QtGui.QLineEdit(self.pageMicro)
        self.editSamplesMSto.setGeometry(QtCore.QRect(220, 115, 100, 21))
        self.editSamplesMSto.setCursorPosition(1)
        self.editSamplesMSto.setObjectName(_fromUtf8("editSamplesMSto"))


        self.labelThreshold = QtGui.QLabel(self.pageMicro)
        self.labelThreshold.setGeometry(QtCore.QRect(60, 150, 110, 16))
        self.labelThreshold.setObjectName(_fromUtf8("labelThreshold"))
        self.editThreshold = QtGui.QLineEdit(self.pageMicro)
        self.editThreshold.setGeometry(QtCore.QRect(220, 150, 100, 21))
        self.editThreshold.setCursorPosition(1)
        self.editThreshold.setObjectName(_fromUtf8("editThreshold"))
        self.labelNumSim = QtGui.QLabel(self.pageMicro)
        self.labelNumSim.setGeometry(QtCore.QRect(60, 185, 130, 16))
        self.labelNumSim.setObjectName(_fromUtf8("labelNumSim"))
        self.editNumSim = QtGui.QLineEdit(self.pageMicro)
        self.editNumSim.setGeometry(QtCore.QRect(220, 185, 100, 21))
        self.editNumSim.setCursorPosition(1)
        self.editNumSim.setObjectName(_fromUtf8("editNumSim"))
        self.pbRunMicro = QtGui.QPushButton(self.pageMicro)
        self.pbRunMicro.setGeometry(QtCore.QRect(30, 220, 311, 40))
        self.pbRunMicro.setObjectName(_fromUtf8("pbRunMicro"))

        self.labelMicro = QtGui.QLabel(self.pageMicro)
        self.labelMicro.setGeometry(QtCore.QRect(60, 285, 130, 16))
        self.labelMicro.setObjectName(_fromUtf8("labelMicro"))
        self.editMicro = QtGui.QLineEdit(self.pageMicro)
        self.editMicro.setGeometry(QtCore.QRect(220, 285, 100, 21))
        self.editMicro.setCursorPosition(1)
        self.editMicro.setObjectName(_fromUtf8("editMicro"))
        
        
        self.pbViewMicro = QtGui.QPushButton(self.pageMicro)
        self.pbViewMicro.setGeometry(QtCore.QRect(30, 320, 311, 40))
        self.pbViewMicro.setObjectName(_fromUtf8("pbViewMicro"))
        
        
        
        self.toolBoxMicro.addItem(self.pageMicro, _fromUtf8(""))
        self.pageNetwork = QtGui.QWidget()
        self.pageNetwork.setObjectName(_fromUtf8("pageNetwork"))
        self.labelFrequency = QtGui.QLabel(self.pageNetwork)
        self.labelFrequency.setGeometry(QtCore.QRect(10, 30, 91, 16))
        self.labelFrequency.setObjectName(_fromUtf8("labelFrequency"))
        self.editFreqLow = QtGui.QLineEdit(self.pageNetwork)
        self.editFreqLow.setGeometry(QtCore.QRect(140, 30, 61, 21))
        self.editFreqLow.setCursorPosition(4)
        self.editFreqLow.setObjectName(_fromUtf8("editFreqLow"))
        self.pbRunNetwork = QtGui.QPushButton(self.pageNetwork)
        self.pbRunNetwork.setGeometry(QtCore.QRect(0, 70, 320, 50))
        self.pbRunNetwork.setObjectName(_fromUtf8("pbRunNetwork"))
        self.pbViewNetwork = QtGui.QPushButton(self.pageNetwork)
        self.pbViewNetwork.setGeometry(QtCore.QRect(0, 130, 320, 50))
        self.pbViewNetwork.setObjectName(_fromUtf8("pbViewNetwork"))
        self.editFreqHigh = QtGui.QLineEdit(self.pageNetwork)
        self.editFreqHigh.setGeometry(QtCore.QRect(230, 30, 61, 21))
        self.editFreqHigh.setCursorPosition(4)
        self.editFreqHigh.setObjectName(_fromUtf8("editFreqHigh"))
        self.label_13 = QtGui.QLabel(self.pageNetwork)
        self.label_13.setGeometry(QtCore.QRect(210, 30, 16, 20))
        self.label_13.setObjectName(_fromUtf8("label_13"))
        self.toolBoxMicro.addItem(self.pageNetwork, _fromUtf8(""))
        self.pageOthers = QtGui.QWidget()
        self.pageOthers.setGeometry(QtCore.QRect(0, 0, 391, 375))
        self.pageOthers.setObjectName(_fromUtf8("pageOthers"))
        self.labelTimeSeries = QtGui.QLabel(self.pageOthers)
        self.labelTimeSeries.setGeometry(QtCore.QRect(20, 20, 91, 16))
        self.labelTimeSeries.setObjectName(_fromUtf8("labelTimeSeries"))
        self.labelAR = QtGui.QLabel(self.pageOthers)
        self.labelAR.setGeometry(QtCore.QRect(20, 50, 91, 16))
        self.labelAR.setObjectName(_fromUtf8("labelAR"))
        self.labelPCA = QtGui.QLabel(self.pageOthers)
        self.labelPCA.setGeometry(QtCore.QRect(20, 80, 91, 16))
        self.labelPCA.setObjectName(_fromUtf8("labelPCA"))
        self.labelICA = QtGui.QLabel(self.pageOthers)
        self.labelICA.setGeometry(QtCore.QRect(20, 110, 91, 16))
        self.labelICA.setObjectName(_fromUtf8("labelICA"))
        self.labelDCM = QtGui.QLabel(self.pageOthers)
        self.labelDCM.setGeometry(QtCore.QRect(20, 140, 91, 16))
        self.labelDCM.setObjectName(_fromUtf8("labelDCM"))
        self.labelGranger = QtGui.QLabel(self.pageOthers)
        self.labelGranger.setGeometry(QtCore.QRect(20, 170, 91, 16))
        self.labelGranger.setObjectName(_fromUtf8("labelGranger"))
        self.toolBoxMicro.addItem(self.pageOthers, _fromUtf8(""))

        ##---------------------------------------------------------------------
        self.pbMFF.clicked.connect(self.pbMFF_Callback)
        self.pbHeadModel.clicked.connect(self.pbHeadModel_Callback) 
        self.cbMethod.activated.connect(self.cbMethod_Callback)
        self.pbEstimate.clicked.connect(self.pbEstimate_Callback) 

        self.cbCategory.activated.connect(self.cbCategory_Callback)
        self.pbRunMicro.clicked.connect(self.pbRunMicro_Callback) 
        self.pbViewMicro.clicked.connect(self.pbViewMicro_Callback) 
        
        
        self.pbClose.clicked.connect(self.pbClose_Callback)
        self.pbClose.clicked.connect(FormMicro.close)
        
        
        self.retranslateUi(FormMicro)
        self.toolBoxMicro.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(FormMicro)

    ##-------------------------------------------------------------------------
    def retranslateUi(self, FormMicro):
        FormMicro.setWindowTitle(_translate("FormMicro", "Microstates Analysis", None))
        self.pbClose.setText(_translate("FormMicro", "Close", None))
        self.rbScalp.setText(_translate("FormMicro", "Scalp", None))
        self.rbSource.setText(_translate("FormMicro", "Source", None))
        self.rbIndividual.setText(_translate("FormMicro", "Individual (First Level)", None))
        self.rbGroup.setText(_translate("FormMicro", "Group (Second level)", None))
        self.pbMFF.setText(_translate("FormMicro", "MFF", None))
        self.pbHeadModel.setText(_translate("FormMicro", "Head Model", None))
        self.labelMethod.setText(_translate("FormMicro", "Method : ", None))
        self.cbMethod.setItemText(0, _translate("FormMicro", "MN", None))
        self.cbMethod.setItemText(1, _translate("FormMicro", "sMN", None))
        self.cbMethod.setItemText(2, _translate("FormMicro", "CSL", None))
        self.cbMethod.setItemText(3, _translate("FormMicro", "LORETA", None))
        self.pbEstimate.setText(_translate("FormMicro", "Estimate", None))
        self.toolBoxMicro.setItemText(self.toolBoxMicro.indexOf(self.pageData), _translate("FormMicro", "Data", None))
        self.labelNumCluster.setText(_translate("FormMicro", "Num of Cluster", None))
        self.editNumCluster.setText(_translate("FormMicro", "5", None))
        self.labelSamplesMS.setText(_translate("FormMicro", "Time (ms) from", None))
        self.editSamplesMS.setText(_translate("FormMicro", "0", None))
        self.labelSamplesMSto.setText(_translate("FormMicro", "Time (ms) to", None))
        self.editSamplesMSto.setText(_translate("FormMicro", "0", None))

        self.labelThreshold.setText(_translate("FormMicro", "CI (%)", None))
        self.editThreshold.setText(_translate("FormMicro", "95", None))
        self.labelNumSim.setText(_translate("FormMicro", "Num of Simulation", None))
        self.editNumSim.setText(_translate("FormMicro", "1000", None))
        self.pbRunMicro.setText(_translate("FormMicro", "Run", None))         
        self.labelMicro.setText(_translate("FormMicro", "Microstate", None))
        self.editMicro.setText(_translate("FormMicro", "1", None))        
        self.pbViewMicro.setText(_translate("FormMicro", "View in Reciprocity", None))

        self.labelCategory.setText(_translate("FormMicro", "Category", None))
        self.toolBoxMicro.setItemText(self.toolBoxMicro.indexOf(self.pageMicro), _translate("FormMicro", "Microstates", None))
        self.labelFrequency.setText(_translate("FormMicro", "Frequency", None))
        self.editFreqLow.setInputMask(_translate("FormMicro", "0.05", None))
        self.editFreqLow.setText(_translate("FormMicro", "0.05", None))
        self.pbRunNetwork.setText(_translate("FormMicro", "Run", None))
        self.pbViewNetwork.setText(_translate("FormMicro", "View in Reciprocity", None))
        self.editFreqHigh.setInputMask(_translate("FormMicro", "0.05", None))
        self.editFreqHigh.setText(_translate("FormMicro", "0.05", None))
        self.label_13.setText(_translate("FormMicro", "~", None))
        self.toolBoxMicro.setItemText(self.toolBoxMicro.indexOf(self.pageNetwork), _translate("FormMicro", "Network Analysis", None))
        self.labelTimeSeries.setText(_translate("FormMicro", "Time Series", None))
        self.labelAR.setText(_translate("FormMicro", "AR", None))
        self.labelPCA.setText(_translate("FormMicro", "PCA", None))
        self.labelICA.setText(_translate("FormMicro", "ICA", None))
        self.labelDCM.setText(_translate("FormMicro", "DCM", None))
        self.labelGranger.setText(_translate("FormMicro", "Granger", None))
                 
        self.toolBoxMicro.setItemText(self.toolBoxMicro.indexOf(self.pageOthers), _translate("FormMicro", "Others", None))

    ##-------------------------------------------------------------------------                
    def pbMFF_Callback(self):

        global s, nC, nSamples, nSamplesPre, nTrials, srate, trialsName 
        global categoryName, nCategory  

        import config 
        from mff import read_mff_header, read_mff_data 
        import numpy as np
 
        options = QtGui.QFileDialog.Options()
        mfffilenames = QtGui.QFileDialog.getOpenFileNames(self.pbMFF, "Select a mff file", "", "mff files (*.mff)", options= options) 
        filePath = mfffilenames[0] #[0] 
        print('%s' %filePath)
        
        hdr = read_mff_header.read_mff_header(filePath)       
        
        nC = hdr['nChans'] ; config.nC = nC
        nSamples = hdr['nSamples'] ; config.nSamples = nSamples
        nSamplesPre = hdr['nSamplesPre'] ; config.nSamplesPre = nSamplesPre
        nTrials = hdr['nTrials'] ; config.nTrials = nTrials
        srate = hdr['Fs'] ; config.srate = srate
        baseline = (nSamplesPre * 1000 / srate) ; config.baseline = baseline
        msSamples = np.arange(0, nSamples,1) * 1000/srate  - baseline 
        config.msSamples = msSamples
        summaryInfo = hdr['orig'] ; # config.nC = nC
        trialsName = summaryInfo['epochLabels'] ; config.trialsName = trialsName
        categoryName = list(set(trialsName)) ; config.categoryName = categoryName
        nCategory = len(categoryName) ; config.nCategory = nCategory        

        self.cbCategory.clear() 
        for i in range(nCategory):
            self.cbCategory.addItem(categoryName[i])

        self.editSamplesMSto.setText(_translate("FormMicro", "%d" %msSamples[-1],  None))

        data = read_mff_data.read_mff_data(filePath, 'epoch', 1, hdr['nTrials'], hdr)    

        SegStatus = hdr['orig']['epochSegStatus']  
        GoodSeg = []; BadSeg = []; #epochLabels = []; 
        for i in range(len(SegStatus)):
            if SegStatus[i] == 'bad' :  
                BadSeg.append(i)
            else :
                GoodSeg.append(i)
        
        nTrials = len(GoodSeg)

        ## average reference 
        H = np.identity(nC) - np.ones((nC, nC))/nC  
        s = np.zeros(data.shape) 
        if nTrials > 1 :
            for i in range(nTrials):
                s[:,:,i] = np.dot(H, data[:,:, GoodSeg[i]])
        else :
            s = np.dot(H, data)  

        whichCat = np.zeros((nTrials, nCategory))
        for j in range(nCategory):
            for i in range(nTrials):
                if (trialsName[i] == categoryName[j]):
                    whichCat[i,j] = 1
        whichCat = np.array(whichCat, dtype='i')

        config.s = s 
        config.whichCat = whichCat
        config.categoryName = categoryName
        
        print('nC:%d'  %nC)
        print('nSamples:%d'  %nSamples)
        print('nTrials:%d'  %nTrials)
        print('nCategory:%d' %nCategory)
        for i in range(nCategory):
            print('categoryName:%s' %categoryName[i])
            
        print("EEG data is loaded.")
        
        return s  

    ##-------------------------------------------------------------------------                
    def pbHeadModel_Callback(self):
        
        options = QtGui.QFileDialog.DontResolveSymlinks | QtGui.QFileDialog.ShowDirsOnly
        HMdir = QtGui.QFileDialog.getExistingDirectory(self.pbHeadModel,
                "Select a Head Model Directory.",
                self.pbHeadModel.text(), options=options)    
                
        import os 
        from braink import read_lfm                 
        import config                
        import glob 
    
        os.system("open /Applications/EAV/Reciprocity.app")

        lfm_fname = glob.glob(HMdir+'/*.lfm')
        if os.path.exists(HMdir+'/fdmForwardMatrixOriented'):          
            K  = read_lfm.forward(HMdir+'/fdmForwardMatrixOriented', config.nE)  
            print("fdmForwardMatrixOriented is loaded. ")            
        elif len(lfm_fname) > 0:  # elif os.path.exists(HMdir+'/Leadfield.lfm'): 
            K  = read_lfm.lfm(lfm_fname[0], config.nE)                          
            print("Leadfield.lfm is loaded. ")      
        else: 
            print("Either fdmForwardMatrixOriented or *.lfm are not loaded. ")            
        
        config.K = K    
        config.nV = K.shape[1]    
        
        return K

    ##-------------------------------------------------------------------------                
    def cbMethod_Callback(self, MethodIndex):
        import config 
        
        if MethodIndex == 0:
            config.MethodVal = 0
            print('MN is selected. ')
        elif MethodIndex == 1:
            config.MethodVal = 1
            print('sMN is selected. ')
        else:
            print("Inverse Method is not selected. ")
        return             
                       
    ##-------------------------------------------------------------------------                
    def pbEstimate_Callback(self):
        import config 
        from swcm import Imatrix, Lcurve
        import numpy as np      
        MethodVal = config.MethodVal 
        K = config.K
        s = config.s
        nC, nSamples, nTrials = s.shape 
        nV = K.shape[1]
                
#        logalpha = Lcurve.MN(s[:,0,0], K)
#        poweralpha = logalpha[0]
#        alpha = pow(10, poweralpha)                      
#        print("from Lcurve alpha=%3.2e" %alpha)

        alpha = 0.001
        if MethodVal == 0:
            Imat = Imatrix.MN(alpha, K) 
        elif MethodVal == 1:
            Imat = Imatrix.sMN(alpha, K) 
        else:
            print("Inverse Method is not selected. ")
            
        if nTrials > 1 :    
            sden = np.zeros((nV, nSamples, nTrials)) 
            for i in range(nTrials):             
                sden[:,:,i] = np.dot(Imat, s[:,:,i]) 
        else :
            sden = np.dot(Imat, s) 
            
        print("Source is estimated")                    
        config.sden = sden 
        
        return  sden 
            
    ##-------------------------------------------------------------------------                
    def cbCategory_Callback(self, Index):
        import config
        print('nCategory= %d'%config.nCategory)
        for i in range(config.nCategory): 
            if Index == i:
                config.CatVal = i
                print('Category[%d] (%s) is selected.' %(config.CatVal, config.categoryName[i]))


    def pbRunMicro_Callback(self):        
        import numpy as np 
        import config 
#        import configMicrostates 
        from microstates import Kmeans  
        import matplotlib.pyplot as plt
        plt.ion()
        from matplotlib import gridspec

        from pandas import DataFrame 
        import pandas as pd         
        gsn257 = pd.read_table('/Users/jesong1126/Python27/PermuStat/microstates/GSN257ToPy2.dat')  
        frame = DataFrame(gsn257 , columns=['ChLabel', 'X3', 'Y3', 'Z3','X2','Y2'])
        x2 = frame.values[:,4] 
        y2 = frame.values[:,5] 

        msSamples = config.msSamples
        s = config.s
        
        #---------------------------------------------------------------------                
        catName = categoryName[config.CatVal] 
        print('%s is selected' %catName) 

        NumK = self.editNumCluster.text()
        NumK = int(NumK)

        SamplesMS = self.editSamplesMS.text()
        SamplesMS = int(SamplesMS)
        nSamplesPre = np.where(msSamples ==  SamplesMS)[0]

        SamplesMSto = self.editSamplesMSto.text()
        SamplesMSto = int(SamplesMSto)
        nSamplesPost = np.where(msSamples ==  SamplesMSto)[0]
        
        Threshold = self.editThreshold.text()
        Threshold = float(Threshold)
 
        NumSim = self.editNumSim.text()
        NumSim = int(NumSim)

        if config.nTrials > 1:
            X = s[:,:,config.CatVal]  ; 
        else: 
            X = s ; 
        GFP = np.std(X, axis=0) 
        
        ##-----------------------------------------------------------------------------
        Data = X[:, nSamplesPre:nSamplesPost]
        nT = Data.shape[1] 
        KmeansOut = Kmeans.Kmeans(Data, NumK, NumSim, Threshold)
        
        KmeanId = np.array(KmeansOut['KmeanId'], 'i4') 
        KmeanId1 = np.array(KmeansOut['KmeanId1'], 'i4') 
        KmeanId2 = np.array(KmeansOut['KmeanId2'], 'i4')
        Tmaps3 = KmeansOut['Tmaps3'] 
        KmeanId3 = np.array(KmeansOut['KmeanId3'], 'i4') 
                
        Mycolorkeys = ('w','b', 'g', 'r', 'c', 'm', 'y', 'k', 'gray','purple','firebrick','darkgoldenrod','#afeeee','#8EBA42','#7A68A6','#56B4E9','#D55E00','b', 'g', 'r', 'c', 'm', 'y', 'k') 
        KmeanColorId = []  
        for i in range(nT):
            KmeanColorId.append(Mycolorkeys[KmeanId[i]])  
        
        KmeanIdColorId1 = []  
        for i in range(nT):  
            KmeanIdColorId1.append(Mycolorkeys[KmeanId1[i]])  
        
        KmeanIdColorId2 = []  
        for i in range(nT):  
            KmeanIdColorId2.append(Mycolorkeys[KmeanId2[i]])  
        
        KmeanIdColorId3 = [] # [None] * nT 
        for i in range(nT):  
            KmeanIdColorId3.append(Mycolorkeys[KmeanId3[i]])  
        
        startTP2 = np.zeros(NumK)  
        for i in range(NumK):
           startTP2[i] = np.where(KmeanId2==(i+1))[0][0]
        startTP2 = np.array(startTP2, dtype='i')   
            
        startTP3 = np.zeros(NumK) # np.array ([np.arange(k+1), np.zeros(k+1)]).reshape(2, k+1)
        for i in range(NumK):
           startTP3[i] = np.where(KmeanId3==(i+1))[0][0]
        startTP3 = np.array(startTP3, dtype='i')   
         
        plt.figure(figsize=(11.5, 8)) 
        gs = gridspec.GridSpec(6, NumK, height_ratios=[2,1,1,1,1,2], width_ratios = np.ones(NumK)) 
        ax0 = plt.subplot(gs[0,:])
        ax0.plot(msSamples, X.T)
        ax0.set_title('%s' %catName)
        ax0.set_xlim([msSamples[0], msSamples[-1]])
        ax1 = plt.subplot(gs[1,:])
        ax1.plot(msSamples, GFP, 'k')
        ax1.set_xlim([msSamples[0], msSamples[-1]])
        ax1.set_ylabel('Kmean')
        ax1.vlines(msSamples[nSamplesPre:nSamplesPost], [0], GFP[nSamplesPre:nSamplesPost], color=KmeanColorId)   
        ax2 = plt.subplot(gs[2,:])
        ax2.plot(msSamples, GFP, 'k')
        ax2.set_xlim([msSamples[0], msSamples[-1]])
        ax2.set_ylabel('%d ' %Threshold)
        ax2.vlines(msSamples[nSamplesPre:nSamplesPost], [0], GFP[nSamplesPre:nSamplesPost], color= KmeanIdColorId1) 
        ax3 = plt.subplot(gs[3,:])
        ax3.plot(msSamples, GFP, 'k')
        ax3.set_xlim([msSamples[0], msSamples[-1]])
        ax3.set_ylabel('Continuous')
        ax3.vlines(msSamples[nSamplesPre:nSamplesPost], [0], GFP[nSamplesPre:nSamplesPost], color= KmeanIdColorId2)  
        for i in range(NumK):
            ax3.text(msSamples[nSamplesPre+(startTP2[i])], GFP[nSamplesPre+(startTP2[i])]+0.1, ('%d' %(i+1)), fontsize=15)    
        ax4 = plt.subplot(gs[4,:])
        ax4.plot(msSamples, GFP, 'k')
        ax4.set_xlim([msSamples[0], msSamples[-1]])
        ax4.set_ylabel('Rename')
        ax4.vlines(msSamples[nSamplesPre:nSamplesPost], [0], GFP[nSamplesPre:nSamplesPost], color= KmeanIdColorId3)#, alpha= 0.2)    
        for i in range(NumK):
            ax4.text(msSamples[nSamplesPre+(startTP3[i])], GFP[nSamplesPre+(startTP3[i])]+0.1, ('%d' %(i+1)), fontsize=15)            
        for ii in range(NumK):
            axes = plt.subplot(gs[5, ii])
            axes.scatter(x2, y2, c=Tmaps3[:,ii], s=30, cmap=plt.get_cmap('seismic'), alpha= .5)
            axes.set_alpha(0.75)
            axes.set_xlabel(('MicroState %d' %(ii+1)))
            axes.set_xticks([])
            axes.set_yticks([])           
        plt.show()    
        

        from swcm import Imatrix #, Lcurve
        import numpy as np      
        MethodVal = config.MethodVal 
        K = config.K
        nV = K.shape[1]
                
#        logalpha = Lcurve.MN(Tmaps3[:,0], K)
#        poweralpha = logalpha[0]
#        alpha = pow(10, poweralpha)                      
#        print("from Lcurve alpha=%3.2e" %alpha)
       
        alpha = 0.001         

        if MethodVal == 0:
            Imat = Imatrix.MN(alpha, K) 
        elif MethodVal == 1:
            Imat = Imatrix.sMN(alpha, K) 
        else:
            print("Inverse Method is not selected. ")
            
        sdenTmaps = np.zeros((nV, NumK)) 
        for i in range(NumK):             
            sdenTmaps[:,i] = np.dot(Imat, Tmaps3[:,i]) 

        config.Tmaps = Tmaps3
        config.sdenTmaps = sdenTmaps

    ##-------------------------------------------------------------------------              
    def pbViewMicro_Callback(self):      
        import config 
#        import os 
#        os.system("open /Applications/EAV/Reciprocity.app")

        Tmaps = config.Tmaps; 
        sdenTmaps = config.sdenTmaps; 
        
        nC, NumK = Tmaps.shape  
        micI = self.editMicro.text()
        micI = int(micI)-1

        from RabbitMQ import Connection
        c = Connection.Connection('localhost')
        c.connect()
        c.sendEEGData(Tmaps[:, micI], "Microstates", "3:00pm")
#        if config.Orientedval == True:           
        c.sendOrientedData(abs(sdenTmaps[:, micI]), "Oriented", "2:00pm")        
#        else: 
#            c.sendTripleData(sdenTmaps[:, micI], "Triples", "2:00pm")        

        c.disconnect()       


    def pbClose_Callback(self):
        import matplotlib.pyplot as plt
        plt.close("all")


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    FormMicro = QtGui.QWidget()
    ui = Ui_FormMicro()
    ui.setupUi(FormMicro)
    FormMicro.show()
    sys.exit(app.exec_())

