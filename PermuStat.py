# -*- coding: utf-8 -*-

# Form implementation generated from reading ui file 'PermuStat.ui'
# Created: Mon Apr  6 15:58:15 2015
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

class Ui_FormStat(object):
    def setupUi(self, FormStat):
        FormStat.setObjectName(_fromUtf8("FormStat"))
        FormStat.resize(432, 650)
        
        ##---------------------------------------------------------------------
        self.tabPermuStat = QtGui.QTabWidget(FormStat)
        self.tabPermuStat.setEnabled(True)
        self.tabPermuStat.setGeometry(QtCore.QRect(20, 10, 390, 340))
        self.tabPermuStat.setMaximumSize(QtCore.QSize(600, 16777215))
        self.tabPermuStat.setFocusPolicy(QtCore.Qt.NoFocus)
        self.tabPermuStat.setTabPosition(QtGui.QTabWidget.North)
        self.tabPermuStat.setObjectName(_fromUtf8("tabPermuStat"))
        ##---------------------------------------------------------------------
        self.tabData = QtGui.QWidget()
        self.tabData.setObjectName(_fromUtf8("tabData"))
        self.pbMFF = QtGui.QPushButton(self.tabData)
        self.pbMFF.setGeometry(QtCore.QRect(20, 10, 341, 32))
        self.pbMFF.setObjectName(_fromUtf8("pbMFF"))

        self.frSource = QtGui.QFrame(self.tabData)
        self.frSource.setGeometry(QtCore.QRect(30, 50, 320, 250))
        self.frSource.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frSource.setFrameShadow(QtGui.QFrame.Raised)
        self.frSource.setObjectName(_fromUtf8("frSource"))
    
        self.pbHeadModel = QtGui.QPushButton(self.frSource)
        self.pbHeadModel.setGeometry(QtCore.QRect(10, 10, 300, 32))
        self.pbHeadModel.setObjectName(_fromUtf8("pbHeadModel"))
        
        self.frMethod = QtGui.QFrame(self.frSource)
        self.frMethod.setGeometry(QtCore.QRect(20, 50, 280, 40))
        self.frMethod.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frMethod.setFrameShadow(QtGui.QFrame.Raised)
        self.frMethod.setObjectName(_fromUtf8("frMethod"))
        
        self.label_Method = QtGui.QLabel(self.frMethod)
        self.label_Method.setGeometry(QtCore.QRect(20, 10, 71, 16))
        self.label_Method.setObjectName(_fromUtf8("label_Method"))
        self.cbMethod = QtGui.QComboBox(self.frMethod)
        self.cbMethod.setGeometry(QtCore.QRect(150, 8, 120, 26))
        self.cbMethod.setObjectName(_fromUtf8("cbMethod"))
        self.cbMethod.addItems(["MN", "sMN"])
#        self.cbMethod.addItem(_fromUtf8("MN"))
#        self.cbMethod.addItem(_fromUtf8("sMN"))
        
        self.frReg = QtGui.QFrame(self.frSource)
        self.frReg.setGeometry(QtCore.QRect(20, 100, 280, 71))
        self.frReg.setFrameShape(QtGui.QFrame.StyledPanel)
        self.frReg.setFrameShadow(QtGui.QFrame.Raised)
        self.frReg.setObjectName(_fromUtf8("frReg"))
        
        self.rbTikhonov = QtGui.QRadioButton(self.frReg)
        self.rbTikhonov.setGeometry(QtCore.QRect(20, 10, 91, 20))
        self.rbTikhonov.setObjectName(_fromUtf8("rbTikh"))
        self.rbTikhonov.setChecked(True)

        self.rbLcurve = QtGui.QRadioButton(self.frReg)
        self.rbLcurve.setGeometry(QtCore.QRect(20, 40, 81, 20))
        self.rbLcurve.setObjectName(_fromUtf8("rbLcurve")) 
        self.label_alpha_Tikh = QtGui.QLabel(self.frReg)
        self.label_alpha_Tikh.setGeometry(QtCore.QRect(110, 10, 81, 16))
        self.label_alpha_Tikh.setObjectName(_fromUtf8("label_alpha_Tikh"))
        self.leTikhonov = QtGui.QLineEdit(self.frReg)
        self.leTikhonov.setGeometry(QtCore.QRect(200, 10, 60, 21))
        self.leTikhonov.setObjectName(_fromUtf8("leTikhonov"))
        
        self.pbEstimate = QtGui.QPushButton(self.frSource)
        self.pbEstimate.setGeometry(QtCore.QRect(10, 180, 300, 32))
        self.pbEstimate.setObjectName(_fromUtf8("pbEstimate"))
        
        self.EstimateBar = QtGui.QProgressBar(self.frSource)
        self.EstimateBar.setGeometry(15, 215, 290, 25)
        self.EstimateBar.setObjectName(_fromUtf8("EstimateBar"))
        self.EstimateBar.setMinimum(0)
        self.EstimateBar.setMaximum(100)   
        
        self.tabPermuStat.addTab(self.tabData, _fromUtf8(""))
        
        
        ##---------------------------------------------------------------------
        self.tabSource = QtGui.QWidget()
        self.tabSource.setObjectName(_fromUtf8("tabSource"))
        self.rbScalp = QtGui.QRadioButton(self.tabSource)
        self.rbScalp.setGeometry(QtCore.QRect(100, 30, 170, 20))
        self.rbScalp.setChecked(True)
        self.rbScalp.setObjectName(_fromUtf8("rbScalp"))
        self.rbSource = QtGui.QRadioButton(self.tabSource)
        self.rbSource.setGeometry(QtCore.QRect(100, 70, 170, 20))
        self.rbSource.setObjectName(_fromUtf8("rbSource"))
        self.tabPermuStat.addTab(self.tabSource, _fromUtf8(""))
        ##---------------------------------------------------------------------
        self.tabLevel = QtGui.QWidget()
        self.tabLevel.setObjectName(_fromUtf8("tabLevel"))
        self.rbIndividual = QtGui.QRadioButton(self.tabLevel)
        self.rbIndividual.setGeometry(QtCore.QRect(80, 30, 170, 20))
        self.rbIndividual.setChecked(True)
        self.rbIndividual.setObjectName(_fromUtf8("rbIndividual"))
        self.rbGroup = QtGui.QRadioButton(self.tabLevel)
        self.rbGroup.setGeometry(QtCore.QRect(80, 80, 170, 20))
        self.rbGroup.setObjectName(_fromUtf8("rbGroup"))
        self.tabPermuStat.addTab(self.tabLevel, _fromUtf8(""))
        ##---------------------------------------------------------------------
        self.tabParametric = QtGui.QWidget()
        self.tabParametric.setObjectName(_fromUtf8("tabParametric"))
        self.rbPara = QtGui.QRadioButton(self.tabParametric)
        self.rbPara.setGeometry(QtCore.QRect(100, 30, 170, 20))
        self.rbPara.setChecked(True)
        self.rbPara.setObjectName(_fromUtf8("rbPara"))
        self.rbNonPara = QtGui.QRadioButton(self.tabParametric)
        self.rbNonPara.setGeometry(QtCore.QRect(100, 70, 170, 20))
        self.rbNonPara.setObjectName(_fromUtf8("rbNonPara"))
        self.tabPermuStat.addTab(self.tabParametric, _fromUtf8(""))
        ##---------------------------------------------------------------------
        self.tabStat = QtGui.QWidget()
        self.tabStat.setObjectName(_fromUtf8("tabStat"))
        self.label_SigLev = QtGui.QLabel(self.tabStat)
        self.label_SigLev.setGeometry(QtCore.QRect(40, 30, 110, 16))
        self.label_SigLev.setObjectName(_fromUtf8("label"))
        self.label_NumSim = QtGui.QLabel(self.tabStat)
        self.label_NumSim.setGeometry(QtCore.QRect(40, 110, 110, 16))
        self.label_NumSim.setObjectName(_fromUtf8("label_NumSim"))
        self.label_Output = QtGui.QLabel(self.tabStat)
        self.label_Output.setGeometry(QtCore.QRect(40, 70, 110, 16))
        self.label_Output.setObjectName(_fromUtf8("label_Output"))
        self.editSigLevel = QtGui.QLineEdit(self.tabStat)
        self.editSigLevel.setGeometry(QtCore.QRect(170, 30, 100, 21))
        self.editSigLevel.setObjectName(_fromUtf8("editSigLevel"))
        self.editOutput = QtGui.QLineEdit(self.tabStat)
        self.editOutput.setGeometry(QtCore.QRect(170, 70, 100, 21))
        self.editOutput.setObjectName(_fromUtf8("editOutput"))
        self.editNumSim = QtGui.QLineEdit(self.tabStat)
        self.editNumSim.setGeometry(QtCore.QRect(170, 110, 100, 21))
        self.editNumSim.setObjectName(_fromUtf8("editNumSim"))
        self.tabPermuStat.addTab(self.tabStat, _fromUtf8(""))
        ##---------------------------------------------------------------------
        self.tabTest = QtGui.QWidget()
        self.tabTest.setObjectName(_fromUtf8("tabTest"))
        self.rbOneT = QtGui.QRadioButton(self.tabTest)
        self.rbOneT.setGeometry(QtCore.QRect(20, 20, 141, 20))
        self.rbOneT.setChecked(True)
        self.rbOneT.setCheckable(True)
        self.rbOneT.setObjectName(_fromUtf8("rbOneT"))
        self.rbTwoT = QtGui.QRadioButton(self.tabTest)
        self.rbTwoT.setGeometry(QtCore.QRect(20, 80, 141, 20))
        self.rbTwoT.setObjectName(_fromUtf8("rbTwoT"))
        self.rbTwoT.setCheckable(True)
        self.rbPairT = QtGui.QRadioButton(self.tabTest)
        self.rbPairT.setGeometry(QtCore.QRect(20, 170, 141, 20))
        self.rbPairT.setObjectName(_fromUtf8("rbPairT"))
        self.rbPairT.setCheckable(True)
        
        self.rbAnova = QtGui.QRadioButton(self.tabTest)
        self.rbAnova.setGeometry(QtCore.QRect(20, 260, 141, 20))
        self.rbAnova.setObjectName(_fromUtf8("rbAnova"))
        self.rbAnova.setCheckable(True)
        
        self.label_OneT = QtGui.QLabel(self.tabTest)
        self.label_OneT.setGeometry(QtCore.QRect(70, 50, 81, 16))
        self.label_OneT.setObjectName(_fromUtf8("label_OneT"))
        self.label_TwoT1 = QtGui.QLabel(self.tabTest)
        self.label_TwoT1.setGeometry(QtCore.QRect(70, 105, 81, 16))
        self.label_TwoT1.setObjectName(_fromUtf8("label_TwoT1"))
        self.label_TwoT2 = QtGui.QLabel(self.tabTest)
        self.label_TwoT2.setGeometry(QtCore.QRect(70, 135, 81, 16))
        self.label_TwoT2.setObjectName(_fromUtf8("label_TwoT2"))
        self.label_PairT1 = QtGui.QLabel(self.tabTest)
        self.label_PairT1.setGeometry(QtCore.QRect(70, 195, 81, 16))
        self.label_PairT1.setObjectName(_fromUtf8("label_PairT1"))
        self.label_PairT2 = QtGui.QLabel(self.tabTest)
        self.label_PairT2.setGeometry(QtCore.QRect(70, 225, 81, 16))
        self.label_PairT2.setObjectName(_fromUtf8("label_PairT2"))
        self.cbOneTCat = QtGui.QComboBox(self.tabTest)
        self.cbOneTCat.setGeometry(QtCore.QRect(180, 45, 141, 26))
        self.cbOneTCat.setObjectName(_fromUtf8("cbOneTCat"))
        self.cbTwoTCat1 = QtGui.QComboBox(self.tabTest)
        self.cbTwoTCat1.setGeometry(QtCore.QRect(180, 100, 141, 26))
        self.cbTwoTCat1.setObjectName(_fromUtf8("cbTwoTCat1"))
        self.cbTwoTCat2 = QtGui.QComboBox(self.tabTest)
        self.cbTwoTCat2.setGeometry(QtCore.QRect(180, 130, 141, 26))
        self.cbTwoTCat2.setObjectName(_fromUtf8("cbTwoTCat2"))
        self.cbPairTCat1 = QtGui.QComboBox(self.tabTest)
        self.cbPairTCat1.setGeometry(QtCore.QRect(180, 190, 141, 26))
        self.cbPairTCat1.setObjectName(_fromUtf8("cbPairTCat1"))
        self.cbPairTCat2 = QtGui.QComboBox(self.tabTest)
        self.cbPairTCat2.setGeometry(QtCore.QRect(180, 220, 141, 26))
        self.cbPairTCat2.setObjectName(_fromUtf8("cbPairTCat2"))
        self.tabPermuStat.addTab(self.tabTest, _fromUtf8(""))
        ##---------------------------------------------------------------------                
        self.pbRunStat = QtGui.QPushButton(FormStat)
        self.pbRunStat.setGeometry(QtCore.QRect(20, 355, 390, 45))
        self.pbRunStat.setObjectName(_fromUtf8("pbRunStat"))
        ##---------------------------------------------------------------------                
        self.pbSaveStat = QtGui.QPushButton(FormStat)
        self.pbSaveStat.setGeometry(QtCore.QRect(20, 400, 390, 45))
        self.pbSaveStat.setObjectName(_fromUtf8("pbSaveStat"))
        ##---------------------------------------------------------------------
        self.formGroupBox = QtGui.QGroupBox(FormStat)
        self.formGroupBox.setGeometry(QtCore.QRect(25, 450, 385, 145))
        self.formGroupBox.setObjectName(_fromUtf8("formGroupBox"))
        self.label_Time = QtGui.QLabel(self.formGroupBox)
        self.label_Time.setGeometry(QtCore.QRect(80, 15, 100, 16))
        self.label_Time.setObjectName(_fromUtf8("label_Time"))
        self.editTime = QtGui.QLineEdit(self.formGroupBox)
        self.editTime.setGeometry(QtCore.QRect(200, 15, 100, 21))
        self.editTime.setObjectName(_fromUtf8("editTime"))
        self.label_Channel = QtGui.QLabel(self.formGroupBox)
        self.label_Channel.setGeometry(QtCore.QRect(80, 45, 100, 16))
        self.label_Channel.setObjectName(_fromUtf8("label_Channel"))
        self.editChannel = QtGui.QLineEdit(self.formGroupBox)
        self.editChannel.setGeometry(QtCore.QRect(200, 45, 100, 21))
        self.editChannel.setObjectName(_fromUtf8("editChannel"))
        self.label_Dipole = QtGui.QLabel(self.formGroupBox)
        self.label_Dipole.setGeometry(QtCore.QRect(80, 75, 100, 16))
        self.label_Dipole.setObjectName(_fromUtf8("label_Dipole"))
        self.editDipole = QtGui.QLineEdit(self.formGroupBox)
        self.editDipole.setGeometry(QtCore.QRect(200, 75, 100, 21))
        self.editDipole.setObjectName(_fromUtf8("editDipole"))
        self.pbEAVStat = QtGui.QPushButton(self.formGroupBox)
        self.pbEAVStat.setGeometry(QtCore.QRect(40, 110, 300, 32))
        self.pbEAVStat.setObjectName(_fromUtf8("pbEAVStat"))         
        
        ##---------------------------------------------------------------------                
        self.pbCloseStat = QtGui.QPushButton(FormStat)
        self.pbCloseStat.setGeometry(QtCore.QRect(20, 600, 390, 45))
        self.pbCloseStat.setObjectName(_fromUtf8("pbCloseStat"))

        ##---------------------------------------------------------------------                            
        self.pbMFF.clicked.connect(self.pbMFF_Callback)
        self.pbHeadModel.clicked.connect(self.pbHeadModel_Callback) 
        self.cbMethod.activated.connect(self.cbMethod_Callback)
        
        self.pbEstimate.clicked.connect(self.pbEstimate_Callback) 
#        self.pbEstimate.clicked.connect(self.doAction)

        self.cbOneTCat.activated.connect(self.cbOneTCat_Callback)    
        self.cbTwoTCat1.activated.connect(self.cbTwoTCat1_Callback)  
        self.cbTwoTCat2.activated.connect(self.cbTwoTCat2_Callback)
        self.cbPairTCat1.activated.connect(self.cbPairTCat1_Callback)
        self.cbPairTCat2.activated.connect(self.cbPairTCat2_Callback)
        
        self.pbRunStat.clicked.connect(self.pbRunStat_Callback)
        self.pbSaveStat.clicked.connect(self.pbSaveStat_Callback)
        self.pbEAVStat.clicked.connect(self.pbEAVStat_Callback)
        self.pbCloseStat.clicked.connect(self.pbCloseStat_Callback)
        self.pbCloseStat.clicked.connect(FormStat.close)
        
        
        self.retranslateUi(FormStat)
        self.tabPermuStat.setCurrentIndex(0)
        QtCore.QMetaObject.connectSlotsByName(FormStat)

    ##-------------------------------------------------------------------------                
    def retranslateUi(self, FormStat):
        FormStat.setWindowTitle(_translate("FormStat", "PermuStat", None))
        self.pbMFF.setText(_translate("FormStat", "MFF", None))
        self.pbHeadModel.setText(_translate("FormStat", "Head Model", None))
        self.label_Method.setText(_translate("FormStat", "Method : ", None))
        self.rbTikhonov.setText(_translate("FormStat", "Tikhonov :", None))
        self.rbLcurve.setText(_translate("FormStat", "L-curve ", None))
        self.label_alpha_Tikh.setText(_translate("FormStat", "alpha = 10^", None))
        self.leTikhonov.setText(_translate("FormStat", "-2", None))
        self.cbMethod.setItemText(0, _translate("FormStat", "MN", None))
        self.cbMethod.setItemText(1, _translate("FormStat", "sMN", None))
        self.pbEstimate.setText(_translate("FormStat", "Estimate", None))
        self.tabPermuStat.setTabText(self.tabPermuStat.indexOf(self.tabData), _translate("FormStat", "Data", None))
        self.rbScalp.setText(_translate("FormStat", "Scalp", None))
        self.rbSource.setText(_translate("FormStat", "Source", None))
        self.tabPermuStat.setTabText(self.tabPermuStat.indexOf(self.tabSource), _translate("FormStat", "Source", None))
        self.rbIndividual.setText(_translate("FormStat", "Individual (First Level)", None))
        self.rbGroup.setText(_translate("FormStat", "Group (Second level)", None))
        self.tabPermuStat.setTabText(self.tabPermuStat.indexOf(self.tabLevel), _translate("FormStat", "Level", None))
        self.rbPara.setText(_translate("FormStat", "Parametric", None))
        self.rbNonPara.setText(_translate("FormStat", "NonParametric", None))
        self.tabPermuStat.setTabText(self.tabPermuStat.indexOf(self.tabParametric), _translate("FormStat", "Parametric", None))
        self.label_SigLev.setText(_translate("FormStat", "Significant Level", None))
        self.label_NumSim.setText(_translate("FormStat", "Num. Simulation", None))
        self.label_Output.setText(_translate("FormStat", "Output Suffix", None))
        self.editSigLevel.setText(_translate("FormStat", "0.05", None))
        self.editOutput.setText(_translate("FormStat", "out", None))
        self.editNumSim.setText(_translate("FormStat", "100", None))
        self.tabPermuStat.setTabText(self.tabPermuStat.indexOf(self.tabStat), _translate("FormStat", "Stat", None))
        self.rbOneT.setText(_translate("FormStat", "One Sample T-test", None))
        self.rbTwoT.setText(_translate("FormStat", "Two Sample T-test", None))
        self.rbPairT.setText(_translate("FormStat", "Paired T-test", None))
        self.rbAnova.setText(_translate("FormStat", "ANOVA", None))
        self.label_OneT.setText(_translate("FormStat", "Category", None))
        self.label_TwoT1.setText(_translate("FormStat", "Category 1", None))
        self.label_TwoT2.setText(_translate("FormStat", "Category 2", None))
        self.label_PairT1.setText(_translate("FormStat", "Category 1", None))
        self.label_PairT2.setText(_translate("FormStat", "Category 2", None))
        self.tabPermuStat.setTabText(self.tabPermuStat.indexOf(self.tabTest), _translate("FormStat", "Test", None))
        self.pbRunStat.setText(_translate("FormStat", "Run", None))
        self.pbSaveStat.setText(_translate("FormStat", "Reciprocity", None))
        self.label_Time.setText(_translate("FormStat", "Time (ms)", None))
        self.label_Channel.setText(_translate("FormStat", "Channel", None))
        self.label_Dipole.setText(_translate("FormStat", "Dipole", None))
        self.editTime.setText(_translate("FormStat", "0", None))
        self.editChannel.setText(_translate("FormStat", "0", None))
        self.editDipole.setText(_translate("FormStat", "0", None))
        self.pbEAVStat.setText(_translate("FormStat", "View", None))
        self.pbCloseStat.setText(_translate("FormStat", "Close", None))
 
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

        for i in range(nCategory):
            self.cbOneTCat.addItem(categoryName[i])
            self.cbTwoTCat1.addItem(categoryName[i])
            self.cbTwoTCat2.addItem(categoryName[i])
            self.cbPairTCat1.addItem(categoryName[i])
            self.cbPairTCat2.addItem(categoryName[i])

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
   
        lfm_fname = glob.glob(HMdir+'/*.lfm')
        if os.path.exists(HMdir+'/fdmForwardMatrixOriented'):          
            K  = read_lfm.forward(HMdir+'/fdmForwardMatrixOriented', config.nE)  
            print("fdmForwardMatrixOriented is loaded. ")            
        elif len(lfm_fname) > 0:  # elif os.path.exists(HMdir+'/Leadfield.lfm'): 
            K  = read_lfm.lfm(lfm_fname[0], config.nE)                          
            print("Leadfield.lfm is loaded. ")      
        else: 
            print("Either fdmForwardMatrixOriented or *.lfm are not loaded. ")            
        
#        nV = K.shape[1]
        config.K = K    
        config.nV = K.shape[1]    
        
        return K

    ##-------------------------------------------------------------------------                
    def cbMethod_Callback(self, MethodIndex):
        import config 
        
        if MethodIndex == 0:
            config.MethodVal = 0
#            print(MethodIndex)
            print('MN is selected. ')
        elif MethodIndex == 1:
            config.MethodVal = 1
#            print(MethodIndex)
            print('sMN is selected. ')
        else:
            print("Inverse Method is not selected. ")

        return             
                       
    ##-------------------------------------------------------------------------                
    def pbEstimate_Callback(self):
        import config 
        from swcm import Imatrix, Lcurve
        import numpy as np      
#        import pdb
#        pdb.pm() 
        MethodVal = config.MethodVal 
        K = config.K
        s = config.s
        nC, nSamples, nTrials = s.shape 
        nV = K.shape[1]
                
        if self.rbTikhonov.isChecked():
            poweralpha = self.leTikhonov.text()
            poweralpha = int(poweralpha)
            alpha = pow(10, poweralpha)  
            print("alpha = %3.2e "  %alpha)
        else :         
            logalpha = Lcurve.MN(s[:,0,0], K)
            poweralpha = logalpha[0]
            alpha = pow(10, poweralpha)                      
            print("from Lcurve alpha=%3.2e" %alpha)

        if MethodVal == 0:
            Imat = Imatrix.MN(alpha, K) 
        elif MethodVal == 1:
            Imat = Imatrix.sMN(alpha, K) 
        else:
            print("Inverse Method is not selected. ")
            
        sden = np.zeros((nV, nSamples, nTrials)) 
        for i in range(nTrials):             
            sden[:,:,i] = np.dot(Imat, s[:,:,i]) 
        
        print("Source is estimated")                    
        config.sden = sden 
        
        return  sden 
            
#    def timerEvent(self, event):  
#        self.timer = QtCore.QBasicTimer()        
#        self.step = 0 
#        if self.step >= 100:            
#            self.timer.stop()
#            self.pbEstimate.setText('Finished')
#            return            
#        self.step = self.step + 1
#        self.EstimateBar.setValue(self.step)

#    def doAction(self, val):      
#        self.timer = QtCore.QBasicTimer()        
#        self.step = 0 
#        self.EstimateBar.setValue(val)
##        perct = "{0}%".format(val)
#        if self.timer.isActive():
#            self.timer.stop()
##            self.EstimateBar.setText('Start')            
#        else:
#            self.timer.start(100, self.EstimateBar)
##            self.EstimateBar.setText('Stop')
            

    ##-------------------------------------------------------------------------                
    def cbOneTCat_Callback(self, Index):
        import config
        print('nCategory= %d'%config.nCategory)
        for i in range(config.nCategory): 
            if Index == i:
                config.OneTCatVal = i
                print('Category[%d] (%s) is selected.' %(config.OneTCatVal, config.categoryName[i]))

    ##-------------------------------------------------------------------------                
    def cbTwoTCat1_Callback(self, Index):
        import config
        for i in range(config.nCategory): 
            if Index == i:
                config.TwoTCatVal1 = i
                print('Category[%d] (%s) is selected.' %(config.TwoTCatVal1, config.categoryName[i]))

    ##-------------------------------------------------------------------------                
    def cbTwoTCat2_Callback(self, Index):
        import config
        for i in range(config.nCategory): 
            if Index == i:
                config.TwoTCatVal2 = i    
                print('Category[%d] (%s) is selected.' %(config.TwoTCatVal2, config.categoryName[i]))

    ##-------------------------------------------------------------------------
    def cbPairTCat1_Callback(self, Index):
        import config
        for i in range(config.nCategory): 
            if Index == i:
                config.PairTCatVal1 = i
                print('Category[%d] (%s) is selected.' %(config.PairTCatVal1, config.categoryName[i]))

    ##-------------------------------------------------------------------------                
    def cbPairTCat2_Callback(self, Index):
        import config
        for i in range(config.nCategory): 
            if Index == i:
                config.PairTCatVal2 = i    
                print('Category[%d] (%s) is selected.' %(config.PairTCatVal2, config.categoryName[i]))

    ##-------------------------------------------------------------------------                
    def pbRunStat_Callback(self):        
        import numpy as np 
        import config 
        from permustat import permustat_fun, perm_plot 

        msSamples = config.msSamples
        s = config.s
        sden = config.sden
        categoryName = config.categoryName
        whichCat = config.whichCat 
        nCategory = config.nCategory
        
        #---------------------------------------------------------------------                
        SigLevel = self.editSigLevel.text()
        SigLevel = float(SigLevel)
        print("Significant level is %3.2f" %SigLevel)      

        NumSim = self.editNumSim.text()
        NumSim = int(NumSim)
         
        if self.rbPara.isChecked():
            config.ParaVal = 0
            print("Parametric is selected")
        else:
            config.ParaVal = 1
            print("Non-Parametric is selected and Num Sim is %d" %NumSim)      

        ##---------------------------------------------------------------------
        if self.rbOneT.isChecked():
            catName = categoryName[config.OneTCatVal] 
            print("One sample T-test:%s" %catName)
            if self.rbScalp.isChecked():
                X = s[:,:, (whichCat[:, config.OneTCatVal]>0)]        
            else: 
                X = sden[:,:, (whichCat[:, config.OneTCatVal]>0)]    
            XMean = np.mean(X, axis=2)
            XSigma = np.std(X, axis=2)   
            statOut = permustat_fun.one_sample_t(X, SigLevel)         
            perm_plot.oneT(XMean, XSigma, statOut['T'], statOut['nSigS'], catName, msSamples)
            config.Stat = statOut['T']                  

        elif self.rbTwoT.isChecked():
            catName1 = categoryName[config.TwoTCatVal1]
            catName2 = categoryName[config.TwoTCatVal2]
            print('Two sample T-test: %s vs %s' %(catName1, catName2)) 
            if self.rbScalp.isChecked():
                X1 = s[:,:, (whichCat[:, config.TwoTCatVal1]>0)]        
                X2 = s[:,:, (whichCat[:, config.TwoTCatVal2]>0)]        
            else: 
                X1 = sden[:,:, (whichCat[:, config.TwoTCatVal1]>0)]        
                X2 = sden[:,:, (whichCat[:, config.TwoTCatVal2]>0)]        
                
            XMean1 = np.mean(X1, axis=2); XSigma1 = np.std(X1, axis=2);
            XMean2 = np.mean(X2, axis=2); XSigma2 = np.std(X2, axis=2);   
            statOut1 = permustat_fun.one_sample_t(X1, SigLevel/2) 
            statOut2 = permustat_fun.one_sample_t(X2, SigLevel/2) 
            
            if self.rbPara.isChecked():
                statOut = permustat_fun.two_sample_t(X1, X2, SigLevel) 
                perm_plot.twoT(XMean1, XSigma1, statOut1['T'], statOut1['nSigS'], catName1, 
                               XMean2, XSigma2, statOut2['T'], statOut2['nSigS'], catName2, statOut['T'], statOut['nSigS'], msSamples)
            else:
                statOut = permustat_fun.two_sample_permu(X1, X2, SigLevel/2, NumSim)
                logPvalST = np.log(statOut['Pval'])
                perm_plot.permT(statOut['T'], logPvalST, msSamples) 
            
            config.Stat = statOut['T']                  

        elif self.rbPairT.isChecked():
            catName1 = categoryName[config.PairTCatVal1]
            catName2 = categoryName[config.PairTCatVal2]
            print('Pair T-test: %s vs %s' %(catName1, catName2)) 
            if self.rbScalp.isChecked():
                X1 = s[:,:, (whichCat[:, config.PairTCatVal1]>0)]        
                X2 = s[:,:, (whichCat[:, config.PairTCatVal2]>0)]        
            else: 
                X1 = sden[:,:, (whichCat[:, config.PairTCatVal1]>0)]        
                X2 = sden[:,:, (whichCat[:, config.PairTCatVal2]>0)]        
            
            XMean1 = np.mean(X1, axis=2); XSigma1 = np.std(X1, axis=2);
            XMean2 = np.mean(X2, axis=2); XSigma2 = np.std(X2, axis=2);   
            statOut1 = permustat_fun.one_sample_t(X1, SigLevel/2) 
            statOut2 = permustat_fun.one_sample_t(X2, SigLevel/2) 
                
            if self.rbPara.isChecked():
                statOut = permustat_fun.pairt(X1, X2, SigLevel) 
                perm_plot.pairT(XMean1, XSigma1, statOut1['T'], statOut1['nSigS'], catName1, 
                               XMean2, XSigma2, statOut2['T'], statOut2['nSigS'], catName2, statOut['T'], statOut['nSigS'], msSamples)
            else:
                statOut = permustat_fun.pairt_permu(X1, X2, SigLevel/2, NumSim)
                logPvalST = np.log(statOut['Pval'])
                perm_plot.permPairT(statOut['T'], logPvalST, msSamples) 
            
            config.Stat = statOut['T']                  
            
        elif self.rbAnova.isChecked():
            if self.rbScalp.isChecked():
                statOut = permustat_fun.anova_f(s, whichCat, nCategory, SigLevel)
            else: 
                statOut = permustat_fun.anova_f(sden, whichCat, nCategory, SigLevel)
            perm_plot.anova(statOut['F'], statOut['Pval'], statOut['Sig'], msSamples)              
            config.Stat = statOut['F']

        config.Pval = statOut['Pval']  
        config.Sig = statOut['Sig'] 
        config.StatOut = statOut

        return statOut

    ##-------------------------------------------------------------------------              
    def pbSaveStat_Callback(self):      
        import os
        os.system("open /Applications/EAV/Mimir.app")
#        os.system("open /Applications/EAV/Reciprocity.app")

#        import config         
#        Stat = config.Stat         
#        Sig = config.Sig                 
#        Pval = config.Pval                 
#        
#        OutputSuffix = self.editOutput.text()
#        print("OutputSuffix is %s" %OutputSuffix)      
#
#        f=open("Stat_"+OutputSuffix+".bin","wb")
#        f.write(Stat)
#        f.close() 
#        
#        f=open("PvalST_"+OutputSuffix+".bin","wb")
#        f.write(Pval)
#        f.close() 
#        
#        f=open("SigST_"+OutputSuffix+".bin","wb")
#        f.write(Sig)
#        f.close() 

    ##-------------------------------------------------------------------------
    def pbEAVStat_Callback(self):      
        import config         
        import matplotlib.pyplot as plt 
        from RabbitMQ import Connection 
        import numpy as np
        from pandas import DataFrame 
        import pandas as pd 
        pd.options.display.float_format = '{:,.3f}'.format

        statOut = config.StatOut
        XT = config.Stat; Sig = config.Sig; # Pval = config.Pval                
        msSamples = config.msSamples  
        
        EAVtime = self.editTime.text(); EAVtime = int(EAVtime)
#        print("Time is %d" %EAVtime)      
        
        SampleTime = np.where(msSamples == EAVtime)[0]
#        print("SampleTime is %d" %SampleTime)
        
        EAVch = self.editChannel.text(); EAVch = int(EAVch)
#        print("Channel is %d" %EAVch)      
        
        EAVdipole = self.editDipole.text(); EAVdipole = int(EAVdipole)
#        print("Dipole is %d" %EAVdipole)      
        
        ##---------------------------------------------------------------------      
        if self.rbScalp.isChecked():
            config.SourceVal = 0
            print("Scalp is selected")
            i = EAVch ; 
            figCh = plt.figure()
            ax = figCh.add_subplot(111) 
            ax.plot(msSamples, XT[EAVch, :])
            ax.axvline(EAVtime, color='r')
            ax.set_xlim([msSamples[0],msSamples[-1]])            
            ax.set_title('Electrode %d' %EAVch)
            plt.show() 
        else:
            config.SourceVal = 1
            print("Source is selected") 
            i = EAVdipole ; 
            figCh = plt.figure()
            ax = figCh.add_subplot(111) 
            ax.plot(msSamples, XT[EAVdipole, :])
            ax.axvline(EAVtime, color='r')
            ax.set_xlim([msSamples[0],msSamples[-1]])            
            ax.set_title('Dipole %d' %EAVdipole)
            plt.show()
        
        j = SampleTime;  
        if self.rbAnova.isChecked():
            tbl_one = [[statOut['df_b'], float(statOut['SS_b'][i,j]), float(statOut['MS_b'][i,j]), float(statOut['F'][i,j]), 
                      float(statOut['Pval'][i,j]), ('***' if statOut['Sig'][i,j] == 1 else '')], 
                      [statOut['df_w'], float(statOut['SS_w'][i,j]), float(statOut['MS_w'][i,j]), '', '', ''], 
                      [statOut['df_b']+statOut['df_w'], float(statOut['SS_t'][i,j]), '', '', '', '']]
            anova = DataFrame(tbl_one, index=['Between', 'Within', 'Total'], columns=[ 'DF' ,'SS' ,'MS', 'F' , 'Pval', 'Sig'])
            print anova
        
        c = Connection.Connection('localhost')
        c.connect()
        if self.rbScalp.isChecked():
            c.sendEEGData(XT[:,SampleTime], "StatEEG", "")
            c.sendEEGData(Sig[:,SampleTime], "SigEEG", "")
        else: 
            c.sendOrientedData(XT[:, SampleTime], "StatSource", "")
            c.sendOrientedData(Sig[:, SampleTime], "SigSource", "")         
            
        c.disconnect()
        return 
    
    ##-------------------------------------------------------------------------
    def pbCloseStat_Callback(self):
        import matplotlib.pyplot as plt
        plt.close("all") 


if __name__ == "__main__":
    app = QtGui.QApplication(sys.argv)
    FormStat = QtGui.QWidget()
    ui = Ui_FormStat()
    ui.setupUi(FormStat)
    FormStat.show()
    sys.exit(app.exec_())


