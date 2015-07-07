#import numpy as np 

nE = 256; nC = 257; s = []; K = []; sden = [];
MethodVal=0; alpha = 0.01; CondNum = 1; CondTxt = 1;  

srate = 1000; nSamples = 0; nSamplesPre = 0; baseline = 0; msSamples = 0; 
nTrials = 1; trialsName = []; nCategory = 1; categoryName = [];
whichCat = []; 

LevelVal = 0; TestVal = 0; SourceVal = 0; ParaVal = 0;  

OneTCatVal = 0;  TwoTCatVal1 = 0; TwoTCatVal2 = 1; 
PairTCatVal1 = 0; PairTCatVal2 = 1; 
Stat = []; Pval =[]; Sig =[]; StatOut = [];

CatVal = 0;

#thresholdcutoff = 0 
#startSplTxt = 1 
#endSplTxt = 1 
#MovieVal = 1  
#TopoVal = 1
#HMdir = 'Home'
#nbase = 100 
#nem = 30 
#sGFP = []

#sdenMean = []
#sdenSigma = []
#sdenT = []
#sMean = []
#sSigma = []
#sT = []

#Spl = 1 
#startMS = 0
#endMS = 1 
#EpochVal = 1 
#InfCorVal = 1

#s = np.empty([1,1,1])
#K = np.empty([1,1,1]) #[]
#sden = np.empty([1,1,1])
#RegVal = 1

#HMval = 1
#Atlasval=1
#Orientedval=1 
#Editedval=1 
#Scenarioval=1 
#SensorTypeSuffixval=1 
#scenario = 'MRI' 
#sensorTypeSuffix = 'GPS'
