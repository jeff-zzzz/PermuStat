
def read_mff_header(filePath):
#    import numpy as np
    from xml.dom.minidom import parse
    from mff.mff_getSummaryInfo import mff_getSummaryInfo

    ##---------------------------------------------------------------------------
    summaryInfo = mff_getSummaryInfo(filePath) 
    
    ##---------------------------------------------------------------------------
    # Pull header info from the summary info. 
    nSamplesPre = 0 
    if summaryInfo['epochType'] =='seg':
        nSamples = summaryInfo['epochNumSamps'][0];
        nTrials = len(summaryInfo['epochNumSamps'])  
#        nTrials = summaryInfo['blocks']  
        # if Time0 is the same for all segments...
        if len(set(summaryInfo['epochTime0'])) == 1 :
            nSamplesPre = summaryInfo['epochTime0'][0] 
    else :
        nSamples = sum(summaryInfo['epochNumSamps'])
        nTrials = 1 

    ##---------------------------------------------------------------------------
    # Add the sensor info. 
    sensorLayoutfile = filePath+'/sensorLayout.xml'
    sensorLayoutObj = parse(sensorLayoutfile)    

    sensors = sensorLayoutObj.getElementsByTagName('sensor')
    label = []
    chantype = []
    chanunit = []
    tmpLabel = []
    nChans = 0
    for sensor in sensors:
        sensortype = int(sensor.getElementsByTagName('type')[0].firstChild.data)
        if sensortype == 0 or sensortype == 1 :
            if sensor.getElementsByTagName('name')[0].firstChild == None: 
                sn = sensor.getElementsByTagName('number')[0].firstChild.data.encode() 
                tmpLabel = 'E'+ sn.decode()
#                tmpLabel = 'E'+ sensor.getElementsByTagName('number')[0].firstChild.data.encode() 
            else:
                sn = sensor.getElementsByTagName('name')[0].firstChild.data.encode()
                tmpLabel = sn.decode()
            label.append(tmpLabel)
            chantype.append('eeg')
            chanunit.append('uV')        
            nChans = nChans + 1
    if nChans != summaryInfo['nChans']:
        print("Error. Should never occur.")

    ##----------------------------------------------------------------------------
    if summaryInfo['pibNChans'] > 0 :
        pnsSetfile = filePath + '/pnsSet.xml' 
        pnsSetObj = parse(pnsSetfile)    
        
        pnsSensors = pnsSetObj.getElementsByTagName('sensor')
        for p in range(summaryInfo['pibNChans']):
            tmpLabel = 'pib' + str(p+1)
            label.append(tmpLabel)
            pnsSensorObj = pnsSensors[p] 
            chantype.append(pnsSensorObj.getElementsByTagName('name')[0].firstChild.data.encode())
            chanunit.append(pnsSensorObj.getElementsByTagName('unit')[0].firstChild.data.encode())
           
    nChans = nChans + summaryInfo['pibNChans']
    
    ##-------------------------------------------------------------------------------
    header = {'Fs':summaryInfo['sampRate'], 'nChans':nChans,
       'nSamplesPre':nSamplesPre, 'nSamples':nSamples,'nTrials':nTrials,'label':label,'chantype':chantype, 'chanunit':chanunit, 'orig':summaryInfo }

    return header


 
