##%% mff_getObject.m
#  Matlab File
#  author Colin Davey
#  date 7/22/2014
#  Copyright 2014 EGI. All rights reserved.
#  Support routine for MFF Matlab code. Not intended to be called directly.
#  
# Returns objects of various types: 
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_Info = 7
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_EventTrack = 3
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_Categories = 9
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_Signal = 2
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_Epochs = 4
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_InfoN = 8
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_SensorLayout = 10
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_Coordinates = 11
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_History = 6
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_MFFFile = 1
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_Photogrammetry = 12
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_Subject = 5
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_Any = 0
# com.egi.services.mff.api.MFFResourceType.kMFF_RT_Unknown = -1
##

def mff_getObject(objType, filename, path):
    URI = path 
    if objType != com.egi.services.mff.api.MFFResourceType.kMFF_RT_MFFFile:
        URI = [URI filesep filename]

    delegate =
    javaObject('com.egi.services.mff.api.LocalMFFFactoryDelegate') 
    factory = javaObject('com.egi.services.mff.api.MFFFactory',
    delegate) 
    resourceVal = objType 
    resourceType =
    javaObject('com.egi.services.mff.api.MFFResourceType',
    resourceVal) 
    # fprintf('%s %s\n', char(URI), char(resourceType));
    theObject = factory.openResourceAtURI(URI, resourceType) 
    if ~isempty(theObject): 
        theObject.loadResource();

    return theObject


