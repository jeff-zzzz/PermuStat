# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 15:16:22 2015
@author: jesong1126
"""
>>> dir()
>>> print dir()
['CatchOutputs', 'GeoPy', 'Helper', 'MFF', 'PyEGI', '__builtins__', '__doc__', '__name__', '__package__', 'catchOutError', 'egiData', 'os', 'sys']
>>> print ss
<module 'libEGIPythonAdapter' from '/Applications/EAV/Mimir.app/Contents/MacOS/libEGIPythonAdapter.so'>
>>> print dir(ss)
['Condition', 'CorticalMesh', 'EEGNet', 'LeadFieldMatrix', 'PythonAdapter', 'ScalarArray', '__doc__', '__file__', '__name__', '__package__']
>>> PyEGI = egiData.PythonAdapter()
>>> print(PyEGI)
<libEGIPythonAdapter.PythonAdapter object at 0x113c6d788>
>>> __doc__
>>> print(__doc__)
None
>>> print GeoPy
<module 'EGIPython.GeoPy' from '/Applications/EAV/Mimir.app/Contents/MacOS/EGIPython/GeoPy/__init__.pyc'>
>>> PyEGI.addOrientedCondition("LFM")
>>> c = PyEGI.getCondition("LFM")
>>> nE = 256
>>> for i in range(nE):
...         c.addSensorData(K[i,j])
>>> help(egiData)
Help on module libEGIPythonAdapter:

NAME
    libEGIPythonAdapter

FILE
    /Applications/EAV/Mimir.app/Contents/MacOS/libEGIPythonAdapter.so

CLASSES
    Boost.Python.instance(__builtin__.object)
        Condition
        CorticalMesh
        EEGNet
        LeadFieldMatrix
        PythonAdapter
        ScalarArray
    
    class Condition(Boost.Python.instance)
     |  Condition():
     |  The Condition represens a set of values representing a condition in EEG.  Each condition may have dipole data (either oriented or triples) as well as sensor data.  These data are then used in the Condition View to display source localization results and EEG data at each sensor glyph.
     |  
     |  Method resolution order:
     |      Condition
     |      Boost.Python.instance
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(...)
     |  
     |  __reduce__ = <unnamed Boost.Python function>(...)
     |  
     |  addSensorData(...)
     |      void addSensorData(float value):
     |      Add data to the sensors.  This will add the value to the back of the vector.  This means that data should be added to the Condition sequentially, beginning from sensor 0 and ending with the last sensor.
     |  
     |  clearSensorData(...)
     |      void clearSensorData():
     |      Clears the sensor data associated with the Condition.  This will set the value returned by "numSensors()" to 0.
     |  
     |  getDipoleDirectionX(...)
     |      float getDipoleDirectionX(int dipoleID):
     |      Get the X component of the dipole orientation associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipoleDirectionY(...)
     |      float getDipoleDirectionX(int dipoleID):
     |      Get the Y component of the dipole orientation associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipoleDirectionZ(...)
     |      float getDipoleDirectionX(int dipoleID):
     |      Get the Z component of the dipole orientation associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipolePositionX(...)
     |      float getDipolePositionX(int dipoleID):
     |      Get the X component of the dipole position associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipolePositionY(...)
     |      float getDipolePositionX(int dipoleID):
     |      Get the Y component of the dipole position associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipolePositionZ(...)
     |      float getDipolePositionX(int dipoleID):
     |      Get the Z component of the dipole position associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipoleScalar(...)
     |      float getDipoleScalar(int dipoleID):
     |      Get the scalar value associated with the dipole at dipoleID.  Note:  Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDisplayTime(...)
     |      string getDisplayTime():
     |      Get the string representing the time for the data represented by the Condition.
     |  
     |  getName(...)
     |      string getName():
     |      Get the name of the Condition.  Conditions are indexed based on their name, so care must be taken when defining new conditions.  The convention used is to prepend the string "T: " for conditions containing triples data, "O: " for conditions containing oriented dipole solutions, and "E: " for conditions containing only EEG data.
     |  
     |  getSensorData(...)
     |      float getSensorData(int sensorID):
     |      Get the data associated with sensor sensorID.  If sensorID cannot be found, 0.0 is returned.  Note: Sensors are indexed from 0, rather than 1.
     |  
     |  clearSensorData(...)
     |      void clearSensorData():
     |      Clears the sensor data associated with the Condition.  This will set the value returned by "numSensors()" to 0.
     |  
     |  getDipoleDirectionX(...)
     |      float getDipoleDirectionX(int dipoleID):
     |      Get the X component of the dipole orientation associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipoleDirectionY(...)
     |      float getDipoleDirectionX(int dipoleID):
     |      Get the Y component of the dipole orientation associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipoleDirectionZ(...)
     |      float getDipoleDirectionX(int dipoleID):
     |      Get the Z component of the dipole orientation associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipolePositionX(...)
     |      float getDipolePositionX(int dipoleID):
     |      Get the X component of the dipole position associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipolePositionY(...)
     |      float getDipolePositionX(int dipoleID):
     |      Get the Y component of the dipole position associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipolePositionZ(...)
     |      float getDipolePositionX(int dipoleID):
     |      Get the Z component of the dipole position associated with dipole dipoleID.  Note: Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDipoleScalar(...)
     |      float getDipoleScalar(int dipoleID):
     |      Get the scalar value associated with the dipole at dipoleID.  Note:  Dipoles are indexed beginning from 0 rather than 1.
     |  
     |  getDisplayTime(...)
     |      string getDisplayTime():
     |      Get the string representing the time for the data represented by the Condition.
     |  
     |  getName(...)
     |      string getName():
     |      Get the name of the Condition.  Conditions are indexed based on their name, so care must be taken when defining new conditions.  The convention used is to prepend the string "T: " for conditions containing triples data, "O: " for conditions containing oriented dipole solutions, and "E: " for conditions containing only EEG data.
     |  
     |  getSensorData(...)
     |      float getSensorData(int sensorID):
     |      Get the data associated with sensor sensorID.  If sensorID cannot be found, 0.0 is returned.  Note: Sensors are indexed from 0, rather than 1.
     |  
     |  isEEGOnly(...)
     |      bool isEEGOnly():
     |      Return whether or not the Condition contains both EEG data as w
     |  ----------------------------------------------------------------------
     |  Data and other attributes inherited from Boost.Python.instance:
     |  
     |  __new__ = <built-in method __new__ of Boost.Python.class object>
     |      T.__new__(S, ...) -> a new object with type S, a subtype of T
    
    class CorticalMesh(Boost.Python.instance)
     |  CorticalMesh():
     |  The CorticalMesh represents either a standard or inflated version of the cortical mesh used in the head model.
     |  Construction of CorticalMesh objects should NOT be done programatically.  Instead, the objects should be
     |  first loaded into the environment through either the PyEGI interface or the File drop-down menu.
     |  Once loaded, the PyEGI interface should be used to retrieve the appropriate objects for manipulation.
     |  
     |  Method resolution order:
     |      CorticalMesh
     |      Boost.Python.instance
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(...)
     |  
     |  __reduce__ = <unnamed Boost.Python function>(...)
     |  
     |  addDipoleToSelection(...)
     |      void addDipoleToSelection(int dipole):
     |      Add the given dipole to the currently active selected set.
     |  
     |  clearSelection(...)
     |      void clearSelection():
     |      Clear any selected dipoles that have been identified by user driven interaction.
     |  
     |  disableOrientations(...)
     |      void disableOrientations():
     |      Disable the display of the orientation vectors for selected dipoles on the cortical surface.
     |  
     |  drawOrientations(...)
     |      bool drawOrientations():
     |      Return a boolean based on whether orientation vectors should be displayed for selected dipoles.
     |  
     |  enableOrientations(...)
     |      void enableOrientations():
     |      Enable the display of the orientation vectors for selected dipoles on the cortical surface.
     |  
     |  getDipoleClosestToPoint(...)
     |      int getDipoleClosestToPoint(float x, float y, float z):
     |      Return the dipole closest to the 3D point (x,y,z).  If a suitable dipole cannot be found, -1 is returned.
     |  
     |  getDipoleScalar(...)
     |      float getDipoleScalar(int dipole):
     |      Return the scalar value at a given dipole using the currently displayed scalar field.
     |  
     |  getDipoleScalarByIdAndColor(...)
     |      float getDipoleScalarByIdAndColor(int dipole, ColorMethod method, bool normalize):
     |      Return the scalar value at a given dipole.  The scalar array the value is drawn from is defined by the ColorMethod
     |      passed as an argument.  If normalize is set to True, the value return will be clamped to the range [0,1].
     |      The ColorMethod variable is an enumeration defined as follows:
     |              Total Energy = 0
     |              Radial Energy = 1
     |              EEG Net Color = 2
     |              Global Color = 3
     |              Dipole Color = 4
     |              External Data (Cortical) = 5
     |              External Data (EEG Net) = 6
     |              Charge Density = 7
     |  
     |  setScalarArray(...)
     |      void setScalarArray(ScalarArray):
     |      Set the scalar array for the External Data coloring.  The ScalarArray object must have at least as many entries as the
     |      highest dipole ID resident in the CorticalMesh object.
     |      Note:  This method reproduces functionality defined in the PythonAdapter class.
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __instance_size__ = 784
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Boost.Python.instance:
     |  
     |  __dict__
     |  
     |  __weakref__
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes inherited from Boost.Python.instance:
     |  
     |  __new__ = <built-in method __new__ of Boost.Python.class object>
     |      T.__new__(S, ...) -> a new object with type S, a subtype of T
    
    class EEGNet(Boost.Python.instance)
     |  EEGNet():
     |  The EEGNet represents the EEG Net used during acquisition and then transformed using the GPS system.
     |  Construction of EEGNet objects should NOT be done programatically.  Instead, the objects should be
     |  first loaded into the environment through the PyEGI interface for the File drop-down menu.
     |  Once loaded, the PyEGI interface should be used to retrieve the appropriate data object for manipulation.
     |  
     |  Method resolution order:
     |      EEGNet
     |      Boost.Python.instance
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(...)
     |  
     |  __reduce__ = <unnamed Boost.Python function>(...)
     |  
     |  clearSelections(...)
     |      void clearSelections():
     |      Clear any selected EEG Electrodes that have been identified by the user.
     |  
     |  getLeadFieldMatrix(...)
     |      LeadFieldMatrix getLeadfieldMatrix():
     |      Get the LeadfieldMatrix associated with this EEGNet.
     |  
     |  getNearestElectrode(...)
     |      int getNearestElectrode(float x, float y, float z):
     |      Return the ID of the electrode nearest to the point (x,y,z).
     |      If no appropriate electrode is found, -1 is returned.
     |  
     |  selectElectrode(...)
     |      void selectElectrode(int electrodeID, ElectrodeType type):
     |      Select the electrode with ID electrodeID and assign it the given ElectrodeType.
     |      The ElectrodeType variable is an enumeration defined as follows:
     |              Cathode Electrode = -1
     |              Neutral Electrode = 0
     |              Anode Electrode = 1
     |      Note:  The electrode type is taken into account only when the GTEN environment is enabled.
     |  
     |  setScalarArray(...)
     |      void setScalarArray(ScalarArray scalars):
     |      Set the scalar array for the External Data coloring.  The ScalarArray object must have at least as many entries as the
     |      highest EEG Electrode ID resident in the EEGNet object.
     |      Note:  This method reproduces functionality defined in the PythonAdapter class.
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __instance_size__ = 368
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Boost.Python.instance:
     |  
     |  __dict__
     |  
     |  __weakref__
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes inherited from Boost.Python.instance:
     |  
     |  __new__ = <built-in method __new__ of Boost.Python.class object>
     |      T.__new__(S, ...) -> a new object with type S, a subtype of T
    
    class LeadFieldMatrix(Boost.Python.instance)
     |  LeadFieldMatrix():
     |  The LeadfieldMatrix represents the Leadfield computed in the Physics module representing the forward projection of each dipole activity from the cortical surface to the scalp surface.
     |  Construction of the LeadfieldMatrix should NOT be done programatically.  Instead, the object should be first loaded into the environment by the EAV->loadHeadModel command.  Once loaded, the PyEGI interface should be used to retrieve the matrix for access in Python.
     |  
     |  Method resolution order:
     |      LeadFieldMatrix
     |      Boost.Python.instance
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(...)
     |  
     |  __reduce__ = <unnamed Boost.Python function>(...)
     |  
     |  allocateInverse(...)
     |      void allocateInverse():
     |      Allocate the memory for the inverse of the current matrix.  This inverse will be the shape of the transposed leadfield.
     |  
     |  getCols(...)
     |      int getCols():
     |      Get the number of columns in the matrix.
     |  
     |  getDataElement(...)
     |      float getData(int dipole, int electrode, int dimension):
     |      Get the matrix entry associated with a dipole and electrode.  The dimension input parameter should be 0 in the case of cortical surface solutions, and [0,1,2] in the case of triples.
     |  
     |  getNumDipoles(...)
     |      int getNumDipoles():
     |      Get the number of dipoles (columns) represented in the matrix.
     |  
     |  getNumElectrodes(...)
     |      int getNumElectrodes():
     |      Get the number of electrodes (rows) represented in the matrix.
     |  
          |  getRows(...)
     |      int getRows():
     |      Get the number of rows in the matrix.
     |  
     |  setInverseElement(...)
     |      void setInverseElement(int row, int col, float data):
     |      Set a single entry in the inverse matrix.
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __instance_size__ = 72
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Boost.Python.instance:
     |  
     |  __dict__
     |  
     |  __weakref__
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes inherited from Boost.Python.instance:
     |  
     |  __new__ = <built-in method __new__ of Boost.Python.class object>
     |      T.__new__(S, ...) -> a new object with type S, a subtype of T
    
    class PythonAdapter(Boost.Python.instance)
     |  This interface defines general methods for loading data from files, retrieving existing data objects, and assigning data to data objects.
     |  
     |  Method resolution order:
     |      PythonAdapter
     |      Boost.Python.instance
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(...)
     |  
     |  __reduce__ = <unnamed Boost.Python function>(...)
     |  
     |  addComboBoxOption(...)
     |      void addComboBoxOption(int widgetID, string boxLabel, string option):  Add an option to an existing QComboBox identified by its label inside a Basic Python Widget given its widgetID.  If the QComboBox cannot be found (via the boxLabel), a new QLabel:QComboBox pair will be created and the option will be added to it.
     |  
     |  addComboBoxWithLabel(...)
     |      void addComboBoxWithLabel(int widgetID, string label):  Add a QLabel:QComboBox pair to the Basic Python Widget given by the widgetID.
     |  
     |  addDoubleSpinBoxWithLabel(...)
     |      void addDoubleSpinBoxWithLabel(int widgetID, string label, double minVal, double maxVal, double step, double initVal):  Add a QDoubleSpinBox to the Basic Python Widget with id widgetID.  The QDoubleSpinBox will have the range [minVal, maxVal] and the step size will be step.  Its value will be initialized to initVal
     |  
     |  addIntegerSpinBoxWithLabel(...)
     |      void addIntegerSpinBoxWithLabel(int widgetID, string label, int minVal, int maxVal, int step, int initVal):  Add a QSpinBox to the Basic Python Widget with id widgetID.  The QSpinBox will have the range [minVal, maxVal] and the step size will be step.  Its value will be initialized to initVal
     |  
     |  addLineEditWithLabel(...)
     |      void addLineEditWithLabel(int widgetID, string editLabel, string initVal):  Add a QLabel:QLineEdit pair tot he Basic Python Widget given by the widgetID.  If supplied, the QLineEdit will be populated with the string in initVal.
     |  
     |  addOrientedCondition(...)
     |      void addOrientedCondition(Condition cond):  Add a Condition, cond, to the list of available oriented dipole conditions.  If a condition with cond.getName() alread exists, no change is made.
     |  
     |  addPushButton(...)
     |      void addPushButton(int widgetID, string buttonName):  Add a QPushButton with the given name.  For this button to react to the triggered() signal, a call to setCallback_onButtonPRessed(...) must be made to specify the Python handler for the signal.
     |  
     |  addTriplesCondition(...)
     |      void addTriplesCondition(Condition cond):  Add a Condition, cond, to the list of available triples dipole conditions.  If a condition with cond.getName() alread exists, no change is made.
     |  
     |  clearBasicPythonWidgets(...)
     |      void clearBasicPythonWidgets(int widgetID):  Clear the children widgets inside a Basic Python Widget.  This uses the widgetID displayed on the top of each widget in the Custom View Widget.
     |  
     |  getComboBoxSelection(...)
     |      string getComboBoxSelection(int widgetID, string boxLabel):  Get the currently active selection from the QComboBox with label boxLabel in the Basic Python Widget widgetID.
     |  
     |  getCondition(...)
     |      Condition getCondition(string name):  Retrieve the Condition data object given its name.
     |  
     |  getCorticalMesh(...)
     |      CorticalMesh getCorticalMesh(meshID):  Retrieve the CorticalMesh data object given its meshID.
     |  
     |  getDoubleSpinBoxValue(...)
     |      int getDoubleSpinBoxValue(int widgetID, string spinBoxLabel):  Return the current integer value from the QDoubleSpinBox specified by spinBoxLabel in the Basic Python Widget widgetID.
     |  
     |  getEEGNet(...)
     |      EEGNet getEEGNet(int netID):  Retrieve the EEGNet data object given its netID.
     |  
     |  getIntegerSpinBoxValue(...)
     |      int getIntegerSpinBoxValue(int widgetID, string spinBoxLabel):  Return the current integer value from the QSpinBox specified by spinBoxLabel in the Basic Python Widget widgetID.
     |  
     |  getLineEditValue(...)
     |      void getLineEditValue(int widgetID, string lineEditLabel):  Get the current string from the QLineEdit with label lineEditLabel inside the Basic Python Widget given by widgetID.
     |  
     |  loadCorticalMesh(...)
     |      void loadCorticalMesh(string filename):  Load a cortical mesh from a file with the supplied file name.
     |  
     |  loadEEGNet(...)
     |      void loadEEGNet(string filename):  Load an EEGNet from a file with the supplied file name.
     |  
     |  loadLeadfield(...)
     |      void loadLeadfieldMatrhx(string filename, int netID):  Load a Leadfield Matrix from a file with the supplied file name
     |      and attach it to the EEGNet defined by netID.
     |  
     |  loadMRI(...)
     |      void loadMRI(string filename):  Load an MR Volume Image from a file with the supplied file name.
     |  
     |  loadScalpMesh(...)
     |      void loadScalpMesh(string filename):  Load a scalp mesh from a file with the supplied file name.
     |  
     |  loadTransferMatrix(...)
     |      void loadTransferMatrix(string filename, int meshID):  Load a Transfer Matrix from a file with the supplied file name
     |      and attach it to the cortical mesh defined by meshID.
     |      Note:  This function is only available when GTEN functionality is enabled.
     |  
     |  notifyWindows_OrientedConditionDataChanged(...)
     |      void notifyWindows_OrientedConditionDataChanged(string conditionName):  Notify the display windows that the data associted with the oriented condition, conditionName, has changed.  This will most likely trigger an update of the visualization.
     |  
     |  setCallback_onButtonPressed(...)
     |      void setCallback_onButtonPressed(string buttonName, method func):  Set the function to be called when the button with name = buttonName is pressed (the signal, triggered() is sent).  This callback (specified in Python) will then be run.
     |      Example:
     |      PyEgi.setCallback_onButtonPressed(MyFunction)
     |  
     |  setCallback_onEEGDataChanged(...)
     |      void setCallback_onEEGDataChanged(method func):  Set the function to be called when the application experiences an EEGDataChanged signal.  This callback (specified in Python) will then be called whenever the EEG data is changed (typically in the NetStation Review interface.
     |      Example:
     |      PyEGI.setCallback_onEEGDataChanged(GeoPy.IMatrix.printSelf)
     |  
     |  setCorticalMeshScalars(...)
     |      void setCorticalMeshScalars(int meshID, ScalarArray scalars):  Set the scalar array to display as External Data for a mesh with the supplied meshID.
     |  
     |  setEEGNetScalars(...)
     |      void setEEGNetScalars(int netID, ScalarArray scalars):  Set the scalar array to display as External Data for an EEGNet with the supplied netID.
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __instance_size__ = 56
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Boost.Python.instance:
     |  
     |  __dict__
     |  
     |  __weakref__
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes inherited from Boost.Python.instance:
     |  
     |  __new__ = <built-in method __new__ of Boost.Python.class object>
     |      T.__new__(S, ...) -> a new object with type S, a subtype of T
    
    class ScalarArray(Boost.Python.instance)
     |  ScalarArray(string name, int numValues, float initialValue):
     |  The ScalarArray represents an array of values suitable for mapping onto a CorticalMesh,
     |  ScalpMesh, or EEGNet data object.  Each ScalarArray must have a name, a size, and a
     |  value to initialize the array with.  In order to assign a ScalarArray to a data object, it must have sufficient capacity.
     |  ScalarArray objects smaller than expected will not be accepted for assignment to the data object.
     |  
     |  Method resolution order:
     |      ScalarArray
     |      Boost.Python.instance
     |      __builtin__.object
     |  
     |  Methods defined here:
     |  
     |  __init__(...)
     |  
     |  __reduce__ = <unnamed Boost.Python function>(...)
     |  
     |  clear(...)
     |      void clear():  Clear the contents of the ScalarArray object.  This does not resize the array.
     |  
     |  copy(...)
     |      void copy(ScalarArray other):  Copy the contents of the other array.  This performs a deep copy of the data.
     |  
     |  getData(...)
     |      float getData(int loc):  Return the value of the ScalarArray at the given location, loc.
     |  
     |  setData(...)
     |      void setData(int loc, float val):  Assign the value, val, to the array location, loc.
     |  
     |  size(...)
     |      int size():  Return the current size of the ScalarArray object.
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes defined here:
     |  
     |  __instance_size__ = 80
     |  
     |  ----------------------------------------------------------------------
     |  Data descriptors inherited from Boost.Python.instance:
     |  
     |  __dict__
     |  
     |  __weakref__
     |  
     |  ----------------------------------------------------------------------
     |  Data and other attributes inherited from Boost.Python.instance:
     |  
     |  __new__ = <built-in method __new__ of Boost.Python.class object>
     |      T.__new__(S, ...) -> a new object with type S, a subtype of T



