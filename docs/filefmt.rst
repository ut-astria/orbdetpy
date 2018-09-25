============================
orbdetpy - JSON file formats
============================

orbdetpy uses JSON files to store settings and measurements for both
data simulation and orbit determination. The following sections describe
the formatting requirements for these files.

Note that all strings in JSON files are case sensitive. Timestamps
must be in UTC and given by the format "yyyy-MM-ddThh:mm:ss.fffZ".

Configuration Files
-------------------

Configuration files are needed for both simulation and orbit determination.

"Gravity" : Configure EGM96 gravity field {

 "Degree" : Integer in [2,360].

 "Order" : Integer in [0,degree].

 }

"OceanTides" : Configure FES2004 ocean tide model {

 "Degree" : -1 to disable ocean tides or integer in [1,100].

 "Order" : -1 to disable ocean tides or integer in [0,degree].

 }

"Drag" : Configure atmospheric density model {

 "Model" : "MSISE" for NRL MSISE-00 or "Exponential" for exponential drag.

 "MSISEDisable" : Array of integers in [1,23] to disable various MSISE-00 settings. Click `here <https://www.orekit.org/site-orekit-development/apidocs/org/orekit/forces/drag/atmosphere/NRLMSISE00.html>`_ for details.

 "ExpRho0" : Density constant [kg/m^3] for exponential drag.
 
 "ExpH0": Altitude constant [m] for exponential drag.
 
 "ExpHScale" : Altitude scale factor [m] for exponential drag.

 "Coefficient" : Drag coefficient {
 
    "Value" : Initial value.
    
    "Min" : Minimum value.
    
    "Max" : Maximum value.
    
    "Estimation" : "Estimate" to estimate parameter or any other string to prevent estimation.
    
    }
    
 }

"SolidTides" : Configure solid tides {

 "Sun" : false to disable the Sun's contribution or true to enable.

 "Moon" : false to disable the Moon's contribution or true to enable.

 }

"ThirdBodies" : Configure third body point mass perturbations {

 "Sun" : false to disable the Sun's contribution or true to enable.
 
 "Moon" : false to disable the Moon's contribution or true to enable.

 }

"RadiationPressure" : Configure radiation pressure {

 "Sun" : false to disable the Sun's contribution or true to enable.
 
 "Creflection" : Coefficient of reflection {

  "Value" : Initial value.
 
  "Min" : Minimum value.

  "Max" : Maximum value.

  "Estimation" : "Estimate" to estimate parameter or any other string to prevent estimation.
  
  }

 }

"SpaceObject" : Details about the Resident Space Object {

 "Mass" : Object mass [kg].
    
 "Area" : Object average cross sectional area [m^2].

 }

"Propagation" : Propagation settings {

 "Start" : Epoch for InitialState.

 "End" : End time for RSO state propagation.

 "InitialState" : Initial state vector in J2000. Must be of size six + number of estimated parameters. Units are [m] and [m/s] for first six elements.

 }

"Stations" : One or more ground stations for measurement processing {

 "Station1" : {
 
  "Latitude" : Geodetic latitude [rad].
  
  "Longitude" : Geodetic longitude [rad].
  
  "Altitude" : Height above Mean Sea Level [m].
  
  }
  
 }

"Maneuvers" : One or more constant thrust maneuvers to include during simulation or less commonly with orbit determination [

 {
  "Time" : Time of maneuver.

  "Duration" : Maneuver duration [s].

  "Thrust" : Thrust force [N].

  "Isp" : Engine specific impulse [s].

  "Direction" : Unit vector in the RSO frame specifying thrust direction.
  
 }
 
 ]

"Measurements" : Configure input measurements for orbit determination or output measurements from simulated data {

 "Range" : {

  "TwoWay" : true or false.

  "Enabled" : true or false. Functionality is not currently implemented.

  "Error" : Theoretical measurement error [m].
  
 }

 "RangeRate" : {

  "TwoWay" : true or false.

  "Enabled" : true or false. Functionality is not currently implemented.

  "Error" : Theoretical measurement error [m/s].

 }

 "Azimuth" : {

  "Enabled" : true or false. Functionality is not currently implemented.

  "Error" : Theoretical measurement error [rad].

 }

 "Elevation" : {

  "Enabled" : true or false. Functionality is not currently implemented.

  "Error" : Theoretical measurement error [rad].

 }

 "RightAscension" : {

  "Enabled" : true or false. Functionality is not currently implemented.

  "Error" : Theoretical measurement error [rad].

 }

 "Declination" : {

  "Enabled" : true or false. Functionality is not currently implemented.

  "Error" : Theoretical measurement error [rad].

 }

 }

Valid combinations of measurements are as follows:

1) Range
2) RangeRate
3) Range + RangeRate
4) Azimuth + Elevation
5) RightAscension + Declination
 
"Estimation" : Configure parameters for estimation filters {

 "Covariance" : Diagonal elements of covariance matrix with the same dimensions as InitialState.

 "ProcessNoise" : Diagonal elements of process noise matrix with the same dimensions as InitialState.

 }

Input Files
-----------

Only orbit determination requires input (measurement) files, which must
have the following structure. Each entry in the array corresponds to the
measurement(s) taken at a particular time instant and must conform to the
valid combinations listed  above.

[

 {
 
  "Time" : Measurement time stamp
  
  "Station" : Ground station name(s) from the configuration file's "Stations" array.
  
  "Range" : Optional based on measurements configured in "Measurements" [m].
  
  "RangeRate" : Optional based on measurements configured in "Measurements" [m/s].

  "Azimuth" : Optional based on measurements configured in "Measurements" [rad].

  "Elevation" : Optional based on measurements configured in "Measurements" [rad].

  "RightAscension" : Optional based on measurements configured in "Measurements" [rad].

  "Declination" : Optional based on measurements configured in "Measurements" [rad].

 }

]
