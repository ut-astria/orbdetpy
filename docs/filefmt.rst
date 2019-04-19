============================
orbdetpy - JSON file formats
============================

orbdetpy uses JSON files to store settings and measurements for both
data simulation and orbit determination. The following sections describe
the formatting requirements for these files.

All strings in key names and values are case sensitive. Timestamps must be
in UTC and given using the ISO 8601 format "yyyy-MM-ddThh:mm:ss.fffZ".

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

 "MSISEFlags" : nx2 array of integers where n is in [1,23]. The first column is the flag number and the second column is the value of the corresponding MSISE flag. See `here <https://www.orekit.org/site-orekit-development/apidocs/org/orekit/forces/drag/atmosphere/NRLMSISE00.html>`_ for details. All flags default to 1.

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

 "Start" : Epoch for InitialState; optional if InitialTLE is provided.

 "End" : End time for RSO state propagation.

 "Step" : Time step size [s]; used only for simulation runs.

 "InitialState" : Initial state vector. Units are [m] and [m/s] for position and velocity, respectively.

 "InitialTLE" : Array of 2 TLE strings which may me provided instead of InitialState.

 "InitialStateFrame" : Reference frame for InitialState. Valid values are ["GCRF", "EME2000", "ITRF"];
                       defaults to "EME2000".
 
 "InertialFrame" : Inertial reference frame to use for propagation and all state vector outputs. Valid
                   values are ["GCRF", "EME2000"]; defaults to "EME2000".

 }

"Simulation" : Optional simulation settings; not used for orbit determination runs {

 "SimulateMeasurements" : true to simulate measurements or false to only output state vectors;
                          defaults to true.

 "SkipUnobservable" : true to skip instants when the object is below the local horizon; defaults to true.

 "IncludeExtras" : true to generate additional simulated data; defaults to false.

 }

"Integration" : Optional tolerances for numerical integration {

 "MinTimeStep" : Minimum integration time step [s]; defaults to 1.0E-3 s.

 "MaxTimeStep" : Maximum integration time step [s]; defaults to 300.0 s.

 "AbsTolerance" : Integration absolute tolerance; defaults to 1.0E-14.

 "RelTolerance" : Integration relative tolerance; defaults to 1.0E-12.

}
 
"Stations" : Ground stations for measurements. Not required for "PositionVelocity" measurements {

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

  "Error" : Theoretical measurement error [m].
  
 }

 "RangeRate" : {

  "TwoWay" : true or false.

  "Error" : Theoretical measurement error [m/s].

 }

 "Azimuth" : {

  "Error" : Theoretical measurement error [rad].

 }

 "Elevation" : {

  "Error" : Theoretical measurement error [rad].

 }

 "RightAscension" : {

  "Error" : Theoretical measurement error [rad].

 }

 "Declination" : {

  "Error" : Theoretical measurement error [rad].

 }

 "PositionVelocity" : {

  "Error" : Theoretical measurement error [m, m, m, m/s, m/s, m/s].

  "ReferenceFrame" : Reference frame in which position/velocity vectors are expressed.
                     Valid values are ["GCRF", "EME2000", "ITRF"]; defaults to "EME2000".  

 }
 
 }

Valid combinations of measurements are as follows:

1) Range
2) RangeRate
3) Range + RangeRate
4) Azimuth + Elevation
5) RightAscension + Declination
6) PositionVelocity
 
"Estimation" : Configure parameters for estimation filters {

 "Filter" : Must be either "UKF" or "EKF".

 "Covariance" : Diagonal elements of covariance matrix with dimension 6 plus number of estimated parameters.

 "ProcessNoise" : Diagonal elements of process noise matrix with dimension 6. Not used when DMC is in effect.

 "NoiseTimeDelta" : Delta-T to use for computing the SNC and DMC process noise matrices.

 "DMCCorrTime" : DMC correlation time. Setting this to zero disables DMC.

 "DMCSigmaPert" : Sigma for DMC acceleration. Setting this to zero disables DMC.

 "DMCAcceleration" : DMC acceleration bounds {
 
    "Value" : Initial value [m/s^2].
    
    "Min" : Minimum value [m/s^2].
    
    "Max" : Maximum value [m/s^2].
    
    }

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

  "PositionVelocity" : Optional based on measurements configured in "Measurements" [m, m, m, m/s, m/s, m/s].

 }

]
