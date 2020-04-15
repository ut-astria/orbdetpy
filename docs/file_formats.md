orbdetpy uses JSON files to store settings and measurements for both
data simulation and orbit determination. The following sections describe
the formatting requirements for these files.

All strings in key names and values are case sensitive. Timestamps must be
in UTC and given in the ISO-8601 format "yyyy-MM-ddThh:mm:ss.ffffffZ".

Configuration Files
-------------------

Configuration files are needed for both simulation and orbit determination.

"rsoMass": Object mass [kg].
    
"rsoArea": Object average cross sectional area [m^2].

"gravityDegree": EGM96 gravity field degree; integer in [2, 360].

"gravityOrder": EGM96 gravity field order; integer in [0, degree].

"oceanTidesDegree": FES2004 ocean tide degree; -1 to disable or integer in [1, 100].

"oceanTidesOrder": FES2004 ocean tide order; -1 to disable or integer in [0, degree].

"thirdBodySun": Third body perturbations; true to enable the Sun's contribution.
 
"thirdBodyMoon": Third body perturbations; true to enable the Moon's contribution.

"solidTidesSun": Solid tides; true to enable the Sun's contribution.

"solidTidesMoon": Solid tides; true to enable the Moon's contribution.

"dragModel": "MSISE" for NRL MSISE-00 or "Exponential" for exponential drag.

"dragCoefficient": Drag coefficient {

  "value": Initial value.

  "min": Minimum value.

  "max": Maximum value.

  "estimation": "Estimate" to estimate, "Consider" to only consider the parameter in the UKF, or anything else to do neither.
  
}

"dragMSISEFlags": nx2 array of integers where n is in [1,23]. The first column is the flag number and the second column is the value of the corresponding MSISE flag. See <https://www.orekit.org/site-orekit-development/apidocs/org/orekit/forces/drag/atmosphere/NRLMSISE00.html> for details. All flags default to 1.

"dragExpRho0": Density constant [kg/m^3] for exponential drag.
 
"dragExpH0": Altitude constant [m] for exponential drag.
 
"dragExpHscale": Altitude scale factor [m] for exponential drag.

"rpSun": true to enable solar radiation pressure.
 
"rpCoeffReflection": Coefficient of reflection {

  "value": Initial value.
 
  "min": Minimum value.

  "max": Maximum value.

  "estimation": "Estimate" to estimate, "Consider" to only consider the parameter in the UKF, or anything else to do neither.

}

"cfgManeuvers": One or more maneuvers to include during simulation or less commonly during orbit determination [

 {
 
  "triggerEvent": "DateTime" to trigger at the specified time instant.
                   "LongitudeCrossing" to detect when the given geodetic longitude is crossed.

  "triggerParams": Not required for "DateTime".
  		    [geodetic_longitude_radians] for "LongitudeCrossing" detection.

  "time": Time string for maneuver; required only when TriggerEvent is "DateTime".

  "maneuverType": One of ["ConstantThrust", "NorthSouthStationing", "EastWestStationing", "SemiMajorAxisChange",
                   "PerigeeChange", "EccentricityChange", "InclinationChange", "RAANChange", "ArgPerigeeChange",
		   "StopPropagation"].

  "maneuverParams": [dir_x, dir_y, dir_z, duration_seconds, thrust_Newtons, Isp_seconds] for "ConstantThrust".
    		     [target_geodetic_latitude, deltav_interval, deltav_count] for "NorthSouthStationing".
                     [target_geodetic_longitude, deltav_interval, deltav_count] for "EastWestStationing".
		     [target_SMA, deltav_interval, deltav_count] for "SemiMajorAxisChange".
                     [target_perigee_radius, deltav_interval, deltav_count] for "PerigeeChange".
		     [target_eccentricity, deltav_interval, deltav_count] for "EccentricityChange".
		     [target_inclination, deltav_interval, deltav_count] for "InclinationChange".
		     [target_RAAN, deltav_interval, deltav_count] for "RAANChange".
		     [target_perigee_argument, deltav_interval, deltav_count] for "ArgPerigeeChange".
		     Not required for "StopPropagation".

 }
 
 ]

"propStart": Epoch for InitialState; optional if InitialTLE is provided.

"propEnd": End time for RSO state propagation.

"propStep": Time step size [s]; used only for simulation runs.

"propInitialState": Initial state vector. Units are [m] and [m/s] for position and velocity, respectively.

"propInitialTLE": Array of 2 TLE strings which may me provided instead of InitialState.

"propInertialFrame": Inertial reference frame to use for propagation and all state vector outputs. Valid
                     values are ["GCRF", "EME2000"]; defaults to "EME2000".

"integMinTimeStep": Minimum integration time step [s]; defaults to 1.0E-3 s.

"integMaxTimeStep": Maximum integration time step [s]; defaults to 300.0 s.

"integAbsTolerance": Integration absolute tolerance; defaults to 1.0E-14.

"integRelTolerance": Integration relative tolerance; defaults to 1.0E-12.

"simMeasurements": true to simulate measurements or false to only output state vectors;
                   defaults to true.

"simSkipUnobservable": true to skip instants when the object is below the local horizon; defaults to true.

"simIncludeExtras": true to generate additional simulated data; defaults to false.

"simIncludeStationState": true to include station inertial coordinates; defaults to false.

"cfgStations": Ground stations for measurements. {

 "Station_Name": {
 
  "latitude": Geodetic latitude [rad].
  
  "longitude": Geodetic longitude [rad].
  
  "altitude": Height above mean sea level [m].

  "azimuthBias": Optional sensor bias [rad].

  "elevationBias": Optional sensor bias [rad].

  "rangeBias": Optional sensor bias [m].

  "rangeRateBias": Optional sensor bias [m/s].

  "rightAscensionBias": Optional sensor bias [rad].

  "declinationBias": Optional sensor bias [rad].

  "positionBias": Optional sensor bias [m, m, m].

  "positionVelocityBias": Optional sensor bias [m, m, m, m/s, m/s, m/s].

  "biasEstimation": (UKF only) "Estimate" to estimate, "Consider" to only consider biases, or anything else to do neither.

  }
  
 }

"cfgMeasurements": Configure input measurements for orbit determination or output measurements from simulated data {

 "range": {

  "twoWay": true or false.

  "error": Theoretical measurement error [m].
  
 }

 "rangeRate": {

  "twoWay": true or false.

  "error": Theoretical measurement error [m/s].

 }

 "azimuth": {

  "error": Theoretical measurement error [rad].

 }

 "elevation": {

  "error": Theoretical measurement error [rad].

 }

 "rightAscension": {

  "error": Theoretical measurement error [rad].

 }

 "declination": {

  "error": Theoretical measurement error [rad].

 }

 "position": {

  "error": Theoretical measurement error [m, m, m].

 }

 "positionVelocity": {

  "error": Theoretical measurement error [m, m, m, m/s, m/s, m/s].

 }
 
 }

"estmFilter": Must be either "UKF" or "EKF".

"estmCovariance": Diagonal elements of covariance matrix with dimension 6 plus number of estimated parameters.

"estmProcessNoise": Diagonal elements of process noise matrix with dimension 6. Not used when DMC is in effect.

"estmDMCCorrTime": DMC correlation time. Setting this to zero disables DMC.

"estmDMCSigmaPert": Sigma for DMC acceleration. Setting this to zero disables DMC.

"estmDMCAcceleration": DMC acceleration bounds {
 
    "value": Initial value [m/s^2].
    
    "min": Minimum value [m/s^2].
    
    "max": Maximum value [m/s^2].
    
  }

"estmOutlierSigma": Number of sigmas from innovation for a measurement to be an outlier; positive number enables outlier detection.

"estmOutlierWarmup": Number of measurements to process before enabling outlier detection; positive number enables outlier detection.

Valid combinations of measurements are:

1. range
2. rangeRate
3. range + rangeRate
4. azimuth + elevation
5. rightAscension + declination
6. position
7. positionVelocity

Input Files
-----------

Only orbit determination requires input (measurement) files, which must
have the following structure. Each entry in the array corresponds to the
measurement(s) taken at a particular time instant and must conform to the
valid combinations listed  above.

[

 {
 
  "time": Measurement time stamp
  
  "station": Ground station name(s) from the configuration file's "cfgStations" array.
  
  "range": Optional based on measurements configured in "cfgMeasurements" [m].
  
  "rangeRate": Optional based on measurements configured in "cfgMeasurements" [m/s].

  "azimuth": Optional based on measurements configured in "cfgMeasurements" [rad].

  "elevation": Optional based on measurements configured in "cfgMeasurements" [rad].

  "rightAscension": Optional based on measurements configured in "cfgMeasurements" [rad].

  "declination": Optional based on measurements configured in "cfgMeasurements" [rad].

  "position": Optional based on measurements configured in "cfgMeasurements" [m, m, m].

  "positionVelocity": Optional based on measurements configured in "cfgMeasurements" [m, m, m, m/s, m/s, m/s].

 }

]
