/*
 * Settings.java - Functions to parse OD configuration settings.
 * Copyright (C) 2018-2022 University of Texas
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 */

package org.astria;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.geometry.spherical.twod.S2Point;
import org.hipparchus.geometry.spherical.twod.SphericalPolygonsSet;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.MatrixUtils;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.Attitude;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.attitudes.BodyCenterPointing;
import org.orekit.attitudes.FixedRate;
import org.orekit.attitudes.NadirPointing;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.forces.BoxAndSolarArraySpacecraft;
import org.orekit.forces.ForceModel;
import org.orekit.forces.drag.DragForce;
import org.orekit.forces.drag.DragSensitive;
import org.orekit.forces.drag.IsotropicDrag;
import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
import org.orekit.forces.gravity.NewtonianAttraction;
import org.orekit.forces.gravity.OceanTides;
import org.orekit.forces.gravity.SolidTides;
import org.orekit.forces.gravity.ThirdBodyAttraction;
import org.orekit.forces.gravity.potential.GravityFieldFactory;
import org.orekit.forces.gravity.potential.NormalizedSphericalHarmonicsProvider;
import org.orekit.forces.maneuvers.ConstantThrustManeuver;
import org.orekit.forces.radiation.IsotropicRadiationSingleCoefficient;
import org.orekit.forces.radiation.RadiationSensitive;
import org.orekit.forces.radiation.SolarRadiationPressure;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.LocalOrbitalFrame;
import org.orekit.frames.LOFType;
import org.orekit.frames.Predefined;
import org.orekit.frames.TopocentricFrame;
import org.orekit.geometry.fov.CircularFieldOfView;
import org.orekit.models.earth.EarthITU453AtmosphereRefraction;
import org.orekit.models.earth.atmosphere.Atmosphere;
import org.orekit.models.earth.atmosphere.NRLMSISE00;
import org.orekit.models.earth.atmosphere.SimpleExponentialAtmosphere;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.events.ApsideDetector;
import org.orekit.propagation.events.DateDetector;
import org.orekit.propagation.events.EclipseDetector;
import org.orekit.propagation.events.ElevationDetector;
import org.orekit.propagation.events.GeographicZoneDetector;
import org.orekit.propagation.events.GroundFieldOfViewDetector;
import org.orekit.propagation.events.LatitudeCrossingDetector;
import org.orekit.propagation.events.LatitudeExtremumDetector;
import org.orekit.propagation.events.LongitudeCrossingDetector;
import org.orekit.propagation.events.LongitudeExtremumDetector;
import org.orekit.propagation.events.NodeDetector;
import org.orekit.propagation.events.handlers.EventHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;

public final class Settings
{
    public static class Parameter
    {
	public String name;
	public double value;
	public double min;
	public double max;
	public EstimationType estimation;

	public Parameter(String name, double min, double max, double value, EstimationType estimation)
	{
	    this.name = name;
	    this.min = min;
	    this.max = max;
	    this.value = value;
	    this.estimation = estimation;
	}
    }

    public static class Facet
    {
	public double[] normal;
	public double area;

	public Facet(double[] normal, double area)
        {
	    this.normal = normal;
	    this.area = area;
	}
    }

    public static class Maneuver
    {
	public AbsoluteDate time;
	public ManeuverTrigger triggerEvent;
	public double[] triggerParams;
	public ManeuverType maneuverType;
	public double[] maneuverParams;

	public Maneuver(AbsoluteDate time, ManeuverTrigger triggerEvent, double[] triggerParams,
			ManeuverType maneuverType, double[] maneuverParams)
        {
	    this.time = time;
	    this.triggerEvent = triggerEvent;
	    this.triggerParams = triggerParams;
	    this.maneuverType = maneuverType;
	    this.maneuverParams = maneuverParams;
	}
    }

    public static class Station
    {
	public double latitude;
	public double longitude;
	public double altitude;
	public double[] bias;
	public EstimationType biasEstimation;
	public double fovAzimuth;
	public double fovElevation;
	public double fovAperture;

	public Station(double latitude, double longitude, double altitude, double[] bias, EstimationType biasEstimation,
		       double fovAzimuth, double fovElevation, double fovAperture)
        {
	    this.latitude = latitude;
	    this.longitude = longitude;
	    this.altitude = altitude;
	    this.bias = bias;
	    this.biasEstimation = biasEstimation;
	    this.fovAzimuth = fovAzimuth;
	    this.fovElevation = fovElevation;
	    this.fovAperture = fovAperture;
	}
    }

    public static class Measurement
    {
	public boolean twoWay;
	public double[] error;

	public Measurement(boolean twoWay, double[] error)
        {
	    this.twoWay = twoWay;
	    this.error = error;
	}
    }

    public static enum AttitudeType
    {
	UNDEFINED(0),
	NADIR_POINTING(1),
	BODY_CENTER_POINTING(2),
	FIXED_RATE(3);

	private final int value;

	AttitudeType(int value)
	{
	    this.value = value;
	}
    }

    public static enum DragModel
    {
	UNDEFINED(0),
	EXPONENTIAL(1),
	MSISE2000(2),
	WAM(3);

	private final int value;

	DragModel(int value)
	{
	    this.value = value;
	}
    }

    public static enum EstimationType
    {
	UNDEFINED(0),
	CONSIDER(1),
	ESTIMATE(2);

	private final int value;

	EstimationType(int value)
	{
	    this.value = value;
	}
    }

    public static enum Filter
    {
	EXTENDED_KALMAN(0),
	UNSCENTED_KALMAN(1),
	BATCH_LEAST_SQUARES(2);

	private final int value;

	Filter(int value)
	{
	    this.value = value;
	}
    }

    public static enum ManeuverTrigger
    {
	UNDEFINED(0),
	DATE_TIME(1),
	LONGITUDE_CROSSING(2),
	APSIDE_CROSSING(3),
	LONGITUDE_EXTREMUM(4),
	LATITUDE_CROSSING(5),
	LATITUDE_EXTREMUM(6),
	NODE_CROSSING(7);

	private final int value;

	ManeuverTrigger(int value)
	{
	    this.value = value;
	}
    }

    public static enum ManeuverType
    {
	UNDEFINED(0),
	CONSTANT_THRUST(1),
	NORTH_SOUTH_STATIONING(2),
	EAST_WEST_STATIONING(3),
	SEMI_MAJOR_AXIS_CHANGE(4),
	PERIGEE_CHANGE(5),
	ECCENTRICITY_CHANGE(6),
	INCLINATION_CHANGE(7),
	RAAN_CHANGE(8),
	ARG_PERIGEE_CHANGE(9),
	STOP_PROPAGATION(10);

	private final int value;

	ManeuverType(int value)
	{
	    this.value = value;
	}
    }

    public static enum MeasurementType
    {
	AZIMUTH(0),
	ELEVATION(1),
	RANGE(2),
	RANGE_RATE(3),
	RIGHT_ASCENSION(4),
	DECLINATION(5),
	POSITION(6),
	POSITION_VELOCITY(7);

	private final int value;

	MeasurementType(int value)
	{
	    this.value = value;
	}
    }

    public static final int OUTPUT_ESTM_COV = 1;
    public static final int OUTPUT_PROP_COV = 2;
    public static final int OUTPUT_INNO_COV = 4;
    public static final int OUTPUT_RESIDUALS = 8;
    public static final int OUTPUT_DENSITY = 16;
    public static final int OUTPUT_ECLIPSE = 32;

    public double rsoMass = 5.0;
    public double rsoArea = 0.1;
    public Facet[] rsoFacets;
    public double[] rsoSolarArrayAxis;
    public double rsoSolarArrayArea;
    public AttitudeType rsoAttitudeProvider = AttitudeType.UNDEFINED;
    public double[] rsoSpinVelocity;
    public double[] rsoSpinAcceleration;

    public int gravityDegree = 20;
    public int gravityOrder = 20;
    public int oceanTidesDegree = 20;
    public int oceanTidesOrder = 20;
    public boolean thirdBodySun = true;
    public boolean thirdBodyMoon = true;
    public boolean solidTidesSun = true;
    public boolean solidTidesMoon = true;

    public DragModel dragModel = DragModel.MSISE2000;
    public Parameter dragCoefficient = new Parameter("Cd", 1.0, 3.0, 2.0, EstimationType.ESTIMATE);
    public int[][] dragMSISEFlags;
    public double dragExpRho0 = 3.614E-13;
    public double dragExpH0 = 700000.0;
    public double dragExpHscale = 88667.0;

    public boolean rpSun = true;
    public Parameter rpCoeffReflection = new Parameter("Cr", 1.0, 2.0, 1.5, EstimationType.ESTIMATE);
    public double rpCoeffAbsorption;

    public Maneuver[] cfgManeuvers;

    public AbsoluteDate propStart;
    public AbsoluteDate propEnd;
    public double propStep;
    public double[] propInitialState;
    public String[] propInitialTLE;
    public Frame propInertialFrame = FramesFactory.getFrame(Predefined.GCRF);

    public double integMinTimeStep = 1.0E-3;
    public double integMaxTimeStep = 300.0;
    public double integAbsTolerance = 1.0E-14;
    public double integRelTolerance = 1.0E-12;

    public boolean simMeasurements = false;

    public Map<String, Station> cfgStations;
    public Map<MeasurementType, Measurement> cfgMeasurements;
    public double[] geoZoneLatLon;

    public Filter estmFilter = Filter.UNSCENTED_KALMAN;
    public double[] estmCovariance = {25E6, 25E6, 25E6, 1E2, 1E2, 1E2, 1.00, 0.25, 1E-6, 1E-6, 1E-6};
    public double[] estmProcessNoise;
    public double estmDMCCorrTime = 40.0;
    public double estmDMCSigmaPert = 5.0E-9;
    public Parameter estmDMCAcceleration = new Parameter("DMC", -1E-3, 1E-3, 0.0, EstimationType.ESTIMATE);
    public double estmOutlierSigma = 0.0;
    public int estmOutlierWarmup = 0;
    public int estmSmootherIterations = 3;
    public double estmDetectionProbability = 0.99;
    public double estmGatingProbability = 0.99;
    public double estmGatingThreshold = 5.0;
    public int outputFlags = OUTPUT_ESTM_COV | OUTPUT_PROP_COV | OUTPUT_INNO_COV | OUTPUT_RESIDUALS;

    protected Atmosphere atmModel;
    protected HashMap<String, GroundStation> stations;
    protected ArrayList<ForceModel> forces;
    protected ArrayList<Parameter> parameters;
    protected RealMatrix parameterMatrix;

    public Settings build()
    {
	loadGroundStations();
	loadForces();
	loadParameters();
	return(this);
    }

    private void loadGroundStations()
    {
	stations = new HashMap<String, GroundStation>();
	if (cfgStations != null)
	{
	    for (Map.Entry<String, Station> kv: cfgStations.entrySet())
	    {
		String k = kv.getKey();
		Station v = kv.getValue();
		GroundStation sta = new GroundStation(
		    new TopocentricFrame(DataManager.earthShape, new GeodeticPoint(v.latitude, v.longitude, v.altitude), k));
		sta.getPrimeMeridianOffsetDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
		sta.getPolarOffsetXDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
		sta.getPolarOffsetYDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
		stations.put(k, sta);
	    }
	}
    }

    private void loadForces()
    {
	forces = new ArrayList<ForceModel>();
	if (oceanTidesDegree >= 0 && oceanTidesOrder >= 0)
	    forces.add(new OceanTides(FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP), Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
				      Constants.EGM96_EARTH_MU, oceanTidesDegree, oceanTidesOrder, IERSConventions.IERS_2010,
				      TimeScalesFactory.getUT1(IERSConventions.IERS_2010, false)));

	NormalizedSphericalHarmonicsProvider harmonics = null;
	if (gravityDegree >= 2 && gravityOrder >= 0)
	    harmonics = GravityFieldFactory.getNormalizedProvider(gravityDegree, gravityOrder);
	if ((solidTidesSun || solidTidesMoon) && harmonics != null)
	    forces.add(new SolidTides(FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP), Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
				      Constants.EGM96_EARTH_MU, harmonics.getTideSystem(), IERSConventions.IERS_2010,
				      TimeScalesFactory.getUT1(IERSConventions.IERS_2010, false),
				      CelestialBodyFactory.getSun(), CelestialBodyFactory.getMoon()));

	if (thirdBodySun)
	    forces.add(new ThirdBodyAttraction(CelestialBodyFactory.getSun()));
	if (thirdBodyMoon)
	    forces.add(new ThirdBodyAttraction(CelestialBodyFactory.getMoon()));

	DragSensitive dragsc = null;
	RadiationSensitive radnsc = null;
	if (rsoFacets != null && rsoSolarArrayAxis != null)
	{
	    BoxAndSolarArraySpacecraft.Facet[] facets = new BoxAndSolarArraySpacecraft.Facet[rsoFacets.length];
	    for (int i = 0; i < rsoFacets.length; i++)
		facets[i] = new BoxAndSolarArraySpacecraft.Facet(new Vector3D(rsoFacets[i].normal), rsoFacets[i].area);
	    dragsc = new BoxAndSolarArraySpacecraft(facets, CelestialBodyFactory.getSun(), rsoSolarArrayArea, new Vector3D(rsoSolarArrayAxis),
						    dragCoefficient.value, rpCoeffAbsorption, rpCoeffReflection.value);
	    radnsc = new BoxAndSolarArraySpacecraft(facets, CelestialBodyFactory.getSun(), rsoSolarArrayArea, new Vector3D(rsoSolarArrayAxis),
						    dragCoefficient.value, rpCoeffAbsorption, rpCoeffReflection.value);
	}
	else
	{
	    dragsc = new IsotropicDrag(rsoArea, dragCoefficient.value);
	    radnsc = new IsotropicRadiationSingleCoefficient(rsoArea, rpCoeffReflection.value, rpCoeffReflection.min, rpCoeffReflection.max);
	}

	switch (dragModel)
	{
	case EXPONENTIAL:
	    atmModel = new SimpleExponentialAtmosphere(DataManager.earthShape, dragExpRho0, dragExpH0, dragExpHscale);
	    break;
	case MSISE2000:
	    atmModel = new NRLMSISE00(MSISEInputs.getInstance(), CelestialBodyFactory.getSun(), DataManager.earthShape);
	    if (dragMSISEFlags != null)
	    {
		for (int i = 0; i < dragMSISEFlags.length; i++)
		    atmModel = ((NRLMSISE00)atmModel).withSwitch(dragMSISEFlags[i][0], dragMSISEFlags[i][1]);
	    }
	    break;
	case WAM:
	    atmModel = WAM.getInstance();
	    break;
	}

	if (atmModel != null)
	    forces.add(new DragForce(atmModel, dragsc));
	if (rpSun)
	    forces.add(new SolarRadiationPressure(CelestialBodyFactory.getSun(), Constants.WGS84_EARTH_EQUATORIAL_RADIUS, radnsc));

	if (cfgManeuvers != null)
	{
	    for (Maneuver m: cfgManeuvers)
	    {
		if (m.triggerEvent == ManeuverTrigger.DATE_TIME && m.maneuverType == ManeuverType.CONSTANT_THRUST)
		    forces.add(new ConstantThrustManeuver(m.time, m.maneuverParams[3], m.maneuverParams[4], m.maneuverParams[5],
							  new Vector3D(m.maneuverParams[0], m.maneuverParams[1], m.maneuverParams[2])));
	    }
	}

	if (harmonics != null)
	    forces.add(new HolmesFeatherstoneAttractionModel(FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP), harmonics));
	if (gravityDegree >= 0 && gravityOrder >= 0)
	    forces.add(new NewtonianAttraction(Constants.EGM96_EARTH_MU));
    }

    private void loadParameters()
    {
	double bias;
	final int[] counts = {0, 0};
	final EstimationType[] ops = {EstimationType.ESTIMATE, EstimationType.CONSIDER};

	parameters = new ArrayList<Parameter>();
	for (int i = 0; i < ops.length; i++)
	{
	    if (dragCoefficient.estimation == ops[i])
	    {
		counts[i]++;
		parameters.add(new Parameter(DragSensitive.DRAG_COEFFICIENT, dragCoefficient.min, dragCoefficient.max, dragCoefficient.value, ops[i]));
	    }

	    if (rpCoeffReflection.estimation == ops[i])
	    {
		counts[i]++;
		parameters.add(new Parameter(RadiationSensitive.REFLECTION_COEFFICIENT, rpCoeffReflection.min,
					     rpCoeffReflection.max, rpCoeffReflection.value, ops[i]));
	    }

	    for (int j = 0; j < 3 && i == 0 && estmDMCCorrTime > 0.0 && estmDMCSigmaPert > 0.0; j++)
	    {
		counts[0]++;
		parameters.add(new Parameter(Estimation.DMC_ACC_ESTM[j], estmDMCAcceleration.min,
					     estmDMCAcceleration.max, estmDMCAcceleration.value, ops[0]));
	    }

	    if (cfgStations == null || cfgMeasurements == null || estmFilter != Filter.UNSCENTED_KALMAN)
		continue;

	    MeasurementType[] measNames = cfgMeasurements.keySet().toArray(new Settings.MeasurementType[0]);
	    Arrays.sort(measNames);
	    for (Map.Entry<String, Station> skv: cfgStations.entrySet())
	    {
		final Station sv = skv.getValue();
		if (sv.biasEstimation == ops[i] && sv.bias != null && sv.bias.length > 0)
		{
		    final String sk = skv.getKey();
		    for (MeasurementType mk: measNames)
		    {
			if (mk == MeasurementType.AZIMUTH || mk == MeasurementType.RANGE || mk == MeasurementType.RIGHT_ASCENSION)
			    bias = sv.bias[0];
			else
			    bias = sv.bias[sv.bias.length - 1];
			counts[i]++;
			final String name = new StringBuilder(sk).append(mk).toString();
			parameters.add(new Parameter(name, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, bias, ops[i]));
		    }
		}
	    }
	}

	parameterMatrix = MatrixUtils.createRealIdentityMatrix(parameters.size() + 6);
	if (counts[1] > 0)
	{
	    final RealMatrix zeros = MatrixUtils.createRealMatrix(counts[1], counts[1]);
	    parameterMatrix.setSubMatrix(zeros.getData(), counts[0] + 6, counts[0] + 6);
	}
    }

    public double[] getInitialState()
    {
	double[] state0 = propInitialState;
	if (state0 == null)
	{
	    TLE parser = new TLE(propInitialTLE[0], propInitialTLE[1]);
	    if (propStart == null)
		propStart = parser.getDate();
	    PVCoordinates pvc = TLEPropagator.selectExtrapolator(parser).getPVCoordinates(propStart, propInertialFrame);
	    Vector3D p = pvc.getPosition();
	    Vector3D v = pvc.getVelocity();
	    state0 = new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()};
	}

	double[] X0 = new double[parameters.size() + 6];
	for (int i = 0; i < X0.length; i++)
	{
	    if (i < 6)
		X0[i] = state0[i];
	    else
		X0[i] = parameters.get(i - 6).value;
	}
	return(X0);
    }

    public RealMatrix getInitialCovariance()
    {
	int states = parameters.size() + 6;
	RealMatrix cov = MatrixUtils.createRealIdentityMatrix(states);
	for (int i = 0, k = 0; i < states; i++)
	{
	    if (2*estmCovariance.length == states*(states + 1))
	    {
		// Initialize with the given lower triangular entries
		for (int j = 0; j <= i; j++, k++)
		{
		    cov.setEntry(i, j, estmCovariance[k]); 
		    cov.setEntry(j, i, estmCovariance[k]); 
		}
	    }
	    else
	    {
		// Initialize with the given diagonal entries
		if (i < estmCovariance.length)
		    cov.setEntry(i, i, estmCovariance[i]); 
	    }
	}
	return(cov);
    }

    public AttitudeProvider getAttitudeProvider()
    {
	AttitudeProvider att = null;
	switch (rsoAttitudeProvider)
	{
	case NADIR_POINTING:
	    att = new NadirPointing(propInertialFrame, DataManager.earthShape);
	    break;
	case BODY_CENTER_POINTING:
	    att = new BodyCenterPointing(propInertialFrame, DataManager.earthShape);
	    break;
	case FIXED_RATE:
	    double[] X0 = propInitialState;
	    KeplerianPropagator prop = new KeplerianPropagator(new CartesianOrbit(new PVCoordinates(new Vector3D(X0[0], X0[1], X0[2]),
												    new Vector3D(X0[3], X0[4], X0[5])),
										  propInertialFrame, propStart, Constants.EGM96_EARTH_MU));
	    LocalOrbitalFrame lof = new LocalOrbitalFrame(propInertialFrame, LOFType.VVLH, prop, "");
	    att = new FixedRate(new org.orekit.attitudes.Attitude(propStart, lof, Rotation.IDENTITY,
								  new Vector3D(rsoSpinVelocity[0], rsoSpinVelocity[1], rsoSpinVelocity[2]),
								  new Vector3D(rsoSpinAcceleration[0], rsoSpinAcceleration[1], rsoSpinAcceleration[2])));
	    break;
	}
	return(att);
    }

    public RealMatrix getProcessNoiseMatrix(double t)
    {
	int i;
	t = FastMath.abs(t);
	final double t2 = t*t;
	final double t3 = t2*t;
	final double t4 = t3*t;
	final double[][] Q = new double[parameters.size() + 6][parameters.size() + 6];

	if (estmProcessNoise != null && estmProcessNoise.length == 6)
	{
	    final double[] P = estmProcessNoise;
	    for (i = 0; i < 3; i++)
	    {
		Q[i][i] = 0.25*t4*P[i];
		Q[i][i + 3] = 0.5*t3*P[i];
	    }

	    for (i = 3; i < 6; i++)
	    {
		Q[i][i] = t2*P[i];
		Q[i][i - 3] = 0.5*t3*P[i];
	    }
	    return(new Array2DRowRealMatrix(Q));
	}

	if (estmDMCCorrTime < 1E-6)
	    throw(new IllegalArgumentException());
	final int N = parameters.size() - 3;
	final double b = 1.0/estmDMCCorrTime;
	final double b2 = b*b;
	final double b3 = b2*b;
	final double b4 = b3*b;
	final double b5 = b4*b;
	final double et = FastMath.exp(-b*t);
	final double e2t = et*et;
	final double s2 = estmDMCSigmaPert*estmDMCSigmaPert;
	final double Q00 = s2*(t3/(3*b2) - t2/b3 + t*(1 - 2*et)/b4 + 0.5*(1 - e2t)/b5); // pos-pos
	final double Q01 = s2*(0.5*t2/b2 - t*(1 - et)/b3 + (1 - et)/b4 - 0.5*(1 - e2t)/b4); // pos-vel
	final double Q02 = s2*(0.5*(1 - e2t)/b3 - t*et/b2); // pos-acc
	final double Q11 = s2*(t/b2 - 2*(1 - et)/b3 + 0.5*(1 - e2t)/b3); // vel-vel
	final double Q12 = s2*(0.5*(1 + e2t)/b2 - et/b2); // vel-acc
	final double Q22 = 0.5*s2*(1 - e2t)/b; // acc-acc

	for (i = 0; i < 3; i++)
	{
	    Q[i][i] = Q00;
	    Q[i][i + 3] = Q01;
	    Q[i][i + N + 6] = Q02;
	}

	for (i = 3; i < 6; i++)
	{
	    Q[i][i] = Q11;
	    Q[i][i - 3] = Q01;
	    Q[i][i + N + 3] = Q12;
	}

	for (i = N + 6; i < N + 9; i++)
	{
	    Q[i][i] = Q22;
	    Q[i][i - N - 6] = Q02;
	    Q[i][i - N - 3] = Q12;
	}
	return(new Array2DRowRealMatrix(Q));
    }

    public ArrayList<EventHandling> addEventHandlers(Propagator prop, SpacecraftState initialState)
    {
	ArrayList<EventHandling> handles = new ArrayList<EventHandling>();
	if (geoZoneLatLon != null && geoZoneLatLon.length >= 6)
	{
	    S2Point[] vertices = new S2Point[(int)(geoZoneLatLon.length/2)];
	    for (int i = 0; i <= geoZoneLatLon.length - 2; i += 2)
		vertices[(int)(i/2)] = new S2Point(geoZoneLatLon[i+1], 0.5*FastMath.PI - geoZoneLatLon[i]);
	    SphericalPolygonsSet set = new SphericalPolygonsSet(1E-10, vertices);
	    GeographicZoneDetector detector = new GeographicZoneDetector(DataManager.earthShape, set, FastMath.toRadians(0.5));
	    EventHandler<GeographicZoneDetector> handler = new EventHandling<GeographicZoneDetector>(
		ManeuverType.UNDEFINED, 0, EventHandling.GEO_ZONE_NAME, detector.g(initialState) < 0.0);
	    prop.addEventDetector(detector.withHandler(handler));
	    handles.add((EventHandling)handler);
	}

	if ((outputFlags & Settings.OUTPUT_ECLIPSE) != 0)
	{
	    EclipseDetector detector = new EclipseDetector(CelestialBodyFactory.getSun(), 695700E3, DataManager.earthShape);
	    EventHandler<EclipseDetector> handler = new EventHandling<EclipseDetector>(
		ManeuverType.UNDEFINED, 0, EventHandling.UMBRA, detector.g(initialState) < 0.0);
	    prop.addEventDetector(detector.withHandler(handler));
	    handles.add((EventHandling)handler);
	    detector = detector.withPenumbra();
	    handler = new EventHandling<EclipseDetector>(ManeuverType.UNDEFINED, 0, EventHandling.PENUMBRA, detector.g(initialState) < 0.0);
	    prop.addEventDetector(detector.withHandler(handler));
	    handles.add((EventHandling)handler);
	}

	if (cfgStations != null)
	{
	    for (Map.Entry<String, Station> kv: cfgStations.entrySet())
	    {
		Station stn = kv.getValue();
		if (stn.fovAperture <= 1E-6)
		{
		    ElevationDetector detector = new ElevationDetector(stations.get(kv.getKey()).getBaseFrame())
			.withRefraction(new EarthITU453AtmosphereRefraction(stn.altitude)).withConstantElevation(FastMath.toRadians(5.0));
		    EventHandler<ElevationDetector> handler = new EventHandling<ElevationDetector>(
			ManeuverType.UNDEFINED, 0, kv.getKey(), detector.g(initialState) > 0.0);
		    prop.addEventDetector(detector.withHandler(handler));
		    handles.add((EventHandling)handler);
		}
		else
		{
		    Vector3D center = new Vector3D(FastMath.cos(stn.fovElevation)*FastMath.sin(stn.fovAzimuth),
						   FastMath.cos(stn.fovElevation)*FastMath.cos(stn.fovAzimuth), FastMath.sin(stn.fovElevation));
		    CircularFieldOfView fov = new CircularFieldOfView(center, 0.5*stn.fovAperture, 1E-6);
		    GroundFieldOfViewDetector detector = new GroundFieldOfViewDetector(stations.get(kv.getKey()).getBaseFrame(), fov);
		    EventHandler<GroundFieldOfViewDetector> handler = new EventHandling<GroundFieldOfViewDetector>(
			ManeuverType.UNDEFINED, 0, kv.getKey(), detector.g(initialState) < 0.0);
		    prop.addEventDetector(detector.withHandler(handler));
		    handles.add((EventHandling)handler);
		}
	    }
	}

	if (cfgManeuvers == null)
	    return(handles);
	for (Maneuver m: cfgManeuvers)
	{
	    if (m.maneuverType == ManeuverType.CONSTANT_THRUST)
		continue;
	    double delta = m.maneuverParams.length > 0 ? m.maneuverParams[0] : 0.0;
	    switch (m.triggerEvent)
	    {
	    case DATE_TIME:
		EventHandler<DateDetector> time = new EventHandling<DateDetector>(m.maneuverType, delta, null, false);
		prop.addEventDetector(new DateDetector(m.time).withHandler(time));
		handles.add((EventHandling)time);
		break;
	    case LONGITUDE_CROSSING:
		EventHandler<LongitudeCrossingDetector> lonc = new EventHandling<LongitudeCrossingDetector>(m.maneuverType, delta, null, false);
		prop.addEventDetector(new LongitudeCrossingDetector(DataManager.earthShape, m.triggerParams[0]).withHandler(lonc));
		handles.add((EventHandling)lonc);
		break;
	    case APSIDE_CROSSING:
		EventHandler<ApsideDetector> apse = new EventHandling<ApsideDetector>(m.maneuverType, delta, null, false);
		prop.addEventDetector(new ApsideDetector(initialState.getOrbit()).withHandler(apse));
		handles.add((EventHandling)apse);
		break;
	    case LONGITUDE_EXTREMUM:
		EventHandler<LongitudeExtremumDetector> lone = new EventHandling<LongitudeExtremumDetector>(m.maneuverType, delta, null, false);
		prop.addEventDetector(new LongitudeExtremumDetector(DataManager.earthShape).withHandler(lone));
		handles.add((EventHandling)lone);
		break;
	    case LATITUDE_CROSSING:
		EventHandler<LatitudeCrossingDetector> latc = new EventHandling<LatitudeCrossingDetector>(m.maneuverType, delta, null, false);
		prop.addEventDetector(new LatitudeCrossingDetector(DataManager.earthShape, m.triggerParams[0]).withHandler(latc));
		handles.add((EventHandling)latc);
		break;
	    case LATITUDE_EXTREMUM:
		EventHandler<LatitudeExtremumDetector> late = new EventHandling<LatitudeExtremumDetector>(m.maneuverType, delta, null, false);
		prop.addEventDetector(new LatitudeExtremumDetector(DataManager.earthShape).withHandler(late));
		handles.add((EventHandling)late);
		break;
	    case NODE_CROSSING:
		String type = (m.triggerParams.length == 0 || m.triggerParams[0] == 0.0) ? EventHandling.ASC_NODE : EventHandling.DES_NODE;
		EventHandler<NodeDetector> node = new EventHandling<NodeDetector>(m.maneuverType, delta, type, false);
		prop.addEventDetector(new NodeDetector(initialState.getOrbit(), propInertialFrame).withHandler(node));
		handles.add((EventHandling)node);
		break;
	    default:
		throw(new RuntimeException("Invalid maneuver trigger"));
	    }
	}
	return(handles);
    }
}
