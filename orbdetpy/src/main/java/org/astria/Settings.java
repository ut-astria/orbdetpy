/*
 * Settings.java - Functions to parse OD configuration settings.
 * Copyright (C) 2018-2019 University of Texas
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
import java.util.Map;
import java.util.HashMap;
import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.linear.RealMatrix;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.Attitude;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.attitudes.BodyCenterPointing;
import org.orekit.attitudes.FixedRate;
import org.orekit.attitudes.NadirPointing;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
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
import org.orekit.forces.maneuvers.ImpulseManeuver;
import org.orekit.forces.radiation.IsotropicRadiationSingleCoefficient;
import org.orekit.forces.radiation.RadiationSensitive;
import org.orekit.forces.radiation.SolarRadiationPressure;
import org.orekit.frames.Frame;
import org.orekit.frames.LocalOrbitalFrame;
import org.orekit.frames.LOFType;
import org.orekit.frames.TopocentricFrame;
import org.orekit.frames.Transform;
import org.orekit.models.earth.atmosphere.Atmosphere;
import org.orekit.models.earth.atmosphere.NRLMSISE00;
import org.orekit.models.earth.atmosphere.SimpleExponentialAtmosphere;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.events.DateDetector;
import org.orekit.propagation.events.LongitudeCrossingDetector;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.UT1Scale;
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
	public String estimation;

	public Parameter()
	{
	}

	public Parameter(String name, double min, double max, double value, String estimation)
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

	public Facet()
        {
	}

	public Facet(double[] normal, double area)
        {
	    this.normal = normal;
	    this.area = area;
	}
    }

    public static class Maneuver
    {
	public String time;
	public String triggerEvent;
	public double[] triggerParams;
	public String maneuverType;
	public double[] maneuverParams;

	public Maneuver()
        {
	}

	public Maneuver(String time, String triggerEvent, double[] triggerParams, String maneuverType, double[] maneuverParams)
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
	public double azimuthBias;
	public double elevationBias;
	public double rangeBias;
	public double rangeRateBias;
	public double rightAscensionBias;
	public double declinationBias;
	public double[] positionBias;
	public double[] positionVelocityBias;
	public String biasEstimation;

	public Station()
        {
	}

	public Station(double latitude, double longitude, double altitude, double azimuthBias, double elevationBias,
		       double rangeBias, double rangeRateBias, double rightAscensionBias, double declinationBias,
		       double[] positionBias, double[] positionVelocityBias, String biasEstimation)
        {
	    this.latitude = latitude;
	    this.longitude = longitude;
	    this.altitude = altitude;
	    this.azimuthBias = azimuthBias;
	    this.elevationBias = elevationBias;
	    this.rangeBias = rangeBias;
	    this.rangeRateBias = rangeRateBias;
	    this.rightAscensionBias = rightAscensionBias;
	    this.declinationBias = declinationBias;
	    this.positionBias = positionBias;
	    this.positionVelocityBias = positionVelocityBias;
	    this.biasEstimation = biasEstimation;
	}
    }

    public static class Measurement
    {
	public boolean twoWay;
	public double[] error;

	public Measurement()
        {
	}

	public Measurement(boolean twoWay, double[] error)
        {
	    this.twoWay = twoWay;
	    this.error = error;
	}
    }

    public double rsoMass = 1.0;
    public double rsoArea = 1.0;
    public Facet[] rsoFacets;
    public double[] rsoSolarArrayAxis;
    public double rsoSolarArrayArea;
    public String rsoAttitudeProvider;
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

    public String dragModel = "MSISE";
    public Parameter dragCoefficient = new Parameter("Cd", 1.0, 3.0, 2.0, "Estimate");
    public int[][] dragMSISEFlags;
    public double dragExpRho0;
    public double dragExpH0;
    public double dragExpHscale;

    public boolean rpSun = true;
    public Parameter rpCoeffReflection = new Parameter("Cr", 1.0, 2.0, 1.5, "Estimate");
    public double rpCoeffAbsorption;

    public Maneuver[] cfgManeuvers;

    public String propStart;
    public String propEnd;
    public double propStep;
    public double[] propInitialState;
    public String[] propInitialTLE;
    public String propInertialFrame = "EME2000";
    public String propStepHandlerStartTime;
    public String propStepHandlerEndTime;

    public double integMinTimeStep = 1.0E-3;
    public double integMaxTimeStep = 300.0;
    public double integAbsTolerance = 1.0E-14;
    public double integRelTolerance = 1.0E-12;

    public boolean simMeasurements = true;
    public boolean simSkipUnobservable = true;
    public boolean simIncludeExtras = false;
    public boolean simIncludeStationState = false;

    public Map<String, Station> cfgStations;
    public Map<String, Measurement> cfgMeasurements;

    public String estmFilter = "UKF";
    public double[] estmCovariance = new double[]{25E6, 25E6, 25E6, 1E-2, 1E-2, 1E-2, 1.00, 0.25, 1E-6, 1E-6, 1E-6};
    public double[] estmProcessNoise;
    public double estmDMCCorrTime = 40.0;
    public double estmDMCSigmaPert = 5.0E-9;
    public Parameter estmDMCAcceleration = new Parameter("DMC", -1E-3, 1E-3, 0.0, "Estimate");

    protected Atmosphere atmModel;
    protected HashMap<String, GroundStation> stations;
    protected ArrayList<ForceModel> forces;
    protected ArrayList<Parameter> parameters;
    protected Frame propFrame;

    public Settings build()
    {
	propFrame = DataManager.getFrame(propInertialFrame);
	loadGroundStations();
	loadForces();
	loadParameters();
	return(this);
    }

    private void loadGroundStations()
    {
	stations = new HashMap<String, GroundStation>();
	if (cfgStations == null)
	    return;

	for (Map.Entry<String, Station> kv : cfgStations.entrySet())
	{
	    String k = kv.getKey();
	    Station v = kv.getValue();
	    GroundStation sta = new GroundStation(new TopocentricFrame(new OneAxisEllipsoid(
									   Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
									   Constants.WGS84_EARTH_FLATTENING, DataManager.getFrame("ITRF")),
								       new GeodeticPoint(v.latitude, v.longitude, v.altitude), k));
	    sta.getPrimeMeridianOffsetDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
	    sta.getPolarOffsetXDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
	    sta.getPolarOffsetYDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
    	    stations.put(k, sta);
	}
    }

    private void loadForces()
    {
	forces = new ArrayList<ForceModel>();
	NormalizedSphericalHarmonicsProvider grav = null;
	if (gravityDegree >= 2 && gravityOrder >= 0)
	{
	    grav = GravityFieldFactory.getNormalizedProvider(gravityDegree, gravityOrder);
	    forces.add(new HolmesFeatherstoneAttractionModel(DataManager.getFrame("ITRF"), grav));
	}
	else
	    forces.add(new NewtonianAttraction(Constants.EGM96_EARTH_MU));

	if (oceanTidesDegree >= 0 && oceanTidesOrder >= 0)
	    forces.add(new OceanTides(DataManager.getFrame("ITRF"), Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.EGM96_EARTH_MU,
				      oceanTidesDegree,	oceanTidesOrder, IERSConventions.IERS_2010,
				      (UT1Scale) DataManager.getTimeScale("UT1")));

	if ((solidTidesSun || solidTidesMoon) && grav != null)
	    forces.add(new SolidTides(DataManager.getFrame("ITRF"), Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.EGM96_EARTH_MU,
				      grav.getTideSystem(), IERSConventions.IERS_2010, (UT1Scale) DataManager.getTimeScale("UT1"),
				      CelestialBodyFactory.getSun(), CelestialBodyFactory.getMoon()));

	if (thirdBodySun)
	    forces.add(new ThirdBodyAttraction(CelestialBodyFactory.getSun()));
	if (thirdBodyMoon)
	    forces.add(new ThirdBodyAttraction(CelestialBodyFactory.getMoon()));

	DragSensitive dragsc = null;
	RadiationSensitive radnsc = null;
	if (rsoFacets != null && rsoSolarArrayAxis != null && rsoSolarArrayArea > 0.0)
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

	if (dragModel.equals("Exponential"))
	{
	    atmModel = new SimpleExponentialAtmosphere(new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
									    Constants.WGS84_EARTH_FLATTENING, DataManager.getFrame("ITRF")),
						       dragExpRho0, dragExpH0, dragExpHscale);
	}
	else if (dragModel.equals("MSISE"))
	{
	    int apflag = 1;
	    if (dragMSISEFlags != null)
	    {
		for (int i = 0; i < dragMSISEFlags.length; i++)
		{
		    if (dragMSISEFlags[i][0] == 9)
			apflag = dragMSISEFlags[i][1];
		}
	    }

	    atmModel = new NRLMSISE00(new MSISEInputs(DataManager.msiseData.minDate, DataManager.msiseData.maxDate,
						      DataManager.msiseData.data, apflag),
				      CelestialBodyFactory.getSun(), new OneAxisEllipsoid(
					  Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
					  Constants.WGS84_EARTH_FLATTENING, DataManager.getFrame("ITRF")));
	    if (dragMSISEFlags != null)
	    {
		for (int i = 0; i < dragMSISEFlags.length; i++)
		    atmModel = ((NRLMSISE00) atmModel).withSwitch(dragMSISEFlags[i][0], dragMSISEFlags[i][1]);
	    }
	}

	if (atmModel != null)
	    forces.add(new DragForce(atmModel, dragsc));

	if (rpSun)
	    forces.add(new SolarRadiationPressure(CelestialBodyFactory.getSun(), Constants.WGS84_EARTH_EQUATORIAL_RADIUS, radnsc));

	if (cfgManeuvers == null)
	    return;
	for (Maneuver m : cfgManeuvers)
	{
	    if (m.triggerEvent.equals("DateTime") && m.maneuverType.equals("ConstantThrust"))
		forces.add(new ConstantThrustManeuver(DataManager.parseDateTime(m.time), m.maneuverParams[3], m.maneuverParams[4], m.maneuverParams[5],
						      new Vector3D(m.maneuverParams[0], m.maneuverParams[1], m.maneuverParams[2])));
	}
    }

    private void loadParameters()
    {
	int[] counts = new int[]{0, 0};
	String[] ops = {"Estimate", "Consider"};
	parameters = new ArrayList<Parameter>();
	for (int i = 0; i < ops.length; i++)
	{
	    if (dragCoefficient.estimation != null && dragCoefficient.estimation.equals(ops[i]))
	    {
		counts[i]++;
		parameters.add(new Parameter(DragSensitive.DRAG_COEFFICIENT, dragCoefficient.min,
					     dragCoefficient.max, dragCoefficient.value, ops[i]));
	    }

	    if (rpCoeffReflection.estimation != null && rpCoeffReflection.estimation.equals(ops[i]))
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

	    if (cfgStations == null || cfgMeasurements == null || !estmFilter.equals("UKF"))
		continue;

	    for (Map.Entry<String, Station> skv : cfgStations.entrySet())
	    {
		Station sv = skv.getValue();
		if (sv.biasEstimation != null && sv.biasEstimation.equals(ops[i]))
		{
		    String sk = skv.getKey();
		    for (String mk : cfgMeasurements.keySet())
		    {
			double bias = 0.0;
			if (mk.equals("Azimuth"))
			    bias = sv.azimuthBias;
			else if (mk.equals("Elevation"))
			    bias = sv.elevationBias;
			else if (mk.equals("Range"))
			    bias = sv.rangeBias;
			else if (mk.equals("RangeRate"))
			    bias = sv.rangeRateBias;
			else if (mk.equals("RightAscension"))
			    bias = sv.rightAscensionBias;
			else if (mk.equals("Declination"))
			    bias = sv.declinationBias;
			else
			    continue;

			counts[i]++;
			String name = new StringBuilder(sk).append(mk).toString();
			parameters.add(new Parameter(name, Double.NEGATIVE_INFINITY, Double.POSITIVE_INFINITY, bias, ops[i]));
		    }
		}
	    }
	}
    }

    public double[] getInitialState()
    {
	PVCoordinates topv;
	double[] state0 = propInitialState;

	if (state0 == null)
	{
	    AbsoluteDate epoch;
	    TLE parser = new TLE(propInitialTLE[0], propInitialTLE[1]);
	    TLEPropagator prop = TLEPropagator.selectExtrapolator(parser);

	    if (propStart != null)
		epoch = DataManager.parseDateTime(propStart);
	    else
	    {
		epoch = parser.getDate().shiftedBy(0.0);
		propStart = DataManager.getUTCString(epoch);
	    }
	    topv = prop.getPVCoordinates(epoch, propFrame);
	}
	else
	    topv = new PVCoordinates(new Vector3D(state0[0], state0[1], state0[2]),
				     new Vector3D(state0[3], state0[4], state0[5]));

	Vector3D p = topv.getPosition();
	Vector3D v = topv.getVelocity();
	state0 = new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()};
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

    public AttitudeProvider getAttitudeProvider()
    {
	if (rsoAttitudeProvider == null)
	    return(null);

	AttitudeProvider attpro = null;
	OneAxisEllipsoid shape = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
						      Constants.WGS84_EARTH_FLATTENING, DataManager.getFrame("ITRF"));

	if (rsoAttitudeProvider.equals("NadirPointing"))
	    attpro = new NadirPointing(propFrame, shape);
	if (rsoAttitudeProvider.equals("BodyCenterPointing"))
	    attpro = new BodyCenterPointing(propFrame, shape);
	if (rsoAttitudeProvider.equals("FixedRate") && rsoSpinVelocity != null && rsoSpinAcceleration != null)
	{
	    double[] X0 = propInitialState;
	    AbsoluteDate t0 = DataManager.parseDateTime(propStart);
	    KeplerianPropagator prop = new KeplerianPropagator(new CartesianOrbit(new PVCoordinates(new Vector3D(X0[0], X0[1], X0[2]),
												    new Vector3D(X0[3], X0[4], X0[5])),
										  propFrame, t0, Constants.EGM96_EARTH_MU));
	    LocalOrbitalFrame lof = new LocalOrbitalFrame(propFrame, LOFType.VVLH, prop, "");
	    attpro = new FixedRate(new org.orekit.attitudes.Attitude(t0, lof, Rotation.IDENTITY,
								     new Vector3D(rsoSpinVelocity[0], rsoSpinVelocity[1], rsoSpinVelocity[2]),
								     new Vector3D(rsoSpinAcceleration[0], rsoSpinAcceleration[1], rsoSpinAcceleration[2])));
	}

	return(attpro);
    }

    public RealMatrix getProcessNoiseMatrix(double t)
    {
	int i;
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
		Q[i][i+3] = 0.5*t3*P[i];
	    }

	    for (i = 3; i < 6; i++)
	    {
		Q[i][i] = t2*P[i];
		Q[i][i-3] = 0.5*t3*P[i];
	    }

	    return(new Array2DRowRealMatrix(Q));
	}

	final int N = parameters.size() - 3;
	final double b = 1.0/estmDMCCorrTime;
	final double b2 = b*b;
	final double b3 = b2*b;
	final double b4 = b3*b;
	final double b5 = b4*b;
	final double et = FastMath.exp(-1.0*b*t);
	final double e2t = et*et;
	final double s2 = estmDMCSigmaPert*estmDMCSigmaPert;

	final double Q00 = s2*(t3/(3*b2)-t2/b3+t*(1-2*et)/b4+0.5*(1-e2t)/b5); // pos-pos
	final double Q01 = s2*(0.5*t2/b2-t*(1-et)/b3+(1-et)/b4-0.5*(1-e2t)/b4); // pos-vel
	final double Q02 = s2*(0.5*(1-e2t)/b3-t*et/b2); // pos-acc
	final double Q11 = s2*(t/b2-2*(1-et)/b3+0.5*(1-e2t)/b3); // vel-vel
	final double Q12 = s2*(0.5*(1+e2t)/b2-et/b2); // vel-acc
	final double Q22 = 0.5*s2*(1-e2t)/b; // acc-acc

	for (i = 0; i < 3; i++)
	{
	    Q[i][i] = Q00;
	    Q[i][i+3] = Q01;
	    Q[i][i+N+6] = Q02;
	}

	for (i = 3; i < 6; i++)
	{
	    Q[i][i] = Q11;
	    Q[i][i-3] = Q01;
	    Q[i][i+N+3] = Q12;
	}

	for (i = N+6; i < N+9; i++)
	{
	    Q[i][i] = Q22;
	    Q[i][i-N-6] = Q02;
	    Q[i][i-N-3] = Q12;
	}

	return(new Array2DRowRealMatrix(Q));
    }

    public void addEventHandlers(NumericalPropagator prop, SpacecraftState state)
    {
	if (cfgManeuvers == null)
	    return;
	for (Maneuver m : cfgManeuvers)
	{
	    if (m.triggerEvent.equals("DateTime"))
	    {
		if (m.maneuverType.equals("ConstantThrust"))
		    continue;

		EventHandling<DateDetector> handler = new EventHandling<DateDetector>(m.triggerEvent, m.maneuverType,
										      m.maneuverParams[0], (int) m.maneuverParams[2]);
		AbsoluteDate time = DataManager.parseDateTime(m.time);
		for (int i = 0; i < m.maneuverParams[2]; i++)
		{
		    prop.addEventDetector(new DateDetector(time).withHandler(handler));
		    time = new AbsoluteDate(time, m.maneuverParams[1]);
		}
	    }

	    if (m.triggerEvent.equals("LongitudeCrossing"))
	    {
		OneAxisEllipsoid body = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
							     Constants.WGS84_EARTH_FLATTENING, DataManager.getFrame("ITRF"));
		EventHandling<LongitudeCrossingDetector> handler = new EventHandling<LongitudeCrossingDetector>(
		    m.triggerEvent, m.maneuverType, 0.0, 1);
		prop.addEventDetector(new LongitudeCrossingDetector(body, m.triggerParams[0]).withHandler(handler));
	    }
	}
    }
}
