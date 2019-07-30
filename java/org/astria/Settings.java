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

import com.google.gson.Gson;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.stream.Stream;
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
import org.orekit.forces.drag.atmosphere.Atmosphere;
import org.orekit.forces.drag.atmosphere.NRLMSISE00;
import org.orekit.forces.drag.atmosphere.SimpleExponentialAtmosphere;
import org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel;
import org.orekit.forces.gravity.NewtonianAttraction;
import org.orekit.forces.gravity.OceanTides;
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
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.KeplerianPropagator;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.events.DateDetector;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;

public class Settings
{
    class JSONParameter
    {
	double Value;
	double Min;
	double Max;
	String Estimation;
    }

    class JSONGravity
    {
	int Degree;
	int Order;
    }

    class JSONOceanTides
    {
	int Degree;
	int Order;
    }

    class JSONDrag
    {
	String Model;
	int[][] MSISEFlags;
	double ExpRho0;
	double ExpH0;
	double ExpHScale;
	JSONParameter Coefficient;
    }

    class JSONSolidTides
    {
	boolean Sun;
	boolean Moon;
    }

    class JSONThirdBodies
    {
	boolean Sun;
	boolean Moon;
    }

    class JSONRadiationPressure
    {
	boolean Sun;
	JSONParameter Creflection;
	Double Cabsorption;
    }

    class JSONFacet
    {
	double[] Normal;
	double Area;
    }

    class JSONSolarArray
    {
	double[] Axis;
	double Area;
    }

    class JSONAttitude
    {
	String Provider;
	Double[] SpinVelocity;
	Double[] SpinAcceleration;
    }

    class JSONSpaceObject
    {
	double Mass;
	double Area;
	JSONFacet[] Facets;
	JSONSolarArray SolarArray;
	JSONAttitude Attitude;
    }

    class JSONPropagation
    {
	String Start;
	String End;
	double Step;
	double[] InitialState;
	String[] InitialTLE;
	String InitialStateFrame;
	String InertialFrame;
    }

    class JSONIntegration
    {
	Double MinTimeStep;
	Double MaxTimeStep;
	Double AbsTolerance;
	Double RelTolerance;

	public JSONIntegration()
	{
	    if (MinTimeStep == null)
		MinTimeStep = 1E-3;
	    if (MaxTimeStep == null)
		MaxTimeStep = 300.0;
	    if (AbsTolerance == null)
		AbsTolerance = 1E-14;
	    if (RelTolerance == null)
		RelTolerance = 1E-12;
	}
    }
    
    class JSONManeuver
    {
	String TriggerEvent;
	String ManeuverType;
	String Time;
	double Duration;
	double Thrust;
	double Isp;
	double[] Params;
    }

    class JSONStation
    {
	double Latitude;
	double Longitude;
	double Altitude;
	double AzimuthBias;
	double ElevationBias;
	double RangeBias;
	double RangeRateBias;
	double RightAscensionBias;
	double DeclinationBias;
	double[] PositionVelocityBias;
    }

    class JSONMeasurement
    {
	boolean TwoWay;
	double[] Error;
	String ReferenceFrame;
    }

    class JSONEstimation
    {
	String Filter;
	double[] Covariance;
	double[] ProcessNoise;
	double DMCCorrTime;
	double DMCSigmaPert;
	JSONParameter DMCAcceleration;
    }

    class JSONSimulation
    {
	Boolean SimulateMeasurements;
	Boolean SkipUnobservable;
	Boolean IncludeExtras;
	Boolean IncludeStationState;
    }

    class EstimatedParameter
    {
	String name;
	double min;
	double max;
	double value;

	public EstimatedParameter(String n, double mi, double ma, double v)
	{
	    name = n;
	    min = mi;
	    max = ma;
	    value = v;
	}
    }

    JSONGravity Gravity;
    JSONOceanTides OceanTides;
    JSONDrag Drag;
    JSONSolidTides SolidTides;
    JSONThirdBodies ThirdBodies;
    JSONRadiationPressure RadiationPressure;
    JSONSpaceObject SpaceObject;
    JSONPropagation Propagation;
    JSONIntegration Integration;
    JSONManeuver[] Maneuvers;
    Map<String, JSONStation> Stations;
    Map<String, JSONMeasurement> Measurements;
    JSONEstimation Estimation;
    JSONSimulation Simulation;

    Atmosphere atmmodel = null;

    HashMap<String, GroundStation> stations;
    ArrayList<ForceModel> forces;
    ArrayList<Settings.EstimatedParameter> estparams;

    Frame propframe;

    public static Settings loadJSON(String json)
    {
	Settings set = new Gson().fromJson(json, Settings.class);
	if (set.Integration == null)
	    set.Integration = set.new JSONIntegration();
	if (set.Propagation.InertialFrame != null && set.Propagation.InertialFrame.equals("GCRF"))
	    set.propframe = DataManager.gcrf;
	else
	    set.propframe = DataManager.eme2000;

	set.loadGroundStations();
	set.loadForces();
	set.loadEstimatedParameters();
	return(set);
    }

    private void loadGroundStations()
    {
	if (Stations != null)
	    stations = new HashMap<String, GroundStation>(Stations.size());
	else
	{
	    stations = new HashMap<String, GroundStation>();
	    return;
	}

	for (Map.Entry<String, JSONStation> kv : Stations.entrySet())
	{
	    String k = kv.getKey();
	    JSONStation v = kv.getValue();
	    GroundStation sta = new GroundStation(new TopocentricFrame(new OneAxisEllipsoid(
									   Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
									   Constants.WGS84_EARTH_FLATTENING, DataManager.itrf),
								       new GeodeticPoint(v.Latitude, v.Longitude, v.Altitude), k));
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
	if (Gravity.Degree >= 2 && Gravity.Order >= 0)
	{
	    grav = GravityFieldFactory.getNormalizedProvider(Gravity.Degree, Gravity.Order);
	    forces.add(new HolmesFeatherstoneAttractionModel(DataManager.itrf, grav));
	}
	else
	    forces.add(new NewtonianAttraction(Constants.EGM96_EARTH_MU));

	if (OceanTides.Degree >= 0 && OceanTides.Order >= 0)
	    forces.add(new OceanTides(DataManager.itrf, Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
				      Constants.EGM96_EARTH_MU, OceanTides.Degree,
				      OceanTides.Order, IERSConventions.IERS_2010, DataManager.ut1scale));

	if ((SolidTides.Sun || SolidTides.Moon) && grav != null)
	    forces.add(new org.orekit.forces.gravity.SolidTides(DataManager.itrf, Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
								Constants.EGM96_EARTH_MU, grav.getTideSystem(),
								IERSConventions.IERS_2010, DataManager.ut1scale,
								CelestialBodyFactory.getSun(), CelestialBodyFactory.getMoon()));

	if (ThirdBodies.Sun)
	    forces.add(new ThirdBodyAttraction(CelestialBodyFactory.getSun()));
	if (ThirdBodies.Moon)
	    forces.add(new ThirdBodyAttraction(CelestialBodyFactory.getMoon()));

	DragSensitive dragsc = null;
	RadiationSensitive radnsc = null;
	if (SpaceObject.Facets != null && SpaceObject.SolarArray != null)
	{
	    BoxAndSolarArraySpacecraft.Facet[] facets = new BoxAndSolarArraySpacecraft.Facet[
		SpaceObject.Facets.length];
	    for (int i = 0; i < SpaceObject.Facets.length; i++)
		facets[i] = new BoxAndSolarArraySpacecraft.Facet(new Vector3D(SpaceObject.Facets[i].Normal),
								 SpaceObject.Facets[i].Area);

	    dragsc = new BoxAndSolarArraySpacecraft(facets, CelestialBodyFactory.getSun(), SpaceObject.SolarArray.Area,
						    new Vector3D(SpaceObject.SolarArray.Axis), Drag.Coefficient.Value,
						    RadiationPressure.Cabsorption, RadiationPressure.Creflection.Value);
	    radnsc = new BoxAndSolarArraySpacecraft(facets, CelestialBodyFactory.getSun(), SpaceObject.SolarArray.Area,
						    new Vector3D(SpaceObject.SolarArray.Axis), Drag.Coefficient.Value,
						    RadiationPressure.Cabsorption, RadiationPressure.Creflection.Value);
	}
	else
	{
	    dragsc = new IsotropicDrag(SpaceObject.Area, Drag.Coefficient.Value);
	    radnsc = new IsotropicRadiationSingleCoefficient(SpaceObject.Area, RadiationPressure.Creflection.Value,
							     RadiationPressure.Creflection.Min, RadiationPressure.Creflection.Max);
	}

	if (Drag.Model.equals("Exponential"))
	{
	    atmmodel = new SimpleExponentialAtmosphere(new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
									    Constants.WGS84_EARTH_FLATTENING, DataManager.itrf),
						       Drag.ExpRho0, Drag.ExpH0, Drag.ExpHScale);
	}
	else if (Drag.Model.equals("MSISE"))
	{
	    int apflag = 1;
	    if (Drag.MSISEFlags != null)
	    {
		for (int i = 0; i < Drag.MSISEFlags.length; i++)
		{
		    if (Drag.MSISEFlags[i][0] == 9)
			apflag = Drag.MSISEFlags[i][1];
		}
	    }

	    atmmodel = new NRLMSISE00(new MSISEInputs(DataManager.msisedata.mindate, DataManager.msisedata.maxdate,
						      DataManager.msisedata.data, apflag),
				      CelestialBodyFactory.getSun(), new OneAxisEllipsoid(
					  Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
					  Constants.WGS84_EARTH_FLATTENING, DataManager.itrf));
	    if (Drag.MSISEFlags != null)
	    {
		for (int i = 0; i < Drag.MSISEFlags.length; i++)
		    atmmodel = ((NRLMSISE00) atmmodel).withSwitch(Drag.MSISEFlags[i][0], Drag.MSISEFlags[i][1]);
	    }
	}

	if (atmmodel != null)
	    forces.add(new DragForce(atmmodel, dragsc));

	if (RadiationPressure.Sun)
	    forces.add(new SolarRadiationPressure(CelestialBodyFactory.getSun(), Constants.WGS84_EARTH_EQUATORIAL_RADIUS, radnsc));

	if (Maneuvers == null)
	    return;
	for (JSONManeuver m : Maneuvers)
	{
	    if (m.TriggerEvent == null)
		m.TriggerEvent = "DateTime";
	    if (m.ManeuverType == null)
		m.ManeuverType = "ConstantThrust";
	    if (m.TriggerEvent.equals("DateTime") && m.ManeuverType.equals("ConstantThrust"))
		forces.add(new ConstantThrustManeuver(new AbsoluteDate(DateTimeComponents.parseDateTime(m.Time), DataManager.utcscale),
						      m.Duration, m.Thrust, m.Isp, new Vector3D(m.Params)));
	}
    }

    private void loadEstimatedParameters()
    {
	estparams = new ArrayList<EstimatedParameter>();
	if (Drag.Coefficient.Estimation != null &&
	    Drag.Coefficient.Estimation.equals("Estimate"))
	    estparams.add(new EstimatedParameter(DragSensitive.DRAG_COEFFICIENT, Drag.Coefficient.Min,
						 Drag.Coefficient.Max, Drag.Coefficient.Value));

	if (RadiationPressure.Creflection.Estimation != null &&
	    RadiationPressure.Creflection.Estimation.equals("Estimate"))
	    estparams.add(new EstimatedParameter(RadiationSensitive.REFLECTION_COEFFICIENT, RadiationPressure.Creflection.Min,
						 RadiationPressure.Creflection.Max, RadiationPressure.Creflection.Value));

	if (Estimation != null && Estimation.DMCCorrTime > 0.0 && Estimation.DMCSigmaPert > 0.0)
	{
	    for (int i = 0; i < 3; i++)
		estparams.add(new EstimatedParameter(org.astria.Estimation.DMC_ACC_ESTM+i, Estimation.DMCAcceleration.Min,
						     Estimation.DMCAcceleration.Max, Estimation.DMCAcceleration.Value));
	}
    }

    public double[] getInitialState()
    {
	PVCoordinates topv;
	double[] state0 = Propagation.InitialState;

	if (state0 == null)
	{
	    TLE parser = new TLE(Propagation.InitialTLE[0], Propagation.InitialTLE[1]);
	    TLEPropagator prop = TLEPropagator.selectExtrapolator(parser);
	    AbsoluteDate epoch;
	    if (Propagation.Start != null)
		epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(Propagation.Start), DataManager.utcscale);
	    else
	    {
		epoch = parser.getDate().shiftedBy(0.0);
		Propagation.Start = new String(epoch.toString()) + "Z";
	    }
	    topv = prop.getPVCoordinates(epoch, propframe);
	}
	else
	{
	    Frame fromframe = DataManager.eme2000;
	    if (Propagation.InitialStateFrame != null)
	    {
		if (Propagation.InitialStateFrame.equals("GCRF"))
		    fromframe = DataManager.gcrf;
		else if (Propagation.InitialStateFrame.equals("ITRF"))
		    fromframe = DataManager.itrf;
	    }

	    AbsoluteDate epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(Propagation.Start), DataManager.utcscale);
	    Transform xfm = fromframe.getTransformTo(propframe, epoch);
	    PVCoordinates frompv = new PVCoordinates(new Vector3D(state0[0], state0[1], state0[2]),
						     new Vector3D(state0[3], state0[4], state0[5]));
	    topv = xfm.transformPVCoordinates(frompv);
	}

	Vector3D p = topv.getPosition();
	Vector3D v = topv.getVelocity();
	state0 = new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()};
	double[] X0 = new double[estparams.size() + 6];
	for (int i = 0; i < X0.length; i++)
	{
	    if (i < 6)
		X0[i] = state0[i];
	    else
		X0[i] = estparams.get(i - 6).value;
	}

	return(X0);
    }

    public AttitudeProvider getAttitudeProvider()
    {
	AttitudeProvider attpro = null;
	if (SpaceObject.Attitude == null)
	    return(attpro);

	OneAxisEllipsoid shape = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
						      Constants.WGS84_EARTH_FLATTENING, DataManager.itrf);

	if (SpaceObject.Attitude.Provider.equals("NadirPointing"))
	    attpro = new NadirPointing(propframe, shape);
	if (SpaceObject.Attitude.Provider.equals("BodyCenterPointing"))
	    attpro = new BodyCenterPointing(propframe, shape);
	if (SpaceObject.Attitude.Provider.equals("FixedRate") && SpaceObject.Attitude.SpinVelocity != null &&
	    SpaceObject.Attitude.SpinAcceleration != null)
	{
	    double[] X0 = Propagation.InitialState;
	    AbsoluteDate t0 = new AbsoluteDate(DateTimeComponents.parseDateTime(Propagation.Start),
					       DataManager.utcscale);
	    KeplerianPropagator prop = new KeplerianPropagator(new CartesianOrbit(new PVCoordinates(new Vector3D(X0[0], X0[1], X0[2]),
												    new Vector3D(X0[3], X0[4], X0[5])),
										  propframe, t0, Constants.EGM96_EARTH_MU));
	    LocalOrbitalFrame lof = new LocalOrbitalFrame(propframe, LOFType.VVLH, prop, "");
	    attpro = new FixedRate(new Attitude(t0, lof, Rotation.IDENTITY,
						new Vector3D(Stream.of(SpaceObject.Attitude.SpinVelocity).
							     mapToDouble(Double::doubleValue).toArray()),
						new Vector3D(Stream.of(SpaceObject.Attitude.SpinAcceleration).
							     mapToDouble(Double::doubleValue).toArray())));
	}

	return(attpro);
    }
    
    public RealMatrix getProcessNoiseMatrix(double t)
    {
	int i;
	double t2 = t*t;
	double t3 = t2*t;
	double t4 = t3*t;
	double[][] Q = new double[estparams.size() + 6][estparams.size() + 6];

	if (Estimation.DMCCorrTime <= 0.0 || Estimation.DMCSigmaPert <= 0.0)
	{
	    double[] P = Estimation.ProcessNoise;
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

	int N = estparams.size() - 3;
	double b = 1.0/Estimation.DMCCorrTime;
	double b2 = b*b;
	double b3 = b2*b;
	double b4 = b3*b;
	double b5 = b4*b;
	double et = FastMath.exp(-1.0*b*t);
	double e2t = et*et;
	double s2 = Estimation.DMCSigmaPert*Estimation.DMCSigmaPert;

	double Q00 = s2*(t3/(3*b2)-t2/b3+t*(1-2*et)/b4+0.5*(1-e2t)/b5); // pos-pos
	double Q01 = s2*(0.5*t2/b2-t*(1-et)/b3+(1-et)/b4-0.5*(1-e2t)/b4); // pos-vel
	double Q02 = s2*(0.5*(1-e2t)/b3-t*et/b2); // pos-acc
	double Q11 = s2*(t/b2-2*(1-et)/b3+0.5*(1-e2t)/b3); // vel-vel
	double Q12 = s2*(0.5*(1+e2t)/b2-et/b2); // vel-acc
	double Q22 = 0.5*s2*(1-e2t)/b; // acc-acc

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
	if (Maneuvers == null)
	    return;
	for (JSONManeuver m : Maneuvers)
	{
	    AbsoluteDate time = new AbsoluteDate(DateTimeComponents.parseDateTime(m.Time), DataManager.utcscale);
	    if (m.TriggerEvent.equals("DateTime"))
	    {
		if (m.ManeuverType.equals("ConstantThrust"))
		    continue;
		prop.addEventDetector(new DateDetector(time).withHandler(new EventHandling<DateDetector>(m)));
	    }
	}
    }
}
