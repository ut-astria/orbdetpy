/*
 * Simulation.java - Functions to simulate noisy measurements.
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

import com.google.gson.GsonBuilder;
import java.util.ArrayList;
import java.util.Map;
import java.util.Random;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.forces.ForceModel;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.measurements.RangeRate;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public class Simulation
{
    private Settings simcfg;

    public Simulation(String cfgjson)
    {
	simcfg = Settings.loadJSON(cfgjson);
    }

    public String simulateMeasurements()
    {
	Random rand = new Random();
	boolean simulmeas = true, skipunobs = true, inclextra = false, inclstapos = false;
	if (simcfg.cfgSimulation != null)
	{
	    if (simcfg.cfgSimulation.SimulateMeasurements != null && !simcfg.cfgSimulation.SimulateMeasurements)
		simulmeas = false;
	    if (simcfg.cfgSimulation.SkipUnobservable != null && !simcfg.cfgSimulation.SkipUnobservable)
		skipunobs = false;
	    if (simcfg.cfgSimulation.IncludeExtras != null && simcfg.cfgSimulation.IncludeExtras)
		inclextra = true;
	    if (simcfg.cfgSimulation.IncludeStationState != null && simcfg.cfgSimulation.IncludeStationState)
		inclstapos = true;
	}

	double[] Xi = simcfg.getInitialState();
	AbsoluteDate tm = new AbsoluteDate(DateTimeComponents.parseDateTime(simcfg.cfgPropagation.Start),
					   DataManager.utcscale);
	AbsoluteDate prend = new AbsoluteDate(DateTimeComponents.parseDateTime(simcfg.cfgPropagation.End),
					      DataManager.utcscale);

	NumericalPropagator prop = new NumericalPropagator(
	    new DormandPrince853Integrator(simcfg.cfgIntegration.MinTimeStep, simcfg.cfgIntegration.MaxTimeStep,
					   simcfg.cfgIntegration.AbsTolerance, simcfg.cfgIntegration.RelTolerance));
	for (ForceModel fm : simcfg.forces)
	    prop.addForceModel(fm);

	AttitudeProvider attpro = simcfg.getAttitudeProvider();
	if (attpro != null)
	    prop.setAttitudeProvider(attpro);

	SpacecraftState sstate0 = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(Xi[0], Xi[1], Xi[2]),
											   new Vector3D(Xi[3], Xi[4], Xi[5])),
									 simcfg.propframe, tm, Constants.EGM96_EARTH_MU),
						      simcfg.cfgSpaceObject.Mass);
	prop.setInitialState(sstate0);
	simcfg.addEventHandlers(prop, sstate0);
	boolean evthandlers = prop.getEventsDetectors().size() > 0;

	double[] obs, azel;
	double[] zeros = new double[]{0.0, 0.0};
	double[] ones = new double[]{1.0, 1.0};
	ObservableSatellite obssat = new ObservableSatellite(0);
	Measurements meas = new Measurements();
	ArrayList<Measurements.SimulatedMeasurement> mall = new ArrayList<Measurements.SimulatedMeasurement>();

	while (true)
	{
	    SpacecraftState[] sta = new SpacecraftState[]{prop.propagate(tm)};
	    AbsoluteDate proptm = sta[0].getDate();
	    TimeStampedPVCoordinates pvc = sta[0].getPVCoordinates();
	    Vector3D pos = pvc.getPosition();
	    Vector3D vel = pvc.getVelocity();
	    Vector3D acc = pvc.getAcceleration();

	    Measurements.SimulatedMeasurement json = meas.new SimulatedMeasurement();
	    json.Time = proptm.toString() + "Z";
	    json.TrueState = meas.new State();
	    json.TrueState.Cartesian = new double[]{pos.getX(), pos.getY(), pos.getZ(), vel.getX(), vel.getY(), vel.getZ(),
						    acc.getX(), acc.getY(), acc.getZ()};

	    if (inclextra)
	    {
		Orbit orb = sta[0].getOrbit();
		KeplerianOrbit keporb = new KeplerianOrbit(orb);
		json.TrueState.Kepler = meas.new Kepler(keporb.getA(), keporb.getE(), keporb.getI(),
							    keporb.getRightAscensionOfAscendingNode(),
							    keporb.getPerigeeArgument(), keporb.getMeanAnomaly());
		json.TrueState.Equinoctial = meas.new Equinoctial(orb.getA(), orb.getEquinoctialEx(),
								      orb.getEquinoctialEy(), orb.getHx(),
								      orb.getHy(), orb.getLM());

		if (simcfg.atmmodel != null)
		    json.AtmDensity = simcfg.atmmodel.getDensity(proptm, pos, simcfg.propframe);
		for (ForceModel fmod : simcfg.forces)
		{
		    double[] facc = fmod.acceleration(sta[0], fmod.getParameters()).toArray();
		    String ftype = fmod.getClass().getSimpleName();
		    if (ftype.equals("HolmesFeatherstoneAttractionModel") || ftype.equals("NewtonianAttraction"))
			json.AccGravity = facc;
		    if (ftype.equals("DragForce"))
			json.AccDrag = facc;
		    if (ftype.equals("SolidTides"))
			json.AccSolidTides = facc;
		    if (ftype.equals("OceanTides"))
			json.AccOceanTides = facc;
		    if (ftype.equals("SolarRadiationPressure"))
			json.AccRadiationPressure = facc;
		    if (ftype.equals("ConstantThrustManeuver"))
			json.AccThrust = facc;
		    if (ftype.equals("ThirdBodyAttraction"))
		    {
			if (json.AccThirdBodies == null)
			    json.AccThirdBodies = facc;
			else
			{
			    for (int ii = 0; ii < 3; ii++)
				json.AccThirdBodies[ii] += facc[ii];
			}
		    }
		}
	    }

	    if (simulmeas)
	    {
		boolean added = false;
		for (Map.Entry<String, GroundStation> kv : simcfg.stations.entrySet())
		{
		    GroundStation gst = kv.getValue();
		    Settings.Station jsn = simcfg.cfgStations.get(kv.getKey());
		    azel = new AngularAzEl(gst, proptm, zeros, zeros, ones, obssat).estimate(0, 0, sta).getEstimatedValue();
		    if (skipunobs && azel[1] <= 5E-6)
			continue;

		    Measurements.SimulatedMeasurement clone = meas.new SimulatedMeasurement(json);
		    clone.Station = kv.getKey();
		    if (inclstapos)
		    {
			pvc = gst.getBaseFrame().getPVCoordinates(proptm, simcfg.propframe);
			pos = pvc.getPosition();
			vel = pvc.getVelocity();
			acc = pvc.getAcceleration();
			clone.StationState = new double[]{pos.getX(), pos.getY(), pos.getZ(), vel.getX(), vel.getY(),
							  vel.getZ(), acc.getX(), acc.getY(), acc.getZ()};
		    }

		    for (Map.Entry<String, Settings.Measurement> nvp : simcfg.cfgMeasurements.entrySet())
		    {
			String name = nvp.getKey();
			Settings.Measurement val = nvp.getValue();

			if (name.equals("Range"))
			{
			    obs = new Range(gst, val.TwoWay, proptm, 0.0, 0.0, 1.0, obssat).estimate(0, 0, sta).getEstimatedValue();
			    clone.Range = obs[0] + rand.nextGaussian()*val.Error[0] + jsn.RangeBias;
			}
			else if (name.equals("RangeRate"))
			{
			    obs = new RangeRate(gst, proptm, 0.0, 0.0, 1.0, val.TwoWay, obssat).estimate(0, 0, sta).getEstimatedValue();
			    clone.RangeRate = obs[0] + rand.nextGaussian()*val.Error[0] + jsn.RangeRateBias;
			}
			else if (name.equals("RightAscension") || name.equals("Declination") && clone.RightAscension == null)
			{
			    obs = new AngularRaDec(gst, simcfg.propframe, proptm, zeros, zeros, ones,
						   obssat).estimate(0, 0, sta).getEstimatedValue();
			    clone.RightAscension = obs[0] + rand.nextGaussian()*val.Error[0] + jsn.RightAscensionBias;
			    clone.Declination = obs[1] + rand.nextGaussian()*val.Error[0] + jsn.DeclinationBias;
			}
			else if (name.equals("Azimuth") || name.equals("Elevation") && clone.Azimuth == null)
			{
			    clone.Azimuth = azel[0] + rand.nextGaussian()*val.Error[0] + jsn.AzimuthBias;
			    clone.Elevation = azel[1] + rand.nextGaussian()*val.Error[0] + jsn.ElevationBias;
			}
			else if (name.equals("Position"))
			{
			    clone.Position = new Double[3];
			    for (int i = 0; i < 3; i++)
				clone.Position[i] = clone.TrueState.Cartesian[i] +
				    rand.nextGaussian()*val.Error[i] + jsn.PositionBias[i];
			}
			else if (name.equals("PositionVelocity"))
			{
			    clone.PositionVelocity = new Double[6];
			    for (int i = 0; i < 6; i++)
				clone.PositionVelocity[i] = clone.TrueState.Cartesian[i] +
				    rand.nextGaussian()*val.Error[i] + jsn.PositionVelocityBias[i];
			}
		    }

		    added = true;
		    mall.add(clone);
		}

		if (!added)
		    mall.add(json);
	    }
	    else
		mall.add(json);

	    double dt = prend.durationFrom(proptm);
	    if (evthandlers && !proptm.equals(tm))
		break;
	    if (simcfg.cfgPropagation.Step >= 0.0)
		tm = new AbsoluteDate(tm, FastMath.min(dt, simcfg.cfgPropagation.Step));
	    else
	    {
		tm = new AbsoluteDate(tm, FastMath.max(dt, simcfg.cfgPropagation.Step));
		dt = -dt;
	    }
	    if (dt <= 0.0)
		break;
	}

	Measurements.Measurement[] rawmeas = mall.toArray(new Measurements.SimulatedMeasurement[0]);
	return(new GsonBuilder().setPrettyPrinting().create().toJson(rawmeas));
    }
}
