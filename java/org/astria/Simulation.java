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
	if (simcfg.Simulation != null)
	{
	    if (simcfg.Simulation.SimulateMeasurements != null && !simcfg.Simulation.SimulateMeasurements)
		simulmeas = false;
	    if (simcfg.Simulation.SkipUnobservable != null && !simcfg.Simulation.SkipUnobservable)
		skipunobs = false;
	    if (simcfg.Simulation.IncludeExtras != null && simcfg.Simulation.IncludeExtras)
		inclextra = true;
	    if (simcfg.Simulation.IncludeStationState != null && simcfg.Simulation.IncludeStationState)
		inclstapos = true;
	}

	double[] Xi = simcfg.getInitialState();
	AbsoluteDate tm = new AbsoluteDate(DateTimeComponents.parseDateTime(simcfg.Propagation.Start),
					   DataManager.utcscale);
	AbsoluteDate prend = new AbsoluteDate(DateTimeComponents.parseDateTime(simcfg.Propagation.End),
					      DataManager.utcscale);

	NumericalPropagator prop = new NumericalPropagator(
	    new DormandPrince853Integrator(simcfg.Integration.MinTimeStep, simcfg.Integration.MaxTimeStep,
					   simcfg.Integration.AbsTolerance, simcfg.Integration.RelTolerance));
	for (ForceModel fm : simcfg.forces)
	    prop.addForceModel(fm);

	AttitudeProvider attpro = simcfg.getAttitudeProvider();
	if (attpro != null)
	    prop.setAttitudeProvider(attpro);

	SpacecraftState sstate0 = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(Xi[0], Xi[1], Xi[2]),
											   new Vector3D(Xi[3], Xi[4], Xi[5])),
									 simcfg.propframe, tm, Constants.EGM96_EARTH_MU),
						      simcfg.SpaceObject.Mass);
	prop.setInitialState(sstate0);
	simcfg.addEventHandlers(prop, sstate0);

	Measurements meas = new Measurements();
	ArrayList<Measurements.JSONSimulatedMeasurement> mall = new ArrayList<Measurements.JSONSimulatedMeasurement>();
	while (true)
	{
	    SpacecraftState[] sta = new SpacecraftState[]{prop.propagate(tm)};
	    TimeStampedPVCoordinates pvc = sta[0].getPVCoordinates();
	    Vector3D pos = pvc.getPosition();
	    Vector3D vel = pvc.getVelocity();
	    Vector3D acc = pvc.getAcceleration();

	    Measurements.JSONSimulatedMeasurement json = meas.new JSONSimulatedMeasurement();
	    json.Time = tm.toString() + "Z";
	    json.TrueState = meas.new JSONState();
	    json.TrueState.Cartesian = new double[]{pos.getX(), pos.getY(), pos.getZ(), vel.getX(), vel.getY(), vel.getZ(),
						    acc.getX(), acc.getY(), acc.getZ()};

	    if (inclextra)
	    {
		Orbit orb = sta[0].getOrbit();
		KeplerianOrbit keporb = new KeplerianOrbit(orb);
		json.TrueState.Kepler = meas.new JSONKepler(keporb.getA(), keporb.getE(), keporb.getI(),
							    keporb.getRightAscensionOfAscendingNode(),
							    keporb.getPerigeeArgument(), keporb.getMeanAnomaly());
		json.TrueState.Equinoctial = meas.new JSONEquinoctial(orb.getA(), orb.getEquinoctialEx(),
								      orb.getEquinoctialEy(), orb.getHx(),
								      orb.getHy(), orb.getLM());

		if (simcfg.atmmodel != null)
		    json.AtmDensity = simcfg.atmmodel.getDensity(tm, pos, simcfg.propframe);

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
		    Settings.JSONStation jsn = simcfg.Stations.get(kv.getKey());
		    double[] obs = new AngularAzEl(gst, tm, new double[]{0.0, 0.0}, new double[]{0.0, 0.0},
						   new double[]{1.0, 1.0}, new ObservableSatellite(0)).
			estimate(1, 1, sta).getEstimatedValue();
		    if (skipunobs && obs[1] <= 5E-6)
			continue;

		    Measurements.JSONSimulatedMeasurement clone = meas.new JSONSimulatedMeasurement(json);
		    clone.Station = kv.getKey();
		    if (inclstapos)
		    {
			pvc = gst.getBaseFrame().getPVCoordinates(tm, simcfg.propframe);
			pos = pvc.getPosition();
			vel = pvc.getVelocity();
			acc = pvc.getAcceleration();
			clone.StationState = new double[]{pos.getX(), pos.getY(), pos.getZ(), vel.getX(), vel.getY(),
							  vel.getZ(), acc.getX(), acc.getY(), acc.getZ()};
		    }

		    for (Map.Entry<String, Settings.JSONMeasurement> nvp : simcfg.Measurements.entrySet())
		    {
			String name = nvp.getKey();
			Settings.JSONMeasurement val = nvp.getValue();

			if (name.equals("Range"))
			{
			    obs = new Range(gst, val.TwoWay, tm, 0.0, 0.0, 1.0, new ObservableSatellite(0)).
				estimate(1, 1, sta).getEstimatedValue();
			    clone.Range = obs[0] + rand.nextGaussian()*val.Error[0] + jsn.RangeBias;
			}
			else if (name.equals("RangeRate"))
			{
			    obs = new RangeRate(gst, tm, 0.0, 0.0, 1.0, val.TwoWay, new ObservableSatellite(0)).
				estimate(1, 1, sta).getEstimatedValue();
			    clone.RangeRate = obs[0] + rand.nextGaussian()*val.Error[0] + jsn.RangeRateBias;
			}
			else if (name.equals("RightAscension") || name.equals("Declination") &&
				 clone.RightAscension == null)
			{
			    obs = new AngularRaDec(gst, simcfg.propframe, tm, new double[]{0.0, 0.0},
						   new double[]{0.0, 0.0}, new double[]{1.0, 1.0},
						   new ObservableSatellite(0)).estimate(1, 1, sta).getEstimatedValue();
			    clone.RightAscension = obs[0] + rand.nextGaussian()*val.Error[0] + jsn.RightAscensionBias;
			    clone.Declination = obs[1] + rand.nextGaussian()*val.Error[0] + jsn.DeclinationBias;
			}
			else if (name.equals("Azimuth") || name.equals("Elevation") && clone.Azimuth == null)
			{
			    obs = new AngularAzEl(gst, tm, new double[]{0.0, 0.0}, new double[]{0.0, 0.0},
						  new double[]{1.0, 1.0}, new ObservableSatellite(0)).
				estimate(1, 1, sta).getEstimatedValue();
			    clone.Azimuth = obs[0] + rand.nextGaussian()*val.Error[0] + jsn.AzimuthBias;
			    clone.Elevation = obs[1] + rand.nextGaussian()*val.Error[0] + jsn.ElevationBias;
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

	    double dt = prend.durationFrom(tm);
	    if (simcfg.Propagation.Step >= 0.0)
		tm = new AbsoluteDate(tm, FastMath.min(dt, simcfg.Propagation.Step));
	    else
	    {
		tm = new AbsoluteDate(tm, FastMath.max(dt, simcfg.Propagation.Step));
		dt = -dt;
	    }
	    if (dt <= 0.0)
		break;
	}

	Measurements.JSONMeasurement[] rawmeas = mall.toArray(new Measurements.JSONSimulatedMeasurement[0]);
	return(new GsonBuilder().setPrettyPrinting().create().toJson(rawmeas));
    }
}
