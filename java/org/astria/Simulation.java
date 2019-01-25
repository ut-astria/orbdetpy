/*
 * Simulation.java - Functions to simulate noisy measurements.
 * Copyright (C) 2018 Shiva Iyer <shiva.iyer AT utexas DOT edu>
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
import org.astria.DataManager;
import org.astria.Measurements;
import org.astria.Settings;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.forces.ForceModel;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.GroundStation;
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

    public Simulation(String cfgjson) throws Exception
    {
	simcfg = Settings.loadJSON(cfgjson);
	if (simcfg.Propagation.Step <= 0.0)
	    simcfg.Propagation.Step = 60.0;
    }

    public String simulateMeasurements()
    {
	Random rand = new Random();

	double[] Xi = simcfg.Propagation.InitialState;
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

	prop.setInitialState(new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(Xi[0], Xi[1], Xi[2]),
										      new Vector3D(Xi[3], Xi[4], Xi[5])),
								    DataManager.eme2000, tm, Constants.EGM96_EARTH_MU),
						 simcfg.SpaceObject.Mass));

	Measurements meas = new Measurements();
	ArrayList<Measurements.JSONMeasurement> mall = new ArrayList<Measurements.JSONMeasurement>();

	while (true)
	{
	    SpacecraftState[] sta = new SpacecraftState[]{prop.propagate(tm)};

	    for (Map.Entry<String, GroundStation> kv : simcfg.stations.entrySet())
	    {
		String str = kv.getKey();
		GroundStation obj = kv.getValue();
		double[] obs = new AngularAzEl(obj, tm, new double[]{0.0, 0.0}, new double[]{5E-6, 5E-6},
					       new double[]{1.0, 1.0}).estimate(1, 1, sta).getEstimatedValue();
		if ((simcfg.Simulation == null || simcfg.Simulation.SkipUnobservable) && obs[1] <= 5E-6)
		    continue;

		Orbit orb = sta[0].getOrbit();
		KeplerianOrbit keporb = new KeplerianOrbit(orb);
		TimeStampedPVCoordinates pvc = sta[0].getPVCoordinates();
		Vector3D pos = pvc.getPosition();
		Vector3D vel = pvc.getVelocity();
		Vector3D acc = pvc.getAcceleration();

		Measurements.JSONMeasurement json = meas.new JSONMeasurement();
		json.Time = tm.toString() + "Z";
		json.Station = str;

		json.TrueState = meas.new JSONState();
		json.TrueState.Cartesian = new double[]{pos.getX(), pos.getY(), pos.getZ(),
							vel.getX(), vel.getY(), vel.getZ(),
							acc.getX(), acc.getY(), acc.getZ()};

		json.TrueState.Kepler = meas.new JSONKepler(keporb.getA(), keporb.getE(), keporb.getI(),
							    keporb.getRightAscensionOfAscendingNode(),
							    keporb.getPerigeeArgument(), keporb.getMeanAnomaly());

		json.TrueState.Equinoctial = meas.new JSONEquinoctial(orb.getA(), orb.getEquinoctialEx(),
								      orb.getEquinoctialEy(), orb.getHx(),
								      orb.getHy(), orb.getLM());

		for (Map.Entry<String, Settings.JSONMeasurement> nvp : simcfg.Measurements.entrySet())
		{
		    String name = nvp.getKey();
		    Settings.JSONMeasurement val = nvp.getValue();

		    if (name.equals("Range"))
		    {
			obs = new Range(obj, tm, 0.0, val.Error[0], 1.0, val.TwoWay).
			    estimate(1, 1, sta).getEstimatedValue();
			json.Range = obs[0] + rand.nextGaussian()*val.Error[0];
		    }
		    else if (name.equals("RangeRate"))
		    {
			obs = new RangeRate(obj, tm, 0.0, val.Error[0], 1.0, val.TwoWay).
			    estimate(1, 1, sta).getEstimatedValue();
			json.RangeRate = obs[0] + rand.nextGaussian()*val.Error[0];
		    }
		    else if (name.equals("RightAscension") || name.equals("Declination"))
		    {
			obs = new AngularRaDec(obj, DataManager.eme2000, tm, new double[]{0.0, 0.0},
					       new double[]{val.Error[0], val.Error[0]}, new double[]{1.0, 1.0}).
			    estimate(1, 1, sta).getEstimatedValue();
			json.RightAscension = obs[0] + rand.nextGaussian()*val.Error[0];
			json.Declination = obs[1] + rand.nextGaussian()*val.Error[0];
		    }
		    else if (name.equals("Azimuth") || name.equals("Elevation"))
		    {
			json.Azimuth = obs[0] + rand.nextGaussian()*val.Error[0];
			json.Elevation = obs[1] + rand.nextGaussian()*val.Error[0];
		    }
		}

		mall.add(json);
		break;
	    }

	    double dt = prend.durationFrom(tm);
	    tm = new AbsoluteDate(tm, Math.min(dt, simcfg.Propagation.Step));
	    if (dt <= 0.0)
		break;
	}

	meas.rawmeas = mall.toArray(new Measurements.JSONMeasurement[0]);
	return(new GsonBuilder().setPrettyPrinting().create().toJson(meas.rawmeas));
    }
}
