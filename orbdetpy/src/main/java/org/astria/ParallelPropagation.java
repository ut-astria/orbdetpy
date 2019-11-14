/*
 * ParallelPropagation.java - Multiple object propagation functions.
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
import java.util.List;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.forces.ForceModel;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.PropagatorsParallelizer;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitStepInterpolator;
import org.orekit.propagation.sampling.MultiSatStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.AbsolutePVCoordinates;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class ParallelPropagation implements MultiSatStepHandler
{
    public static class PropagationOutput
    {
	public String time;
	public ArrayList<double[]> states;

	public PropagationOutput(AbsoluteDate time, int objCount)
	{
	    this.time = DataManager.getUTCString(time);
	    this.states = new ArrayList<double[]>(objCount);
	}

	public void addState(TimeStampedPVCoordinates pva)
	{
	    Vector3D p = pva.getPosition();
	    Vector3D v = pva.getVelocity();
	    this.states.add(new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()});
	}
    }

    private final List<Settings> configObjs;
    private final List<Propagator> props;
    private ArrayList<PropagationOutput> propOutput;

    public ParallelPropagation(List<Settings> configObjs)
    {
	this.configObjs = configObjs;
	this.props = new ArrayList<Propagator>(configObjs.size());
	for (int i = 0; i < configObjs.size(); i++)
	    props.add(buildPropagator(configObjs.get(i)));
    }

    public ArrayList<PropagationOutput> propagate()
    {
	Settings obj0 = configObjs.get(0);
	PropagatorsParallelizer plel = new PropagatorsParallelizer(props, this);
	AbsoluteDate tm = new AbsoluteDate(DateTimeComponents.parseDateTime(obj0.propStart), DataManager.getTimeScale("UTC"));
	AbsoluteDate proptm = new AbsoluteDate(tm, -0.1);
	AbsoluteDate prend = new AbsoluteDate(DateTimeComponents.parseDateTime(obj0.propEnd), DataManager.getTimeScale("UTC"));

	List<SpacecraftState> staList = null;
	propOutput = new ArrayList<PropagationOutput>((int) FastMath.abs(prend.durationFrom(tm)/obj0.propStep) + 2);
	while (true)
	{
	    staList = plel.propagate(proptm, tm);
	    proptm = staList.get(0).getDate();

	    PropagationOutput pout = new PropagationOutput(proptm, configObjs.size());
	    propOutput.add(pout);
	    for (SpacecraftState sta : staList)
		pout.addState(sta.getPVCoordinates(DataManager.getFrame(obj0.propInertialFrame)));

	    double dt = prend.durationFrom(proptm);
	    if (obj0.propStep >= 0.0)
		tm = new AbsoluteDate(tm, FastMath.min(dt, obj0.propStep));
	    else
	    {
		tm = new AbsoluteDate(tm, FastMath.max(dt, obj0.propStep));
		dt = -dt;
	    }
	    if (dt <= 0.0)
		break;
	}

	return(propOutput);
    }

    private Propagator buildPropagator(Settings cfg)
    {
	if (cfg.propInitialTLE != null && cfg.propInitialTLE[0] != null && cfg.propInitialTLE[1] != null)
	{
	    TLE parser = new TLE(cfg.propInitialTLE[0], cfg.propInitialTLE[1]);
	    return(TLEPropagator.selectExtrapolator(parser));
	}

	NumericalPropagator prop = new NumericalPropagator(
	    new DormandPrince853Integrator(cfg.integMinTimeStep, cfg.integMaxTimeStep, cfg.integAbsTolerance, cfg.integRelTolerance));
	for (ForceModel fm : cfg.forces)
	    prop.addForceModel(fm);

	AttitudeProvider attpro = cfg.getAttitudeProvider();
	if (attpro != null)
	    prop.setAttitudeProvider(attpro);

	double[] Xi = cfg.getInitialState();
	AbsoluteDate tm = new AbsoluteDate(DateTimeComponents.parseDateTime(cfg.propStart), DataManager.getTimeScale("UTC"));
	prop.setInitialState(new SpacecraftState(
				 new CartesianOrbit(new PVCoordinates(new Vector3D(Xi[0], Xi[1], Xi[2]),
								      new Vector3D(Xi[3], Xi[4], Xi[5])),
						    cfg.propFrame, tm, Constants.EGM96_EARTH_MU), cfg.rsoMass));
	return(prop);
    }

    public void handleStep(List<OrekitStepInterpolator> interpolators, boolean isLast)
    {
    }
}
