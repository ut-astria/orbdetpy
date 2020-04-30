/*
 * ParallelPropagation.java - Multiple object propagation functions.
 * Copyright (C) 2018-2020 University of Texas
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
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.forces.ForceModel;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitStepInterpolator;
import org.orekit.propagation.sampling.MultiSatStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.AbsolutePVCoordinates;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class ParallelPropagation
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

    private final ArrayList<Settings> configObjs;
    private final ArrayList<Propagator> propagators;

    public ParallelPropagation(ArrayList<Settings> configObjs)
    {
	this.configObjs = configObjs;
	this.propagators = new ArrayList<Propagator>(configObjs.size());
	for (int i = 0; i < configObjs.size(); i++)
	    propagators.add(buildPropagator(configObjs.get(i)));
    }

    public ArrayList<PropagationOutput> propagate() throws Exception
    {
	final Settings obj0 = configObjs.get(0);
	String start = obj0.propStart, end = obj0.propEnd;
	double step = obj0.propStep;
	for (Settings s: configObjs)
	{
	    if (obj0.propStep > 0.0)
	    {
		if (start.compareTo(s.propStart) > 0)
		    start = s.propStart;
		if (end.compareTo(s.propEnd) < 0)
		    end = s.propEnd;
		if (s.propStep != 0.0)
		    step = FastMath.min(step, s.propStep);
	    }
	    else
	    {
		if (start.compareTo(s.propStart) < 0)
		    start = s.propStart;
		if (end.compareTo(s.propEnd) > 0)
		    end = s.propEnd;
		if (s.propStep != 0.0)
		    step = FastMath.max(step, s.propStep);
	    }
	}

	AbsoluteDate tmFinal = DataManager.parseDateTime(start);
	AbsoluteDate tmStart = new AbsoluteDate(tmFinal, -0.1);
	final AbsoluteDate prend = DataManager.parseDateTime(end);
	final ArrayList<PropagationOutput> propOutput = new ArrayList<PropagationOutput>(
	    (int) FastMath.abs(prend.durationFrom(tmFinal)/step) + 2);

	while (true)
	{
	    final AbsoluteDate t0 = tmStart;
	    final AbsoluteDate t1 = tmFinal;
	    final PropagationOutput pout = new PropagationOutput(tmFinal, propagators.size());
	    propOutput.add(pout);

	    for (int i = 0; i < propagators.size(); i++)
	    {
		final Propagator prop = propagators.get(i);
		final SpacecraftState state = DataManager.threadPool.submit(()->prop.propagate(t0, t1)).get();
		tmStart = state.getDate();
		pout.addState(state.getPVCoordinates(FramesFactory.getFrame(Predefined.valueOf(configObjs.get(i).propInertialFrame))));
	    }

	    double dt = prend.durationFrom(tmStart);
	    if (step >= 0.0)
		tmFinal = new AbsoluteDate(tmFinal, FastMath.min(dt, step));
	    else
	    {
		tmFinal = new AbsoluteDate(tmFinal, FastMath.max(dt, step));
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
	    final TLE parser = new TLE(cfg.propInitialTLE[0], cfg.propInitialTLE[1]);
	    if (cfg.propStart == null || cfg.propStart.length() == 0)
		cfg.propStart = DataManager.getUTCString(parser.getDate());
	    if (cfg.propEnd == null || cfg.propEnd.length() == 0)
		cfg.propEnd = cfg.propStart;
	    return(TLEPropagator.selectExtrapolator(parser));
	}

	final NumericalPropagator prop = new NumericalPropagator(
	    new DormandPrince853Integrator(cfg.integMinTimeStep, cfg.integMaxTimeStep,
					   cfg.integAbsTolerance, cfg.integRelTolerance));
	for (ForceModel fm : cfg.forces)
	    prop.addForceModel(fm);

	final AttitudeProvider attpro = cfg.getAttitudeProvider();
	if (attpro != null)
	    prop.setAttitudeProvider(attpro);

	final double[] Xi = cfg.getInitialState();
	final AbsoluteDate tm = DataManager.parseDateTime(cfg.propStart);
	prop.setInitialState(new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(Xi[0], Xi[1], Xi[2]),
										      new Vector3D(Xi[3], Xi[4], Xi[5])),
								    cfg.propFrame, tm, Constants.EGM96_EARTH_MU), cfg.rsoMass));
	return(prop);
    }
}
