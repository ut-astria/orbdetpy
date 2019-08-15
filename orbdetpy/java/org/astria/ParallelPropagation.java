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
import org.astria.Settings;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.orekit.forces.ForceModel;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.PropagatorsParallelizer;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.MultiSatStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.Constants;

public class ParallelPropagation
{
    protected MultiSatStepHandler stepHandler;

    public ParallelPropagation(MultiSatStepHandler hnd)
    {
	stepHandler = hnd;
    }

    public List<SpacecraftState> propagate(String[] cfgjson, String propStart,
					   String propEnd, double propStep)
    {
	List<Propagator> props = new ArrayList<Propagator>(cfgjson.length);
	for (int i = 0; i < cfgjson.length; i++)
	    props.add(buildPropagator(Settings.loadJSON(cfgjson[i])));

	List<SpacecraftState> ssta = null;
	PropagatorsParallelizer plel = new PropagatorsParallelizer(props, stepHandler);

	AbsoluteDate tm = new AbsoluteDate(DateTimeComponents.parseDateTime(propStart),
					   DataManager.utcscale);
	AbsoluteDate tmprev = tm.shiftedBy(-0.1);
	AbsoluteDate prend = new AbsoluteDate(DateTimeComponents.parseDateTime(propEnd),
					      DataManager.utcscale);

	while (true)
	{
	    ssta = plel.propagate(tmprev, tm);

	    double dt = prend.durationFrom(tm);
	    tmprev = tm.shiftedBy(0.0);
	    tm = new AbsoluteDate(tm, FastMath.min(dt, propStep));
	    if (dt <= 0.0)
		break;
	}

	return(ssta);
    }

    protected NumericalPropagator buildPropagator(Settings obj)
    {
	double[] Xi = obj.getInitialState();
	AbsoluteDate tm = new AbsoluteDate(DateTimeComponents.parseDateTime(obj.Propagation.Start),
					   DataManager.utcscale);

	NumericalPropagator prop = new NumericalPropagator(new DormandPrince853Integrator(
		obj.Integration.MinTimeStep, obj.Integration.MaxTimeStep,
		obj.Integration.AbsTolerance, obj.Integration.RelTolerance));
	for (ForceModel fm : obj.forces)
	    prop.addForceModel(fm);

	prop.setInitialState(
	    new SpacecraftState(new CartesianOrbit(
				    new PVCoordinates(
					new Vector3D(Xi[0], Xi[1], Xi[2]),
					new Vector3D(Xi[3], Xi[4], Xi[5])),
				    obj.propframe, tm, Constants.EGM96_EARTH_MU),
				obj.SpaceObject.Mass));

	return(prop);
    }
}
