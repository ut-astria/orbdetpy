/*
 * PropUtil.java - Functions used to speed up propagation in Python.
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

import java.util.ArrayList;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.orekit.bodies.CelestialBody;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.forces.ForceModel;
import org.orekit.forces.gravity.SolidTides;
import org.orekit.forces.gravity.potential.TideSystem;
import org.orekit.frames.Frame;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.UT1Scale;
import org.orekit.utils.Constants;
import org.orekit.utils.IERSConventions;
import org.orekit.utils.PVCoordinates;

public class PropUtil implements OrdinaryDifferentialEquation
{
    private AbsoluteDate epoch;
    private double mass;
    private Frame frame;
    private ForceModel[] forces;
    private ArrayList<double []> params;
    private int stadim;
    private int vecdim;
    private ODEIntegrator ode;
    private Exception except;

    public PropUtil(AbsoluteDate e, double m, Frame r, ForceModel[] f,
		    int sdim) throws Exception
    {
	epoch = e;
	mass = m;
	frame = r;
	forces = f;
	stadim = sdim;
	ode = new DormandPrince853Integrator(1E-3, 300.0, 1E-14, 1E-12);
	params = new ArrayList<>();
	for (int i = 0; i < forces.length; i++)
	    params.add(forces[i].getParameters());
    }

    public static SolidTides solidtides(Frame bdf, double ae, double mu,
					TideSystem tide, IERSConventions conv,
					UT1Scale ut1, boolean sun, boolean moon)
	throws Exception
    {
	ArrayList<CelestialBody> celb = new ArrayList<>();
	if (sun)
	    celb.add(CelestialBodyFactory.getSun());
	if (moon)
	    celb.add(CelestialBodyFactory.getMoon());

	return(new SolidTides(bdf, ae, mu, tide, conv, ut1,
			      celb.toArray(new CelestialBody[celb.size()])));
    }

    public double[] propagate(double t0, double[] X0, double t1)
	throws Exception
    {
	except = null;
	vecdim = X0.length;
	double[] retval = ode.integrate(this, new ODEState(t0, X0), t1).getPrimaryState();
	if (except != null)
	    throw(except);
	return(retval);
    }

    public int getDimension()
    {
	return(vecdim);
    }

    public double[] computeDerivatives(double t, double[] X)
    {
	double[] Xdot = new double[X.length];
	try
	{
	    AbsoluteDate tm = new AbsoluteDate(epoch, t);
	    for (int i = 0; i < X.length; i += stadim)
	    {
		SpacecraftState ss = new SpacecraftState(
		    new CartesianOrbit(new PVCoordinates(
					   new Vector3D(X[i],   X[i+1], X[i+2]),
					   new Vector3D(X[i+3], X[i+4], X[i+5])),
				       frame, tm, Constants.EGM96_EARTH_MU), mass);

		Vector3D acc = Vector3D.ZERO;
		for (int j = 0; j < forces.length; j++)
		    acc = acc.add(forces[j].acceleration(ss, params.get(j)));

		Xdot[i]   = X[i+3];
		Xdot[i+1] = X[i+4];
		Xdot[i+2] = X[i+5];
		Xdot[i+3] = acc.getX();
		Xdot[i+4] = acc.getY();
		Xdot[i+5] = acc.getZ();
	    }
	}
	catch (Exception exc)
	{
	    except = exc;
	}

	return(Xdot);
    }
}
