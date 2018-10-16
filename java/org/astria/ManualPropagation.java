/*
 * ManualPropagation.java - Low level numerical state propagation.
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

import java.util.Arrays;
import org.astria.DataManager;
import org.astria.Settings;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.orekit.forces.ForceModel;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.Constants;
import org.orekit.utils.ParameterDriver;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public class ManualPropagation implements OrdinaryDifferentialEquation
{
    protected Settings odcfg;
    protected int intvecdim;
    protected int statedim;
    protected AbsoluteDate epoch;
    protected double[] Xdot;
    protected Array2DRowRealMatrix ecirot;

    protected ODEIntegrator odeint;

    protected Exception except;

    public ManualPropagation(Settings cfg, int vecdim)
    {
	odcfg = cfg;
	intvecdim = vecdim;
	statedim = odcfg.estparams.size() + 6;

	epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.Start),
				 DataManager.utcscale);
	Xdot = new double[intvecdim];
	ecirot = new Array2DRowRealMatrix(3, 3);

	odeint = new DormandPrince853Integrator(1E-3, 300.0, 1E-14, 1E-12);
    }

    public double[] propagate(double t0, double[] X0, double t1) throws Exception
    {
	except = null;
	double[] retval = odeint.integrate(this, new ODEState(t0, X0), t1).getPrimaryState();

	if (except != null)
	    throw(except);
	return(retval);
    }

    public int getDimension()
    {
	return(intvecdim);
    }

    public double[] computeDerivatives(double t, double[] X)
    {
	try
	{
	    int i,k,l;
	    AbsoluteDate tm = new AbsoluteDate(epoch, t);

	    for (i = 0; i < X.length; i += statedim)
	    {
		SpacecraftState ss = new SpacecraftState(
		    new CartesianOrbit(new PVCoordinates(
					   new Vector3D(X[i],   X[i+1], X[i+2]),
					   new Vector3D(X[i+3], X[i+4], X[i+5])),
				       DataManager.eme2000, tm, Constants.EGM96_EARTH_MU),
		    odcfg.SpaceObject.Mass);

		l = 6;
		Vector3D acc = Vector3D.ZERO;
		for (ForceModel fmod : odcfg.forces)
		{
		    k = 0;
		    double[] fpar = fmod.getParameters();
		    for (ParameterDriver drv : fmod.getParametersDrivers())
		    {
			for (Settings.EstimatedParameter emp : odcfg.estparams)
			{
			    if (drv.getName().equals(emp.name))
				fpar[k] = X[i + l++];
			}
			k++;
		    }

		    acc = acc.add(fmod.acceleration(ss, fpar));
		}

		Xdot[i]   = X[i+3];
		Xdot[i+1] = X[i+4];
		Xdot[i+2] = X[i+5];
		Xdot[i+3] = acc.getX();
		Xdot[i+4] = acc.getY();
		Xdot[i+5] = acc.getZ();

		if (odcfg.Estimation.DMCCorrTime > 0.0 && odcfg.Estimation.DMCSigmaPert > 0.0)
		{
		    TimeStampedPVCoordinates pvc = ss.getPVCoordinates();
		    Vector3D r = pvc.getPosition().normalize();
		    Vector3D h = pvc.getMomentum().normalize();
		    ecirot.setColumn(0, r.toArray());
		    ecirot.setColumn(1, h.crossProduct(r).toArray());
		    ecirot.setColumn(2, h.toArray());
		    double[] rot = ecirot.operate(Arrays.copyOfRange(X, i+statedim-3, i+statedim));

		    Xdot[i+3] += rot[0];
		    Xdot[i+4] += rot[1];
		    Xdot[i+5] += rot[2];
		    Xdot[i+statedim-3] = -X[i+statedim-3]/odcfg.Estimation.DMCCorrTime;
		    Xdot[i+statedim-2] = -X[i+statedim-2]/odcfg.Estimation.DMCCorrTime;
		    Xdot[i+statedim-1] = -X[i+statedim-1]/odcfg.Estimation.DMCCorrTime;
		}
	    }
	}
	catch (Exception exc)
	{
	    except = exc;
	}

	return(Xdot);
    }
}
