/*
 * ManualPropagation.java - Low level numerical state propagation.
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

import java.util.Arrays;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.orekit.attitudes.Attitude;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.attitudes.LofOffset;
import org.orekit.forces.ForceModel;
import org.orekit.frames.Frame;
import org.orekit.frames.LOFType;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.Constants;
import org.orekit.utils.ParameterDriver;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.PVCoordinatesProvider;
import org.orekit.utils.TimeStampedPVCoordinates;

public class ManualPropagation implements OrdinaryDifferentialEquation, PVCoordinatesProvider
{
    protected Settings odcfg;
    protected AttitudeProvider attprov;
    protected int intvecdim;
    protected int statedim;

    protected double[] Xdot;
    protected AbsoluteDate epoch;

    protected ODEIntegrator odeint;

    protected LofOffset loframe;
    protected AbsoluteDate savedtime;
    protected TimeStampedPVCoordinates savedposv;

    public ManualPropagation(Settings cfg, int vecdim)
    {
	odcfg = cfg;
	intvecdim = vecdim;
	statedim = odcfg.estparams.size() + 6;
	attprov = cfg.getAttitudeProvider();
	loframe = new LofOffset(odcfg.propframe, LOFType.VVLH);
	Xdot = new double[intvecdim];
	epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(odcfg.Propagation.Start),
				 DataManager.utcscale);
	odeint = new DormandPrince853Integrator(cfg.Integration.MinTimeStep, cfg.Integration.MaxTimeStep,
						cfg.Integration.AbsTolerance, cfg.Integration.RelTolerance);
    }

    public double[] propagate(double t0, double[] X0, double t1)
    {
	return(odeint.integrate(this, new ODEState(t0, X0), t1).getPrimaryState());
    }

    public Attitude getAttitude(AbsoluteDate time, double[] X)
    {
	savedtime = new AbsoluteDate(time, 0.0);
	savedposv = new TimeStampedPVCoordinates(time, new Vector3D(X[0], X[1], X[2]),
						 new Vector3D(X[3], X[4], X[5]));

	if (attprov != null)
	    return(attprov.getAttitude(this, time, odcfg.propframe));
	else
	    return(loframe.getAttitude(this, time, odcfg.propframe));
    }

    public int getDimension()
    {
	return(intvecdim);
    }

    public double[] computeDerivatives(double t, double[] X)
    {
	Arrays.fill(Xdot, 0.0);
	AbsoluteDate tm = new AbsoluteDate(epoch, t);
	for (int i = 0; i < X.length; i += statedim)
	{
	    SpacecraftState ss = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(X[i],   X[i+1], X[i+2]),
											  new Vector3D(X[i+3], X[i+4], X[i+5])),
									odcfg.propframe, tm, Constants.EGM96_EARTH_MU),
						     getAttitude(tm, Arrays.copyOfRange(X, i, i+6)), odcfg.SpaceObject.Mass);

	    Vector3D acc = Vector3D.ZERO;
	    for (ForceModel fmod : odcfg.forces)
	    {
		double[] fpar = fmod.getParameters();
		if (X.length > statedim)
		{
		    for (int j = 0; j < odcfg.estparams.size(); j++)
		    {
			Settings.EstimatedParameter emp = odcfg.estparams.get(j);
			if (fmod.isSupported(emp.name))
			    fpar[0] = X[i + j + 6];
		    }
		}
		acc = acc.add(fmod.acceleration(ss, fpar));
	    }

	    Xdot[i]   = X[i+3];
	    Xdot[i+1] = X[i+4];
	    Xdot[i+2] = X[i+5];
	    Xdot[i+3] = acc.getX();
	    Xdot[i+4] = acc.getY();
	    Xdot[i+5] = acc.getZ();
	    if (X.length > statedim && odcfg.Estimation.DMCCorrTime > 0.0 && odcfg.Estimation.DMCSigmaPert > 0.0)
	    {
		Xdot[i+3] += X[i+statedim-3];
		Xdot[i+4] += X[i+statedim-2];
		Xdot[i+5] += X[i+statedim-1];
		Xdot[i+statedim-3] = -X[i+statedim-3]/odcfg.Estimation.DMCCorrTime;
		Xdot[i+statedim-2] = -X[i+statedim-2]/odcfg.Estimation.DMCCorrTime;
		Xdot[i+statedim-1] = -X[i+statedim-1]/odcfg.Estimation.DMCCorrTime;
	    }
	}

	return(Xdot);
    }

    public TimeStampedPVCoordinates getPVCoordinates(AbsoluteDate date, Frame frame)
    {
	return(odcfg.propframe.getTransformTo(frame, date).
	       transformPVCoordinates(savedposv.shiftedBy(date.durationFrom(savedtime))));
    }
}
