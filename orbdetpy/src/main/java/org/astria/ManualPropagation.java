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
import java.util.Map;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.ODEStateAndDerivative;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.ode.sampling.StepNormalizer;
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

public final class ManualPropagation implements OrdinaryDifferentialEquation, PVCoordinatesProvider
{
    private final Settings odCfg;
    private final int vectorDim;
    private final int stateDim;

    private final double[] xDot;
    private final AbsoluteDate epoch;
    private final ODEIntegrator odeInt;
    private boolean enableDMC;
    private int[] posParameters;

    private final AttitudeProvider attProvider;
    private final LofOffset lofFrame;
    private AbsoluteDate savedTime;
    private TimeStampedPVCoordinates savedPos;

    public ManualPropagation(Settings cfg, int vecdim, StepNormalizer handler)
    {
	odCfg = cfg;
	vectorDim = vecdim;
	stateDim = odCfg.parameters.size() + 6;
	attProvider = cfg.getAttitudeProvider();

	enableDMC = true;
	xDot = new double[vectorDim];
	epoch = new AbsoluteDate(DateTimeComponents.parseDateTime(cfg.propStart), DataManager.getTimeScale("UTC"));

	posParameters = new int[odCfg.forces.size() + 3];
	for (int i = 0; i < odCfg.forces.size(); i++)
	{
	    posParameters[i] = -1;
	    ForceModel fmod = odCfg.forces.get(i);
	    double[] fpar = fmod.getParameters();
	    for (int j = 0; j < odCfg.parameters.size(); j++)
	    {
		if (fmod.isSupported(odCfg.parameters.get(j).name))
		    posParameters[i] = j + 6;
	    }
	}

	for (int i = 0; i < odCfg.parameters.size(); i++)
	{
	    if (odCfg.parameters.get(i).name.equals(Estimation.DMC_ACC_ESTM[0]))
	    {
		posParameters[posParameters.length - 3] = i + 6;
		posParameters[posParameters.length - 2] = i + 7;
		posParameters[posParameters.length - 1] = i + 8;
		break;
	    }
	}

	odeInt = new DormandPrince853Integrator(cfg.integMinTimeStep, cfg.integMaxTimeStep,
						cfg.integAbsTolerance, cfg.integRelTolerance);
	if (handler != null)
	    odeInt.addStepHandler(handler);

	lofFrame = new LofOffset(odCfg.propFrame, LOFType.VVLH);
    }

    public ODEStateAndDerivative propagate(double t0, double[] X0, double t1)
    {
	return(odeInt.integrate(this, new ODEState(t0, X0), t1));
    }

    public Attitude getAttitude(AbsoluteDate time, double[] X)
    {
	savedTime = new AbsoluteDate(time, 0.0);
	savedPos = new TimeStampedPVCoordinates(time, new Vector3D(X[0], X[1], X[2]),
						new Vector3D(X[3], X[4], X[5]));
	if (attProvider != null)
	    return(attProvider.getAttitude(this, time, odCfg.propFrame));
	else
	    return(lofFrame.getAttitude(this, time, odCfg.propFrame));
    }

    public void setDMCState(boolean enabled)
    {
	enableDMC = enabled;
    }

    @Override public int getDimension()
    {
	return(vectorDim);
    }

    @Override public double[] computeDerivatives(double t, double[] X)
    {
	SpacecraftState ss;
	Vector3D acc;
	ForceModel fmod;
	double[] fpar;

	Arrays.fill(xDot, 0.0);
	AbsoluteDate tm = new AbsoluteDate(epoch, t);
	for (int i = 0; i < X.length; i += stateDim)
	{
	    ss = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(X[i],   X[i+1], X[i+2]),
									  new Vector3D(X[i+3], X[i+4], X[i+5])),
							odCfg.propFrame, tm, Constants.EGM96_EARTH_MU),
				     getAttitude(tm, Arrays.copyOfRange(X, i, i+6)), odCfg.rsoMass);

	    acc = Vector3D.ZERO;
	    for (int j = 0; j < odCfg.forces.size(); j++)
	    {
		fmod = odCfg.forces.get(j);
		fpar = fmod.getParameters();
		if (posParameters[j] != -1)
		    fpar[0] = X[i + posParameters[j]];
		acc = acc.add(fmod.acceleration(ss, fpar));
	    }

	    xDot[i]   = X[i+3];
	    xDot[i+1] = X[i+4];
	    xDot[i+2] = X[i+5];
	    xDot[i+3] = acc.getX();
	    xDot[i+4] = acc.getY();
	    xDot[i+5] = acc.getZ();
	    if (enableDMC && X.length > stateDim && odCfg.estmDMCCorrTime > 0.0 && odCfg.estmDMCSigmaPert > 0.0)
	    {
		xDot[i + 3] += X[i + posParameters[posParameters.length - 3]];
		xDot[i + 4] += X[i + posParameters[posParameters.length - 2]];
		xDot[i + 5] += X[i + posParameters[posParameters.length - 1]];
		xDot[i + posParameters[posParameters.length - 3]] = -X[i + posParameters[posParameters.length - 3]]/odCfg.estmDMCCorrTime;
		xDot[i + posParameters[posParameters.length - 2]] = -X[i + posParameters[posParameters.length - 2]]/odCfg.estmDMCCorrTime;
		xDot[i + posParameters[posParameters.length - 1]] = -X[i + posParameters[posParameters.length - 1]]/odCfg.estmDMCCorrTime;
	    }
	}

	return(xDot);
    }

    @Override public TimeStampedPVCoordinates getPVCoordinates(AbsoluteDate date, Frame frame)
    {
	return(odCfg.propFrame.getTransformTo(frame, date).
	       transformPVCoordinates(savedPos.shiftedBy(date.durationFrom(savedTime))));
    }
}
