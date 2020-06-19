/*
 * ManualPropagation.java - Low level numerical state propagation.
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

import java.util.Arrays;
import java.util.Map;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
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
import org.orekit.utils.Constants;
import org.orekit.utils.ParameterDriver;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.PVCoordinatesProvider;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class ManualPropagation implements PVCoordinatesProvider
{
    private final Settings odCfg;
    private final int stateDim;

    private final Task[] tasks;
    private CountDownLatch taskSignal;
    private double propStart;
    private double propFinal;

    private final AttitudeProvider attProvider;
    private final LofOffset lofFrame;
    private AbsoluteDate savedTime;
    private TimeStampedPVCoordinates savedPos;

    private boolean enableDMC;
    private final int[] posParameters;

    public ManualPropagation(Settings odCfg)
    {
	this.odCfg = odCfg;
	this.stateDim = odCfg.parameters.size() + 6;

	this.tasks = new Task[2*stateDim];
	for (int i = 0; i < tasks.length; i++)
	    tasks[i] = new Task();

	attProvider = odCfg.getAttitudeProvider();
	lofFrame = new LofOffset(odCfg.propInertialFrame, LOFType.VVLH);

	enableDMC = true;
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
	    if (odCfg.parameters.get(i).name.equalsIgnoreCase(Estimation.DMC_ACC_ESTM[0]))
	    {
		posParameters[posParameters.length - 3] = i + 6;
		posParameters[posParameters.length - 2] = i + 7;
		posParameters[posParameters.length - 1] = i + 8;
		break;
	    }
	}
    }

    public Array2DRowRealMatrix propagate(double t0, Array2DRowRealMatrix inp,
					  double t1, Array2DRowRealMatrix out, boolean enableDMC)
    {
	propStart = t0;
	propFinal = t1;
	this.enableDMC = enableDMC;
	taskSignal = new CountDownLatch(inp.getColumnDimension());
	for (int i = 0; i < inp.getColumnDimension(); i++)
	{
	    tasks[i].xStart = inp.getColumn(i);
	    DataManager.threadPool.execute(tasks[i]);
	}

	try
	{
	    if (!taskSignal.await(10L, TimeUnit.MINUTES))
		throw(new RuntimeException("Propagation of UKF sigma points timed out"));
	}
	catch (Exception exc)
	{
	    throw(new RuntimeException(exc));
	}

	for (int i = 0; i < inp.getColumnDimension(); i++)
	{
	    if (tasks[i].xFinal != null)
		out.setColumn(i, tasks[i].xFinal);
	    else
		throw(new RuntimeException(tasks[i].exception));
	}
	return(out);
    }

    public Attitude getAttitude(AbsoluteDate time, double[] X)
    {
	savedTime = new AbsoluteDate(time, 0.0);
	savedPos = new TimeStampedPVCoordinates(time, new Vector3D(X[0], X[1], X[2]),
						new Vector3D(X[3], X[4], X[5]));
	if (attProvider != null)
	    return(attProvider.getAttitude(this, time, odCfg.propInertialFrame));
	else
	    return(lofFrame.getAttitude(this, time, odCfg.propInertialFrame));
    }

    @Override public TimeStampedPVCoordinates getPVCoordinates(AbsoluteDate date, Frame frame)
    {
	return(odCfg.propInertialFrame.getTransformTo(frame, date).
	       transformPVCoordinates(savedPos.shiftedBy(date.durationFrom(savedTime))));
    }

    private final class Task implements Runnable, OrdinaryDifferentialEquation
    {
	private final ODEIntegrator odeInt;
	public double[] xStart;
	public double[] xFinal;
	public Exception exception;

	public Task()
	{
	    odeInt = new DormandPrince853Integrator(odCfg.integMinTimeStep, odCfg.integMaxTimeStep,
						    odCfg.integAbsTolerance, odCfg.integRelTolerance);
	}

	@Override public void run()
	{
	    try
	    {
		exception = null;
		xFinal = odeInt.integrate(this, new ODEState(propStart, xStart), propFinal).getPrimaryState();
	    }
	    catch (Exception exc)
	    {
		xFinal = null;
		exception = exc;
	    }
	    finally
	    {
		taskSignal.countDown();
	    }
	}

	@Override public int getDimension()
	{
	    return(stateDim);
	}

	@Override public double[] computeDerivatives(double t, double[] X)
	{
	    final double[] xDot = new double[stateDim];
	    Arrays.fill(xDot, 0.0);

	    final AbsoluteDate tm = new AbsoluteDate(odCfg.propStart, t);
	    final SpacecraftState ss = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(X[0], X[1], X[2]),
												new Vector3D(X[3], X[4], X[5])),
									      odCfg.propInertialFrame, tm, Constants.EGM96_EARTH_MU),
							   getAttitude(tm, Arrays.copyOfRange(X, 0, 6)), odCfg.rsoMass);

	    Vector3D acc = Vector3D.ZERO;
	    for (int j = 0; j < odCfg.forces.size(); j++)
	    {
		ForceModel fmod = odCfg.forces.get(j);
		double[] fpar = fmod.getParameters();
		if (posParameters[j] != -1)
		    fpar[0] = X[posParameters[j]];
		acc = acc.add(fmod.acceleration(ss, fpar));
	    }

	    xDot[0] = X[3];
	    xDot[1] = X[4];
	    xDot[2] = X[5];
	    xDot[3] = acc.getX();
	    xDot[4] = acc.getY();
	    xDot[5] = acc.getZ();
	    if (enableDMC && odCfg.estmDMCCorrTime > 0.0 && odCfg.estmDMCSigmaPert > 0.0)
	    {
		xDot[3] += X[posParameters[posParameters.length - 3]];
		xDot[4] += X[posParameters[posParameters.length - 2]];
		xDot[5] += X[posParameters[posParameters.length - 1]];
		xDot[posParameters[posParameters.length - 3]] = -X[posParameters[posParameters.length - 3]]/odCfg.estmDMCCorrTime;
		xDot[posParameters[posParameters.length - 2]] = -X[posParameters[posParameters.length - 2]]/odCfg.estmDMCCorrTime;
		xDot[posParameters[posParameters.length - 1]] = -X[posParameters[posParameters.length - 1]]/odCfg.estmDMCCorrTime;
	    }

	    return(xDot);
	}
    }
}
