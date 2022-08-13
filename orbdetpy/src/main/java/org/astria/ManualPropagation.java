/*
 * ManualPropagation.java - Low level numerical state propagation.
 * Copyright (C) 2018-2022 University of Texas
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

import java.util.concurrent.CountDownLatch;
import java.util.concurrent.TimeUnit;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.hipparchus.ode.ODEIntegrator;
import org.hipparchus.ode.ODEState;
import org.hipparchus.ode.OrdinaryDifferentialEquation;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.orekit.forces.ForceModel;
import org.orekit.frames.Frame;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

public final class ManualPropagation
{
    private final Settings odCfg;
    private final int stateDim;
    private boolean enableDMC;
    private final int[] parPos;

    private final Task[] tasks;
    private CountDownLatch taskSignal;
    private double propStart;
    private double propFinal;

    public ManualPropagation(Settings odCfg)
    {
        this.odCfg = odCfg;
        this.stateDim = odCfg.parameters.size() + 6;
        this.enableDMC = true;
        this.parPos = new int[odCfg.forces.size() + 3];
        for (int i = 0; i < odCfg.forces.size(); i++)
        {
            parPos[i] = -1;
            ForceModel fmod = odCfg.forces.get(i);
            double[] fpar = fmod.getParameters();
            for (int j = 0; j < odCfg.parameters.size(); j++)
            {
                if (fmod.isSupported(odCfg.parameters.get(j).name))
                    parPos[i] = j + 6;
            }
        }

        for (int i = 0; i < odCfg.parameters.size(); i++)
        {
            if (odCfg.parameters.get(i).name.equals(Estimation.DMC_ACC_ESTM[0]))
            {
                parPos[parPos.length - 3] = i + 6;
                parPos[parPos.length - 2] = i + 7;
                parPos[parPos.length - 1] = i + 8;
                break;
            }
        }

        int count = FastMath.min(2*stateDim, Runtime.getRuntime().availableProcessors());
        int q = (2*stateDim)/count, r = (2*stateDim)%count;
        this.tasks = new Task[count];
        for (int i = 0; i < count; i++)
        {
            if (i == 0)
                tasks[i] = new Task((q + r)*stateDim);
            else
                tasks[i] = new Task(q*stateDim);
        }
    }

    public Array2DRowRealMatrix propagate(double t0, Array2DRowRealMatrix inp, double t1, Array2DRowRealMatrix out, boolean enableDMC)
    {
        propStart = t0;
        propFinal = t1;
        this.enableDMC = enableDMC;
        taskSignal = new CountDownLatch(tasks.length);

        for (int i = 0, k = 0; i < tasks.length; i++)
        {
            for (int j = 0; j < tasks[i].xStart.length/stateDim; j++, k++)
                System.arraycopy(inp.getColumn(k), 0, tasks[i].xStart, j*stateDim, stateDim);
            DataManager.threadPool.execute(tasks[i]);
        }

        try
        {
            if (!taskSignal.await(10L, TimeUnit.MINUTES))
                throw(new RuntimeException("UKF sigma point propagation timed out"));
        }
        catch (Exception exc)
        {
            throw(new RuntimeException(exc));
        }

        for (int i = 0, k = 0; i < tasks.length; i++)
        {
            if (tasks[i].xFinal == null)
                throw(new RuntimeException(tasks[i].exception));
            for (int j = 0; j < tasks[i].xFinal.length; j++)
            {
                out.setEntry(j % stateDim, k, tasks[i].xFinal[j]);
                if ((j + 1) % stateDim == 0)
                    k++;
            }
        }
        return(out);
    }

    private final class Task implements OrdinaryDifferentialEquation, Runnable
    {
        private final ODEIntegrator odeInt;
        public double[] xStart;
        public double[] xDeriv;
        public double[] xFinal;
        public Exception exception;

        public Task(int colSize)
        {
            xStart = new double[colSize];
            xDeriv = new double[colSize];
            odeInt = new DormandPrince853Integrator(odCfg.integMinTimeStep, odCfg.integMaxTimeStep, odCfg.integAbsTolerance, odCfg.integRelTolerance);
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
            return(xStart.length);
        }

        @Override public double[] computeDerivatives(double time, double[] X)
        {
            AbsoluteDate timeAd = odCfg.propStart.shiftedBy(time);
            for (int i = 0; i < X.length; i += stateDim)
            {
                CartesianOrbit orb = new CartesianOrbit(new PVCoordinates(new Vector3D(X[i], X[i+1], X[i+2]), new Vector3D(X[i+3], X[i+4], X[i+5])),
                                                        odCfg.propInertialFrame, timeAd, Constants.EGM96_EARTH_MU);
                SpacecraftState scState = new SpacecraftState(orb, odCfg.rsoMass);

                Vector3D acc = Vector3D.ZERO;
                for (int j = 0; j < odCfg.forces.size(); j++)
                {
                    ForceModel model = odCfg.forces.get(j);
                    double[] param = model.getParameters();
                    if (parPos[j] != -1)
                        param[0] = X[i + parPos[j]];
                    acc = acc.add(model.acceleration(scState, param));
                }

                double[] accArray = acc.toArray();
                xDeriv[i] = X[i + 3];
                xDeriv[i + 1] = X[i + 4];
                xDeriv[i + 2] = X[i + 5];
                xDeriv[i + 3] = accArray[0];
                xDeriv[i + 4] = accArray[1];
                xDeriv[i + 5] = accArray[2];
                if (enableDMC && odCfg.estmDMCCorrTime > 0.0 && odCfg.estmDMCSigmaPert > 0.0)
                {
                    xDeriv[i + 3] += X[i + parPos[parPos.length - 3]];
                    xDeriv[i + 4] += X[i + parPos[parPos.length - 2]];
                    xDeriv[i + 5] += X[i + parPos[parPos.length - 1]];
                    xDeriv[i + parPos[parPos.length - 3]] = -X[i + parPos[parPos.length - 3]]/odCfg.estmDMCCorrTime;
                    xDeriv[i + parPos[parPos.length - 2]] = -X[i + parPos[parPos.length - 2]]/odCfg.estmDMCCorrTime;
                    xDeriv[i + parPos[parPos.length - 1]] = -X[i + parPos[parPos.length - 1]]/odCfg.estmDMCCorrTime;
                }
            }
            return(xDeriv);
        }
    }
}
