/*
 * PropagatorBuilder.java - Wrapper for Orekit's propagator builder.
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
import org.astria.Estimation;
import org.astria.Settings;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.linear.Array2DRowRealMatrix;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.conversion.NumericalPropagatorBuilder;
import org.orekit.propagation.conversion.ODEIntegratorBuilder;
import org.orekit.propagation.integration.AdditionalEquations;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.ParameterDriversList;
import org.orekit.utils.TimeStampedPVCoordinates;

public class PropagatorBuilder extends NumericalPropagatorBuilder
{
    protected Settings odcfg;

    public PropagatorBuilder(Settings cfg, Orbit orb, ODEIntegratorBuilder ode,
			     PositionAngle ang, double pos)
    {
	super(orb, ode, ang, pos);
	odcfg = cfg;
    }

    @Override public NumericalPropagator buildPropagator(double[] par)
    {
	NumericalPropagator prop = super.buildPropagator(par);
	if (odcfg.Estimation.DMCCorrTime > 0.0 && odcfg.Estimation.DMCSigmaPert > 0.0)
	{
	    prop.addAdditionalEquations(new DMCEquations());
	    ParameterDriversList plst = getPropagationParametersDrivers();
	    prop.setInitialState(
		prop.getInitialState().addAdditionalState(Estimation.DMC_ACC_PROP,
							  plst.findByName(Estimation.DMC_ACC_ESTM+0).getValue(),
							  plst.findByName(Estimation.DMC_ACC_ESTM+1).getValue(),
							  plst.findByName(Estimation.DMC_ACC_ESTM+2).getValue()));
	}

	return(prop);
    }

    class DMCEquations implements AdditionalEquations
    {
	private double[] acceci;
	private Array2DRowRealMatrix ecirot;

	public DMCEquations()
	{
	    acceci = new double[6];
	    ecirot = new Array2DRowRealMatrix(3, 3);
	}

	public void init(SpacecraftState sta, AbsoluteDate tgt)
	{
	}

	public String getName()
	{
	    return(Estimation.DMC_ACC_PROP);
	}

	public double[] computeDerivatives(SpacecraftState sta, double[] pdot)
	{
	    double[] accric = sta.getAdditionalState(Estimation.DMC_ACC_PROP);
	    for (int i = 0; i < 3; i++)
		pdot[i] = -accric[i]/odcfg.Estimation.DMCCorrTime;

	    TimeStampedPVCoordinates pvc = sta.getPVCoordinates();
	    Vector3D r = pvc.getPosition().normalize();
	    Vector3D h = pvc.getMomentum().normalize();
	    ecirot.setColumn(0, r.toArray());
	    ecirot.setColumn(1, h.crossProduct(r).toArray());
	    ecirot.setColumn(2, h.toArray());

	    double[] rot = ecirot.operate(accric);
	    acceci[3] = rot[0];
	    acceci[4] = rot[1];
	    acceci[5] = rot[2];
	    return(acceci);
	}
    }
}
