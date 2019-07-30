/*
 * PropagatorBuilder.java - Wrapper for Orekit's propagator builder.
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
import org.orekit.orbits.Orbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.conversion.NumericalPropagatorBuilder;
import org.orekit.propagation.conversion.ODEIntegratorBuilder;
import org.orekit.propagation.integration.AdditionalEquations;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.ParameterDriversList;

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
	    prop.setInitialState(prop.getInitialState().addAdditionalState(Estimation.DMC_ACC_PROP,
									   plst.findByName(Estimation.DMC_ACC_ESTM+0).getValue(),
									   plst.findByName(Estimation.DMC_ACC_ESTM+1).getValue(),
									   plst.findByName(Estimation.DMC_ACC_ESTM+2).getValue()));
	}
	odcfg.addEventHandlers(prop, prop.getInitialState());

	return(prop);
    }

    class DMCEquations implements AdditionalEquations
    {
	private double[] acceci;

	public DMCEquations()
	{
	    acceci = new double[6];
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
	    double[] acc = sta.getAdditionalState(Estimation.DMC_ACC_PROP);
	    for (int i = 0; i < 3; i++)
	    {
		acceci[i+3] = acc[i];
		pdot[i] = -acc[i]/odcfg.Estimation.DMCCorrTime;
	    }

	    return(acceci);
	}
    }
}
