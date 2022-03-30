/*
 * PropagatorBuilder.java - Wrapper for Orekit's propagator builder.
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

import org.orekit.attitudes.AttitudeProvider;
import org.orekit.forces.ForceModel;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.conversion.DormandPrince853IntegratorBuilder;
import org.orekit.propagation.conversion.NumericalPropagatorBuilder;
import org.orekit.propagation.integration.AdditionalEquations;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.utils.ParameterDriver;
import org.orekit.utils.ParameterDriversList;

public final class PropagatorBuilder extends NumericalPropagatorBuilder
{
    private final Settings odCfg;
    protected boolean enableDMC;

    public PropagatorBuilder(Settings odCfg, Orbit orb, boolean enableDMC)
    {
	super(orb, new DormandPrince853IntegratorBuilder(odCfg.integMinTimeStep, odCfg.integMaxTimeStep, 1.0), PositionAngle.TRUE, 10.0);
	this.odCfg = odCfg;
	this.enableDMC = enableDMC;
	setMass(odCfg.rsoMass);
	for (ForceModel fm: odCfg.forces)
	    addForceModel(fm);

	ParameterDriversList plst = getPropagationParametersDrivers();
	for (Settings.Parameter ep: odCfg.parameters)
	{
	    ParameterDriver pdrv = new ParameterDriver(ep.name, ep.value, 1.0, ep.min, ep.max);
	    pdrv.setReferenceDate(odCfg.propStart);
	    pdrv.setSelected(true);
	    plst.add(pdrv);
	}

	AttitudeProvider attProv = odCfg.getAttitudeProvider();
	if (attProv != null)
	    setAttitudeProvider(attProv);
	if (odCfg.estmDMCCorrTime > 0.0 && odCfg.estmDMCSigmaPert > 0.0)
	    addAdditionalEquations(new DMCDerivatives());
    }

    @Override public NumericalPropagator buildPropagator(double[] par)
    {
	NumericalPropagator prop = super.buildPropagator(par);
	if (odCfg.estmDMCCorrTime > 0.0 && odCfg.estmDMCSigmaPert > 0.0)
	{
	    ParameterDriversList plst = getPropagationParametersDrivers();
	    prop.setInitialState(prop.getInitialState().addAdditionalState(Estimation.DMC_ACC_PROP, plst.findByName(Estimation.DMC_ACC_ESTM[0]).getValue(),
									   plst.findByName(Estimation.DMC_ACC_ESTM[1]).getValue(),
									   plst.findByName(Estimation.DMC_ACC_ESTM[2]).getValue()));
	}
	return(prop);
    }

    private class DMCDerivatives implements AdditionalEquations
    {
	@Override public double[] computeDerivatives(SpacecraftState state, double[] pdot)
	{
	    double[] accEci = new double[6];
            if (enableDMC)
            {
                double[] acc = state.getAdditionalState(Estimation.DMC_ACC_PROP);
                for (int i = 0; i < 3; i++)
		{
		    accEci[i+3] = acc[i];
                    pdot[i] = -acc[i]/odCfg.estmDMCCorrTime;
		}
            }
	    else
	    {
		for (int i = 0; i < 3; i++)
		    pdot[i] = 0.0;
	    }
	    return(accEci);
	}

	@Override public String getName()
	{
	    return(Estimation.DMC_ACC_PROP);
	}
    }
}
