/*
 * ParallelPropagation.java - Multiple object propagation functions.
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

import java.util.ArrayList;
import java.util.Map;
import java.util.Random;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.forces.ForceModel;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.measurements.RangeRate;
import org.orekit.estimation.measurements.modifiers.Bias;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.propagation.sampling.OrekitStepInterpolator;
import org.orekit.propagation.sampling.MultiSatStepHandler;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.AbsolutePVCoordinates;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class ParallelPropagation
{
    private final ArrayList<Settings> configObjs;
    private final ArrayList<Propagator> propagators;

    public ParallelPropagation(ArrayList<Settings> configObjs)
    {
	this.configObjs = configObjs;
	this.propagators = new ArrayList<Propagator>(configObjs.size());
	for (int i = 0; i < configObjs.size(); i++)
	    propagators.add(buildPropagator(configObjs.get(i)));
    }

    public ArrayList<ArrayList<Measurements.Measurement>> propagate() throws Exception
    {
	Settings obj0 = configObjs.get(0);
	AbsoluteDate start = obj0.propStart, end = obj0.propEnd;
	double step = obj0.propStep;
	for (Settings s: configObjs)
	{
	    if (obj0.propStep > 0.0)
	    {
		if (start.durationFrom(s.propStart) > 0.0)
		    start = s.propStart;
		if (end.durationFrom(s.propEnd) < 0.0)
		    end = s.propEnd;
		if (s.propStep != 0.0)
		    step = FastMath.min(step, s.propStep);
	    }
	    else
	    {
		if (start.durationFrom(s.propStart) < 0.0)
		    start = s.propStart;
		if (end.durationFrom(s.propEnd) > 0.0)
		    end = s.propEnd;
		if (s.propStep != 0.0)
		    step = FastMath.max(step, s.propStep);
	    }
	}

	AbsoluteDate tmFinal = start;
	AbsoluteDate tmStart = new AbsoluteDate(tmFinal, -0.1);
	int capacity = (int)FastMath.abs(end.durationFrom(start)/step) + 2;
	ArrayList<ArrayList<Measurements.Measurement>> output = new ArrayList<ArrayList<Measurements.Measurement>>(propagators.size());
	for (int i = 0; i < propagators.size(); i++)
	    output.add(new ArrayList<Measurements.Measurement>(capacity));

	while (true)
	{
	    AbsoluteDate t0 = tmStart;
	    AbsoluteDate t1 = tmFinal;
	    for (int i = 0; i < propagators.size(); i++)
	    {
		Propagator prop = propagators.get(i);
		SpacecraftState state = DataManager.threadPool.submit(()->prop.propagate(t0, t1)).get();
		tmStart = state.getDate();
		Simulation.simulate(configObjs.get(i), state, output.get(i));
	    }

	    double dt = end.durationFrom(tmStart);
	    if (step >= 0.0)
		tmFinal = new AbsoluteDate(tmFinal, FastMath.min(dt, step));
	    else
	    {
		tmFinal = new AbsoluteDate(tmFinal, FastMath.max(dt, step));
		dt = -dt;
	    }
	    if (dt <= 0.0)
		break;
	}

	return(output);
    }

    private Propagator buildPropagator(Settings cfg)
    {
	Propagator prop;
	if (cfg.propInitialTLE != null && cfg.propInitialTLE.length == 2)
	{
	    TLE parser = new TLE(cfg.propInitialTLE[0], cfg.propInitialTLE[1]);
	    prop = TLEPropagator.selectExtrapolator(parser);
	    if (cfg.propStart == null)
		cfg.propStart = parser.getDate();
	    if (cfg.propEnd == null)
		cfg.propEnd = cfg.propStart;
	}
	else
	{
	    double[] X = cfg.getInitialState();
	    NumericalPropagator np = new NumericalPropagator(
		new DormandPrince853Integrator(cfg.integMinTimeStep, cfg.integMaxTimeStep, cfg.integAbsTolerance, cfg.integRelTolerance));
	    for (ForceModel fm: cfg.forces)
		np.addForceModel(fm);
	    np.setInitialState(new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(X[0], X[1], X[2]), new Vector3D(X[3], X[4], X[5])),
								      cfg.propInertialFrame, cfg.propStart, Constants.EGM96_EARTH_MU), cfg.rsoMass));
	    prop = np;
	}

	AttitudeProvider attProv = cfg.getAttitudeProvider();
	if (attProv != null)
	    prop.setAttitudeProvider(attProv);
	cfg.addEventHandlers(prop);
	return(prop);
    }

    private static class Simulation
    {
	private static final Random rand = new Random(1);
	private static final double[] oneOnes = new double[]{1.0};
	private static final double[] oneNegInf = new double[]{Double.NEGATIVE_INFINITY};
	private static final double[] onePosInf = new double[]{Double.POSITIVE_INFINITY};
	private static final double[] twoZeros = new double[]{0.0, 0.0};
	private static final double[] twoOnes = new double[]{1.0, 1.0};
	private static final double[] twoNegInf = new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
	private static final double[] twoPosInf = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
	private static final String[] biasAzEl = new String[]{Settings.MeasurementType.AZIMUTH.name(), Settings.MeasurementType.ELEVATION.name()};
	private static final String[] biasRaDec = new String[]{Settings.MeasurementType.RIGHT_ASCENSION.name(),
							       Settings.MeasurementType.DECLINATION.name()};
	private static final String[] biasRange = new String[]{Settings.MeasurementType.RANGE.name()};
	private static final String[] biasRangeRate = new String[]{Settings.MeasurementType.RANGE_RATE.name()};
	private static final ObservableSatellite obsSat = new ObservableSatellite(0);
	private static final SpacecraftState[] ssStates = new SpacecraftState[1];

	public static void simulate(Settings simCfg, SpacecraftState state, ArrayList<Measurements.Measurement> measOut)
	{
	    TimeStampedPVCoordinates pvc = state.getPVCoordinates(simCfg.propInertialFrame);
	    Measurements.Measurement meas = new Measurements.Measurement(pvc);
	    if (!simCfg.simMeasurements)
	    {
		measOut.add(meas);
		return;
	    }

	    ssStates[0] = state;
	    for (Map.Entry<String, GroundStation> kv: simCfg.stations.entrySet())
	    {
		GroundStation gst = kv.getValue();
		AngularAzEl azel = new AngularAzEl(gst, pvc.getDate(), twoZeros, twoZeros, twoOnes, obsSat);
		if (azel.estimate(0, 0, ssStates).getEstimatedValue()[1] <= 5E-6)
		    continue;

		Measurements.Measurement clone = new Measurements.Measurement(meas);
		clone.station = kv.getKey();

		double[] bias;
		Settings.Station jsn = simCfg.cfgStations.get(clone.station);
		if (jsn.bias != null && jsn.bias.length > 0)
		    bias = jsn.bias;
		else
		    bias = twoZeros;

		for (Map.Entry<Settings.MeasurementType, Settings.Measurement> nvp: simCfg.cfgMeasurements.entrySet())
		{
		    Settings.MeasurementType name = nvp.getKey();
		    Settings.Measurement val = nvp.getValue();
		    if (name == Settings.MeasurementType.RANGE)
		    {
			Range obs = new Range(gst, val.twoWay, pvc.getDate(), 0.0, 0.0, 1.0, obsSat);
			obs.addModifier(new Bias<Range>(biasRange, new double[]{rand.nextGaussian()*val.error[0]+bias[0]},
							oneOnes, oneNegInf, onePosInf));
			if (clone.values == null)
			    clone.values = new double[2];
			clone.values[0] = obs.estimate(0, 0, ssStates).getEstimatedValue()[0];
		    }
		    else if (name == Settings.MeasurementType.RANGE_RATE)
		    {
			RangeRate obs = new RangeRate(gst, pvc.getDate(), 0.0, 0.0, 1.0, val.twoWay, obsSat);
			obs.addModifier(new Bias<RangeRate>(biasRangeRate, new double[]{rand.nextGaussian()*val.error[0]+bias[1]},
							    oneOnes, oneNegInf, onePosInf));
			if (clone.values == null)
			    clone.values = new double[2];
			clone.values[1] = obs.estimate(0, 0, ssStates).getEstimatedValue()[0];
		    }
		    else if (name == Settings.MeasurementType.RIGHT_ASCENSION || name == Settings.MeasurementType.DECLINATION &&
			     clone.values == null)
		    {
			AngularRaDec obs = new AngularRaDec(gst, simCfg.propInertialFrame, pvc.getDate(), twoZeros, twoZeros, twoOnes, obsSat);
			obs.addModifier(new Bias<AngularRaDec>(biasRaDec, new double[]{rand.nextGaussian()*val.error[0]+bias[0],
										       rand.nextGaussian()*val.error[0]+bias[1]},
				twoOnes, twoNegInf, twoPosInf));
			clone.values = obs.estimate(0, 0, ssStates).getEstimatedValue();
		    }
		    else if (name == Settings.MeasurementType.AZIMUTH || name == Settings.MeasurementType.ELEVATION &&
			     clone.values == null)
		    {
			azel.addModifier(new Bias<AngularAzEl>(biasAzEl, new double[]{rand.nextGaussian()*val.error[0]+bias[0],
										      rand.nextGaussian()*val.error[0]+bias[1]},
				twoOnes, twoNegInf, twoPosInf));
			clone.values = azel.estimate(0, 0, ssStates).getEstimatedValue();
		    }
		}

		measOut.add(clone);
	    }
	}
    }
}
