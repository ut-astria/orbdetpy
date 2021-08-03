/*
 * ParallelPropagation.java - Multiple object propagation functions.
 * Copyright (C) 2018-2021 University of Texas
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
import java.util.concurrent.Future;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.nonstiff.DormandPrince853Integrator;
import org.hipparchus.util.FastMath;
import org.orekit.attitudes.AttitudeProvider;
import org.orekit.forces.ForceModel;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.measurements.RangeRate;
import org.orekit.estimation.measurements.modifiers.AngularRadioRefractionModifier;
import org.orekit.estimation.measurements.modifiers.Bias;
import org.orekit.models.earth.EarthITU453AtmosphereRefraction;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.Propagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.tle.TLE;
import org.orekit.propagation.analytical.tle.TLEPropagator;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class ParallelPropagation
{
    private ArrayList<Settings> configObjs;
    private ArrayList<Propagator> propagators;
    private ArrayList<ArrayList<EventHandling>> eventHandlers;
    private ArrayList<AbsoluteDate> stepEnd;

    public ParallelPropagation(ArrayList<Settings> configObjs)
    {
	this.configObjs = configObjs;
    }

    public ArrayList<ArrayList<Measurements.Measurement>> propagate() throws Exception
    {
	propagators = new ArrayList<Propagator>(configObjs.size());
	eventHandlers = new ArrayList<ArrayList<EventHandling>>(configObjs.size());
	stepEnd = new ArrayList<AbsoluteDate>(configObjs.size());
	ArrayList<Future<SpacecraftState>> futures = new ArrayList<Future<SpacecraftState>>(configObjs.size());
	ArrayList<ArrayList<Measurements.Measurement>> output = new ArrayList<ArrayList<Measurements.Measurement>>(configObjs.size());
	for (Settings cfg: configObjs)
	{
	    buildPropagator(cfg);
	    futures.add(null);
	    output.add(new ArrayList<Measurements.Measurement>());
	}

	boolean allDone = false;
	while (!allDone)
	{
	    allDone = true;
	    for (int i = 0; i < configObjs.size(); i++)
	    {
		AbsoluteDate end = stepEnd.get(i);
		if (end != null)
		{
		    allDone = false;
		    Propagator prop = propagators.get(i);
		    futures.set(i, DataManager.threadPool.submit(()->prop.propagate(end)));
		}
		else
		    futures.set(i, null);
	    }

	    for (int i = 0; !allDone && i < configObjs.size(); i++)
	    {
		Future<SpacecraftState> future = futures.get(i);
		if (future == null)
		    continue;

		Settings cfg = configObjs.get(i);
		SpacecraftState state = future.get();
		if (Simulation.simulate(cfg, state, eventHandlers.get(i), output.get(i)) ||
		    (stepEnd.get(i) != null && stepEnd.get(i).equals(cfg.propEnd)))
		{
		    double dt = cfg.propEnd.durationFrom(state.getDate());
		    if (dt != 0.0)
		    {
			if (cfg.propStep > 0.0)
			    stepEnd.set(i, state.getDate().shiftedBy(FastMath.min(dt, cfg.propStep)));
			else
			    stepEnd.set(i, state.getDate().shiftedBy(FastMath.max(dt, cfg.propStep)));
		    }
		    else
			stepEnd.set(i, null);
		}
		else
		    stepEnd.set(i, cfg.propEnd);
	    }
	}
	return(output);
    }

    private void buildPropagator(Settings cfg)
    {
	Propagator prop;
	SpacecraftState state;
	if (cfg.propInitialTLE != null && cfg.propInitialTLE.length == 2)
	{
	    TLE parser = new TLE(cfg.propInitialTLE[0], cfg.propInitialTLE[1]);
	    prop = TLEPropagator.selectExtrapolator(parser);
	    if (cfg.propStart == null)
		cfg.propStart = parser.getDate();
	    if (cfg.propEnd == null)
		cfg.propEnd = cfg.propStart;
	    state = prop.propagate(cfg.propStart);
	}
	else
	{
	    double[] X = cfg.getInitialState();
	    state = new SpacecraftState(new CartesianOrbit(new PVCoordinates(new Vector3D(X[0], X[1], X[2]), new Vector3D(X[3], X[4], X[5])),
							   cfg.propInertialFrame, cfg.propStart, Constants.EGM96_EARTH_MU), cfg.rsoMass);
	    NumericalPropagator np = new NumericalPropagator(
		new DormandPrince853Integrator(cfg.integMinTimeStep, cfg.integMaxTimeStep, cfg.integAbsTolerance, cfg.integRelTolerance));
	    for (ForceModel fm: cfg.forces)
		np.addForceModel(fm);
	    np.setInitialState(state);
	    prop = np;
	}

	AttitudeProvider attProv = cfg.getAttitudeProvider();
	if (attProv != null)
	    prop.setAttitudeProvider(attProv);

	propagators.add(prop);
	eventHandlers.add(cfg.addEventHandlers(prop, state));
	if (cfg.propStep != 0.0)
	    stepEnd.add(cfg.propStart);
	else
	    throw(new RuntimeException("Invalid propagation step size"));
    }

    private static class Simulation
    {
	private static Random rand = new Random(1);
	private static double[] oneOnes = {1.0};
	private static double[] oneNegInf = {Double.NEGATIVE_INFINITY};
	private static double[] onePosInf = {Double.POSITIVE_INFINITY};
	private static double[] twoZeros = {0.0, 0.0};
	private static double[] twoOnes = {1.0, 1.0};
	private static double[] twoNegInf = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
	private static double[] twoPosInf = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
	private static String[] biasAzEl = {Settings.MeasurementType.AZIMUTH.name(), Settings.MeasurementType.ELEVATION.name()};
	private static String[] biasRaDec = {Settings.MeasurementType.RIGHT_ASCENSION.name(), Settings.MeasurementType.DECLINATION.name()};
	private static String[] biasRange = {Settings.MeasurementType.RANGE.name()};
	private static String[] biasRangeRate = {Settings.MeasurementType.RANGE_RATE.name()};
	private static ObservableSatellite obsSat = new ObservableSatellite(0);
	private static SpacecraftState[] ssStates = new SpacecraftState[1];

	public static boolean simulate(Settings simCfg, SpacecraftState state, ArrayList<EventHandling> handlers,
				       ArrayList<Measurements.Measurement> measOut)
	{
	    ArrayList<Double> extraStates = new ArrayList<Double>();
	    TimeStampedPVCoordinates pvc = state.getPVCoordinates(simCfg.propInertialFrame);
	    if (simCfg.atmModel != null && (simCfg.outputFlags & Settings.OUTPUT_DENSITY) != 0)
		extraStates.add(simCfg.atmModel.getDensity(pvc.getDate(), pvc.getPosition(), simCfg.propInertialFrame));
	    if ((simCfg.outputFlags & Settings.OUTPUT_ECLIPSE) != 0)
	    {
		double eclipse = 0.0;
		for (EventHandling hnd: handlers)
		{
		    if (hnd.stationName != null)
		    {
			if (hnd.stationName.equals(EventHandling.UMBRA) && hnd.detected)
			{
			    eclipse = 1.0;
			    break;
			}
			if (hnd.stationName.equals(EventHandling.PENUMBRA) && hnd.detected)
			    eclipse = 0.5;
		    }
		}
		extraStates.add(eclipse);
	    }

	    Measurements.Measurement meas = new Measurements.Measurement(pvc, extraStates);
	    if (!simCfg.simMeasurements)
	    {
		boolean isVisible = true;
		for (EventHandling hnd: handlers)
		{
		    if (hnd.detected != null && hnd.stationName != null && hnd.stationName.equals(EventHandling.GEO_ZONE_NAME))
		    {
			isVisible = hnd.detected;
			break;
		    }
		}

		if (isVisible)
		    measOut.add(meas);
		return(isVisible);
	    }

	    ssStates[0] = state;
	    boolean continueProp = false;
	    for (Map.Entry<String, GroundStation> kv: simCfg.stations.entrySet())
	    {
		boolean isVisible = false;
		String gsName = kv.getKey();
		for (EventHandling hnd: handlers)
		{
		    if (hnd.detected != null && hnd.stationName != null && hnd.stationName.equals(gsName))
		    {
			isVisible = hnd.detected;
			break;
		    }
		}
		if (!isVisible)
		    continue;

		double[] bias;
		continueProp = true;
		Settings.Station jsn = simCfg.cfgStations.get(gsName);
		if (jsn.bias != null && jsn.bias.length > 0)
		    bias = jsn.bias;
		else
		    bias = twoZeros;

		GroundStation gst = kv.getValue();
		Measurements.Measurement clone = new Measurements.Measurement(meas);
		clone.station = gsName;
		measOut.add(clone);

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
			    clone.values = new double[simCfg.cfgMeasurements.size()];
			clone.values[0] = obs.estimate(0, 0, ssStates).getEstimatedValue()[0];
		    }
		    else if (name == Settings.MeasurementType.RANGE_RATE)
		    {
			RangeRate obs = new RangeRate(gst, pvc.getDate(), 0.0, 0.0, 1.0, val.twoWay, obsSat);
			obs.addModifier(new Bias<RangeRate>(biasRangeRate, new double[]{rand.nextGaussian()*val.error[0]+bias[bias.length-1]},
							    oneOnes, oneNegInf, onePosInf));
			if (clone.values == null)
			    clone.values = new double[simCfg.cfgMeasurements.size()];
			clone.values[clone.values.length-1] = obs.estimate(0, 0, ssStates).getEstimatedValue()[0];
		    }
		    else if (name == Settings.MeasurementType.RIGHT_ASCENSION || name == Settings.MeasurementType.DECLINATION)
		    {
			AngularRaDec obs = new AngularRaDec(gst, simCfg.propInertialFrame, pvc.getDate(), twoZeros, twoZeros, twoOnes, obsSat);
			obs.addModifier(new Bias<AngularRaDec>(biasRaDec, new double[]{rand.nextGaussian()*val.error[0]+bias[0],
										       rand.nextGaussian()*val.error[0]+bias[1]},
				twoOnes, twoNegInf, twoPosInf));
			clone.values = obs.estimate(0, 0, ssStates).getEstimatedValue();
			break;
		    }
		    else if (name == Settings.MeasurementType.AZIMUTH || name == Settings.MeasurementType.ELEVATION)
		    {
			AngularAzEl obs = new AngularAzEl(gst, pvc.getDate(), twoZeros, twoZeros, twoOnes, obsSat);
			obs.addModifier(new AngularRadioRefractionModifier(new EarthITU453AtmosphereRefraction(jsn.altitude)));
			obs.addModifier(new Bias<AngularAzEl>(biasAzEl, new double[]{rand.nextGaussian()*val.error[0]+bias[0],
										     rand.nextGaussian()*val.error[0]+bias[1]},
				twoOnes, twoNegInf, twoPosInf));
			clone.values = obs.estimate(0, 0, ssStates).getEstimatedValue();
			break;
		    }
		    else if (name == Settings.MeasurementType.POSITION || name == Settings.MeasurementType.POSITION_VELOCITY)
		    {
			double[] pos = pvc.getPosition().toArray();
			if (name == Settings.MeasurementType.POSITION)
			    clone.values = pos;
			else
			{
			    double[] vel = pvc.getVelocity().toArray();
			    clone.values = new double[]{pos[0], pos[1], pos[2], vel[0], vel[1], vel[2]};
			}
			for (int i = 0; i < clone.values.length && i < val.error.length; i++)
			    clone.values[i] += rand.nextGaussian()*val.error[i];
		    }
		}
	    }
	    return(continueProp);
	}
    }
}
