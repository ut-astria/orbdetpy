/*
 * Simulation.java - Functions to simulate noisy measurements.
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
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.measurements.RangeRate;
import org.orekit.estimation.measurements.modifiers.Bias;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.numerical.NumericalPropagator;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class Simulation
{
    private final Settings simCfg;

    public Simulation(Settings simCfg)
    {
	this.simCfg = simCfg;
    }

    public ArrayList<Measurements.SimulatedMeasurement> simulateMeasurements()
    {
	final NumericalPropagator propagator = new NumericalPropagator(
	    new DormandPrince853Integrator(simCfg.integMinTimeStep, simCfg.integMaxTimeStep,
					   simCfg.integAbsTolerance, simCfg.integRelTolerance));
	for (ForceModel fm : simCfg.forces)
	    propagator.addForceModel(fm);

	final AttitudeProvider attpro = simCfg.getAttitudeProvider();
	if (attpro != null)
	    propagator.setAttitudeProvider(attpro);

	final double[] Xi = simCfg.getInitialState();
	AbsoluteDate tm = DataManager.parseDateTime(simCfg.propStart);
	final AbsoluteDate prend = DataManager.parseDateTime(simCfg.propEnd);

	final SpacecraftState sstate0 = new SpacecraftState(
	    new CartesianOrbit(new PVCoordinates(new Vector3D(Xi[0], Xi[1], Xi[2]), new Vector3D(Xi[3], Xi[4], Xi[5])),
			       simCfg.propFrame, tm, Constants.EGM96_EARTH_MU), simCfg.rsoMass);
	propagator.setInitialState(sstate0);
	simCfg.addEventHandlers(propagator, sstate0);
	final boolean evthandlers = propagator.getEventsDetectors().size() > 0;

	final Random rand = new Random(1);
	final double[] oneOnes = new double[] {1.0};
	final double[] oneNegInf = new double[] {Double.NEGATIVE_INFINITY};
	final double[] onePosInf = new double[] {Double.POSITIVE_INFINITY};
	final double[] twoZeros = new double[]{0.0, 0.0};
	final double[] twoOnes = new double[] {1.0, 1.0};
	final double[] twoNegInf = new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
	final double[] twoPosInf = new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
	final String[] biasAzEl = new String[] {"Az", "El"};
	final String[] biasRaDec = new String[] {"RA", "Dec"};
	final String[] biasRange = new String[] {"range"};
	final String[] biasRangeRate = new String[] {"rangeRate"};
	final ObservableSatellite obsSat = new ObservableSatellite(0);
	final ArrayList<Measurements.SimulatedMeasurement> mall = new ArrayList<Measurements.SimulatedMeasurement>(
	    (int) FastMath.abs(prend.durationFrom(tm)/simCfg.propStep) + 2);

	while (true)
	{
	    final SpacecraftState[] sta = new SpacecraftState[]{propagator.propagate(tm)};
	    final Orbit orb = sta[0].getOrbit();
	    final KeplerianOrbit keporb = new KeplerianOrbit(orb);
	    final AbsoluteDate proptm = sta[0].getDate();
	    TimeStampedPVCoordinates pvc = sta[0].getPVCoordinates();
	    Vector3D pos = pvc.getPosition();
	    Vector3D vel = pvc.getVelocity();
	    Vector3D acc = pvc.getAcceleration();

	    Measurements.SimulatedMeasurement json = new Measurements.SimulatedMeasurement();
	    json.time = DataManager.getUTCString(proptm);
	    json.trueState = new Measurements.State();
	    json.trueState.cartesian = new double[]{pos.getX(), pos.getY(), pos.getZ(), vel.getX(), vel.getY(), vel.getZ(),
						    acc.getX(), acc.getY(), acc.getZ()};

	    if (simCfg.simIncludeExtras)
	    {
		json.trueState.keplerian = new Measurements.KeplerianElements(keporb.getA(), keporb.getE(), keporb.getI(),
									      keporb.getRightAscensionOfAscendingNode(),
									      keporb.getPerigeeArgument(), keporb.getMeanAnomaly());
		json.trueState.equinoctial = new Measurements.EquinoctialElements(orb.getA(), orb.getEquinoctialEx(),
										  orb.getEquinoctialEy(), orb.getHx(),
										  orb.getHy(), orb.getLM());

		if (simCfg.atmModel != null)
		    json.atmDensity = simCfg.atmModel.getDensity(proptm, pos, simCfg.propFrame);

		for (ForceModel fmod : simCfg.forces)
		{
		    double[] facc = fmod.acceleration(sta[0], fmod.getParameters()).toArray();
		    String ftype = fmod.getClass().getSimpleName();
		    if (ftype.equals("HolmesFeatherstoneAttractionModel") || ftype.equals("NewtonianAttraction"))
			json.accGravity = facc;
		    if (ftype.equals("DragForce"))
			json.accDrag = facc;
		    if (ftype.equals("SolidTides"))
			json.accSolidTides = facc;
		    if (ftype.equals("OceanTides"))
			json.accOceanTides = facc;
		    if (ftype.equals("SolarRadiationPressure"))
			json.accRadiationPressure = facc;
		    if (ftype.equals("ConstantThrustManeuver"))
			json.accThrust = facc;
		    if (ftype.equals("ThirdBodyAttraction"))
		    {
			if (json.accThirdBodies == null)
			    json.accThirdBodies = facc;
			else
			{
			    for (int ii = 0; ii < 3; ii++)
				json.accThirdBodies[ii] += facc[ii];
			}
		    }
		}
	    }

	    if (simCfg.simMeasurements)
	    {
		boolean added = false;
		for (Map.Entry<String, GroundStation> kv : simCfg.stations.entrySet())
		{
		    GroundStation gst = kv.getValue();
		    Settings.Station jsn = simCfg.cfgStations.get(kv.getKey());
		    AngularAzEl azel = new AngularAzEl(gst, proptm, twoZeros, twoZeros, twoOnes, obsSat);
		    if (simCfg.simSkipUnobservable && azel.estimate(0, 0, sta).getEstimatedValue()[1] <= 5E-6)
			continue;

		    Measurements.SimulatedMeasurement clone = new Measurements.SimulatedMeasurement(json);
		    clone.station = kv.getKey();
		    if (simCfg.simIncludeStationState)
		    {
			pvc = gst.getBaseFrame().getPVCoordinates(proptm, simCfg.propFrame);
			pos = pvc.getPosition();
			vel = pvc.getVelocity();
			acc = pvc.getAcceleration();
			clone.stationState = new double[]{pos.getX(), pos.getY(), pos.getZ(), vel.getX(), vel.getY(),
							  vel.getZ(), acc.getX(), acc.getY(), acc.getZ()};
		    }

		    for (Map.Entry<String, Settings.Measurement> nvp : simCfg.cfgMeasurements.entrySet())
		    {
			String name = nvp.getKey();
			Settings.Measurement val = nvp.getValue();
			if (name.equalsIgnoreCase("range"))
			{
			    Range obs = new Range(gst, val.twoWay, proptm, 0.0, 0.0, 1.0, obsSat);
			    obs.addModifier(new Bias<Range>(biasRange, new double[] {rand.nextGaussian()*val.error[0] + jsn.rangeBias},
							    oneOnes, oneNegInf, onePosInf));
			    clone.range = obs.estimate(0, 0, sta).getEstimatedValue()[0];
			}
			else if (name.equalsIgnoreCase("rangeRate"))
			{
			    RangeRate obs = new RangeRate(gst, proptm, 0.0, 0.0, 1.0, val.twoWay, obsSat);
			    obs.addModifier(new Bias<RangeRate>(biasRangeRate, new double[] {rand.nextGaussian()*val.error[0] + jsn.rangeRateBias},
								oneOnes, oneNegInf, onePosInf));
			    clone.rangeRate = obs.estimate(0, 0, sta).getEstimatedValue()[0];
			}
			else if (name.equalsIgnoreCase("rightAscension") || name.equalsIgnoreCase("declination") && clone.rightAscension == 0.0)
			{
			    AngularRaDec obs = new AngularRaDec(gst, simCfg.propFrame, proptm, twoZeros, twoZeros, twoOnes,
								obsSat);
			    obs.addModifier(new Bias<AngularRaDec>(biasRaDec, new double[] {rand.nextGaussian()*val.error[0] + jsn.rightAscensionBias,
											    rand.nextGaussian()*val.error[0] + jsn.declinationBias},
				    twoOnes, twoNegInf, twoPosInf));
			    double[] obsVal = obs.estimate(0, 0, sta).getEstimatedValue();
			    clone.rightAscension = obsVal[0];
			    clone.declination = obsVal[1];
			}
			else if (name.equalsIgnoreCase("azimuth") || name.equalsIgnoreCase("elevation") && clone.azimuth == 0.0)
			{
			    azel.addModifier(new Bias<AngularAzEl>(biasAzEl, new double[] {rand.nextGaussian()*val.error[0] + jsn.azimuthBias,
											   rand.nextGaussian()*val.error[0] + jsn.elevationBias},
				    twoOnes, twoNegInf, twoPosInf));
			    double[] obsVal = azel.estimate(0, 0, sta).getEstimatedValue();
			    clone.azimuth = obsVal[0];
			    clone.elevation = obsVal[1];
			}
			else if (name.equalsIgnoreCase("position"))
			{
			    clone.position = new double[3];
			    for (int i = 0; i < 3; i++)
				clone.position[i] = clone.trueState.cartesian[i] + rand.nextGaussian()*val.error[i] + jsn.positionBias[i];
			}
			else if (name.equalsIgnoreCase("positionVelocity"))
			{
			    clone.positionVelocity = new double[6];
			    for (int i = 0; i < 6; i++)
				clone.positionVelocity[i] = clone.trueState.cartesian[i] + rand.nextGaussian()*val.error[i] + jsn.positionVelocityBias[i];
			}
		    }

		    added = true;
		    mall.add(clone);
		}

		if (!added)
		    mall.add(json);
	    }
	    else
		mall.add(json);

	    double dt = prend.durationFrom(proptm);
	    if (evthandlers && !proptm.equals(tm))
		break;
	    if (simCfg.propStep >= 0.0)
		tm = new AbsoluteDate(tm, FastMath.min(dt, simCfg.propStep));
	    else
	    {
		tm = new AbsoluteDate(tm, FastMath.max(dt, simCfg.propStep));
		dt = -dt;
	    }
	    if (dt <= 0.0)
		break;
	}

	return(mall);
    }
}
