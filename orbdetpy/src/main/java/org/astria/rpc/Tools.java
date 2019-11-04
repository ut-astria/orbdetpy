/*
 * Tools.java - RPC utility functions.
 * Copyright (C) 2019 University of Texas
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

package org.astria.rpc;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.astria.Estimation;
import org.astria.Measurements;
import org.astria.ParallelPropagation;
import org.astria.Settings;

public final class Tools
{
    public static Settings buildSettingsFromRequest(Messages.Settings req)
    {
	Settings cfg = new Settings();
	cfg.rsoMass = req.getRsoMass();
	cfg.rsoArea = req.getRsoArea();
	if (req.getRsoFacetsCount() > 0)
	{
	    cfg.rsoFacets = new Settings.Facet[req.getRsoFacetsCount()];
	    for (int i = 0; i < cfg.rsoFacets.length; i++)
	    {
		cfg.rsoFacets[i] = new Settings.Facet(req.getRsoFacets(i).getNormalList().stream().mapToDouble(Double::doubleValue).toArray(),
						      req.getRsoFacets(i).getArea());
	    }
	}
	if (req.getRsoSolarArrayAxisCount() > 0)
	    cfg.rsoSolarArrayAxis = req.getRsoSolarArrayAxisList().stream().mapToDouble(Double::doubleValue).toArray();
	cfg.rsoSolarArrayArea = req.getRsoSolarArrayArea();
	if (req.getRsoAttitudeProvider().length() > 0)
	    cfg.rsoAttitudeProvider = req.getRsoAttitudeProvider();
	if (req.getRsoSpinVelocityCount() > 0)
	    cfg.rsoSpinVelocity = req.getRsoSpinVelocityList().stream().mapToDouble(Double::doubleValue).toArray();
	if (req.getRsoSpinAccelerationCount() > 0)
	    cfg.rsoSpinAcceleration = req.getRsoSpinAccelerationList().stream().mapToDouble(Double::doubleValue).toArray();

	cfg.gravityDegree = req.getGravityDegree();
	cfg.gravityOrder = req.getGravityOrder();
	cfg.oceanTidesDegree = req.getOceanTidesDegree();
	cfg.oceanTidesOrder = req.getOceanTidesOrder();
	cfg.thirdBodySun = req.getThirdBodySun();
	cfg.thirdBodyMoon = req.getThirdBodyMoon();
	cfg.solidTidesSun = req.getSolidTidesSun();
	cfg.solidTidesMoon = req.getSolidTidesMoon();

	cfg.dragModel = req.getDragModel();
	cfg.dragCoefficient = new Settings.Parameter(req.getDragCoefficient().getName(), req.getDragCoefficient().getMin(),
						     req.getDragCoefficient().getMax(), req.getDragCoefficient().getValue(),
						     req.getDragCoefficient().getEstimation());
	cfg.dragMSISEFlags = unpackInteger2DArray(req.getDragMSISEFlagsList());
	cfg.dragExpRho0 = req.getDragExpRho0();
	cfg.dragExpH0 = req.getDragExpH0();
	cfg.dragExpHscale = req.getDragExpHscale();

	cfg.rpSun = req.getRpSun();
	cfg.rpCoeffReflection = new Settings.Parameter(req.getRpCoeffReflection().getName(), req.getRpCoeffReflection().getMin(),
						       req.getRpCoeffReflection().getMax(), req.getRpCoeffReflection().getValue(),
						       req.getRpCoeffReflection().getEstimation());
	cfg.rpCoeffAbsorption = req.getRpCoeffAbsorption();

	if (req.getManeuversCount() > 0)
	{
	    cfg.cfgManeuvers = new Settings.Maneuver[req.getManeuversCount()];
	    for (int i = 0; i < cfg.cfgManeuvers.length; i++)
	    {
		cfg.cfgManeuvers[i] = new Settings.Maneuver(req.getManeuvers(i).getTime(), req.getManeuvers(i).getTriggerEvent(),
							    req.getManeuvers(i).getTriggerParamsList().stream().mapToDouble(Double::doubleValue).toArray(),
							    req.getManeuvers(i).getManeuverType(),
							    req.getManeuvers(i).getManeuverParamsList().stream().mapToDouble(Double::doubleValue).toArray());
	    }
	}

	if (req.getPropStart().length() > 0)
	    cfg.propStart = req.getPropStart();
	if (req.getPropEnd().length() > 0)
	    cfg.propEnd = req.getPropEnd();
	cfg.propStep = req.getPropStep();
	if (req.getPropInitialStateCount() > 0)
	    cfg.propInitialState = req.getPropInitialStateList().stream().mapToDouble(Double::doubleValue).toArray();
	if (req.getPropInitialTLECount() > 0)
	    cfg.propInitialTLE = req.getPropInitialTLEList().toArray(new String[0]);
	if (req.getPropInertialFrame().length() > 0)
	    cfg.propInertialFrame = req.getPropInertialFrame();
	if (req.getPropStepHandlerStartTime().length() > 0)
	    cfg.propStepHandlerStartTime = req.getPropStepHandlerStartTime();
	if (req.getPropStepHandlerEndTime().length() > 0)
	    cfg.propStepHandlerEndTime = req.getPropStepHandlerEndTime();

	cfg.integMinTimeStep = req.getIntegMinTimeStep();
	cfg.integMaxTimeStep = req.getIntegMaxTimeStep();
	cfg.integAbsTolerance = req.getIntegAbsTolerance();
	cfg.integRelTolerance = req.getIntegRelTolerance();

	cfg.simMeasurements = req.getSimMeasurements();
	cfg.simSkipUnobservable = req.getSimSkipUnobservable();
	cfg.simIncludeExtras = req.getSimIncludeExtras();
	cfg.simIncludeStationState = req.getSimIncludeStationState();

	if (req.getStationsCount() > 0)
	{
	    cfg.cfgStations = new HashMap<String, Settings.Station>(req.getStationsCount());
	    for (Map.Entry<String, Messages.Station> kv : req.getStationsMap().entrySet())
	    {
		Messages.Station v = kv.getValue();
		cfg.cfgStations.put(kv.getKey(), new Settings.Station(v.getLatitude(), v.getLongitude(), v.getAltitude(), v.getAzimuthBias(),
								      v.getElevationBias(), v.getRangeBias(), v.getRangeRateBias(),
								      v.getRightAscensionBias(), v.getDeclinationBias(),
								      v.getPositionBiasList().stream().mapToDouble(Double::doubleValue).toArray(),
								      v.getPositionVelocityBiasList().stream().mapToDouble(Double::doubleValue).toArray()));
	    }
	}

	if (req.getMeasurementsCount() > 0)
	{
	    cfg.cfgMeasurements = new HashMap<String, Settings.Measurement>(req.getMeasurementsCount());
	    for (Map.Entry<String, Messages.MeasurementSetting> kv : req.getMeasurementsMap().entrySet())
	    {
		Messages.MeasurementSetting v = kv.getValue();
		cfg.cfgMeasurements.put(kv.getKey(), new Settings.Measurement(v.getTwoWay(),
									      v.getErrorList().stream().mapToDouble(Double::doubleValue).toArray()));
	    }
	}

	if (req.getEstmFilter().length() > 0)
	    cfg.estmFilter = req.getEstmFilter();
	if (req.getEstmCovarianceCount() > 0)
	    cfg.estmCovariance = req.getEstmCovarianceList().stream().mapToDouble(Double::doubleValue).toArray();
	if (req.getEstmProcessNoiseCount() > 0)
	    cfg.estmProcessNoise = req.getEstmProcessNoiseList().stream().mapToDouble(Double::doubleValue).toArray();
	cfg.estmDMCCorrTime = req.getEstmDMCCorrTime();
	cfg.estmDMCSigmaPert = req.getEstmDMCSigmaPert();
	cfg.estmDMCAcceleration = new Settings.Parameter("", req.getEstmDMCAcceleration().getMin(), req.getEstmDMCAcceleration().getMax(),
							 req.getEstmDMCAcceleration().getValue(), req.getEstmDMCAcceleration().getEstimation());
	return(cfg.build());
    }

    public static Measurements buildMeasurementsFromRequest(List<Messages.Measurement> req, Settings config)
    {
	Measurements output = new Measurements();
	output.rawMeas = new Measurements.Measurement[req.size()];
	for (int i = 0; i < output.rawMeas.length; i++)
	{
	    Messages.Measurement min = req.get(i);
	    output.rawMeas[i] = new Measurements.Measurement();
	    if (min.getTime().length() > 0)
		output.rawMeas[i].Time = min.getTime();
	    if (min.getStation().length() > 0)
		output.rawMeas[i].Station = min.getStation();
	    if (min.getAzimuth() != 0.0)
		output.rawMeas[i].Azimuth = min.getAzimuth();
	    if (min.getElevation() != 0.0)
		output.rawMeas[i].Elevation = min.getElevation();
	    if (min.getRange() != 0.0)
		output.rawMeas[i].Range = min.getRange();
	    if (min.getRangeRate() != 0.0)
		output.rawMeas[i].RangeRate = min.getRangeRate();
	    if (min.getRightAscension() != 0.0)
		output.rawMeas[i].RightAscension = min.getRightAscension();
	    if (min.getDeclination() != 0.0)
		output.rawMeas[i].Declination = min.getDeclination();
	    if (min.getPositionCount() > 0)
		output.rawMeas[i].Position = min.getPositionList().toArray(new Double[0]);
	    if (min.getPositionVelocityCount() > 0)
		output.rawMeas[i].PositionVelocity = min.getPositionVelocityList().toArray(new Double[0]);
	}

	return(output.build(config));
    }

    public static List<Messages.PropagationOutput> buildResponseFromPropagation(ArrayList<ParallelPropagation.PropagationOutput> plist)
    {
	ArrayList<Messages.PropagationOutput> output = new ArrayList<Messages.PropagationOutput>(plist.size());
	for (ParallelPropagation.PropagationOutput pin: plist)
	{
	    Messages.PropagationOutput.Builder builder = Messages.PropagationOutput.newBuilder().setTime(pin.time);
	    for (double[] state: pin.states)
		builder = builder.addStates(Messages.DoubleArray.newBuilder()
					    .addAllArray(DoubleStream.of(state).boxed().collect(Collectors.toList())).build());
	    output.add(builder.build());
	}

	return(output);
    }

    public static List<Messages.Measurement> buildResponseFromMeasurements(ArrayList<Measurements.SimulatedMeasurement> mlist)
    {
	ArrayList<Messages.Measurement> output = new ArrayList<Messages.Measurement>(mlist.size());
	for (Measurements.SimulatedMeasurement min: mlist)
	{
	    Messages.Measurement.Builder builder = Messages.Measurement.newBuilder().setTime(min.Time);
	    if (min.Station != null && min.Station.length() > 0)
		builder = builder.setStation(min.Station);
	    if (min.Azimuth != null)
		builder = builder.setAzimuth(min.Azimuth);
	    if (min.Elevation != null)
		builder = builder.setElevation(min.Elevation);
	    if (min.Range != null)
		builder = builder.setRange(min.Range);
	    if (min.RangeRate != null)
		builder = builder.setRangeRate(min.RangeRate);
	    if (min.RightAscension != null)
		builder = builder.setRightAscension(min.RightAscension);
	    if (min.Declination != null)
		builder = builder.setDeclination(min.Declination);
	    if (min.Position != null)
		builder = builder.addAllPosition(Arrays.asList(min.Position));
	    if (min.PositionVelocity != null)
		builder = builder.addAllPositionVelocity(Arrays.asList(min.PositionVelocity));
	    if (min.TrueState != null)
	    {
		builder = builder.addAllTrueStateCartesian(DoubleStream.of(min.TrueState.Cartesian).boxed().collect(Collectors.toList()));
		builder = builder.setTrueStateSma(min.TrueState.Kepler.SMA);
		builder = builder.setTrueStateEcc(min.TrueState.Kepler.Ecc);
		builder = builder.setTrueStateInc(min.TrueState.Kepler.Inc);
		builder = builder.setTrueStateRaan(min.TrueState.Kepler.RAAN);
		builder = builder.setTrueStateArgp(min.TrueState.Kepler.ArgP);
		builder = builder.setTrueStateMeanAnom(min.TrueState.Kepler.MeanAnom);
		builder = builder.setTrueStateEx(min.TrueState.Equinoctial.Ex);
		builder = builder.setTrueStateEy(min.TrueState.Equinoctial.Ey);
		builder = builder.setTrueStateHx(min.TrueState.Equinoctial.Hx);
		builder = builder.setTrueStateHy(min.TrueState.Equinoctial.Hy);
		builder = builder.setTrueStateLm(min.TrueState.Equinoctial.Lm);
	    }
	    if (min.AtmDensity != null)
		builder = builder.setAtmosphericDensity(min.AtmDensity);
	    if (min.AccGravity != null)
		builder = builder.addAllAccelerationGravity(
		    DoubleStream.of(min.AccGravity).boxed().collect(Collectors.toList()));
	    if (min.AccDrag != null)
		builder = builder.addAllAccelerationDrag(
		    DoubleStream.of(min.AccDrag).boxed().collect(Collectors.toList()));
	    if (min.AccOceanTides != null)
		builder = builder.addAllAccelerationOceanTides(
		    DoubleStream.of(min.AccOceanTides).boxed().collect(Collectors.toList()));
	    if (min.AccSolidTides != null)
		builder = builder.addAllAccelerationSolidTides(
		    DoubleStream.of(min.AccSolidTides).boxed().collect(Collectors.toList()));
	    if (min.AccThirdBodies != null)
		builder = builder.addAllAccelerationThirdBodies(
		    DoubleStream.of(min.AccThirdBodies).boxed().collect(Collectors.toList()));
	    if (min.AccRadiationPressure != null)
		builder = builder.addAllAccelerationRadiationPressure(
		    DoubleStream.of(min.AccRadiationPressure).boxed().collect(Collectors.toList()));
	    if (min.AccThrust != null)
		builder = builder.addAllAccelerationThrust(
		    DoubleStream.of(min.AccThrust).boxed().collect(Collectors.toList()));
	    if (min.StationState != null)
		builder = builder.addAllStationState(
		    DoubleStream.of(min.StationState).boxed().collect(Collectors.toList()));
	    output.add(builder.build());
	}

	return(output);
    }

    public static List<Messages.EstimationOutput> buildResponseFromOrbitDetermination(ArrayList<Estimation.EstimationOutput> elist)
    {
	List<Messages.DoubleArray> dub;
	ArrayList<Messages.EstimationOutput> output = new ArrayList<Messages.EstimationOutput>(elist.size());
	for (Estimation.EstimationOutput ein: elist)
	{
	    Messages.EstimationOutput.Builder builder = Messages.EstimationOutput.newBuilder()
		.setTime(ein.Time)
		.addAllEstimatedState(DoubleStream.of(ein.EstimatedState).boxed().collect(Collectors.toList()));
	    if (ein.Station != null && ein.Station.length() > 0)
		builder = builder.setStation(ein.Station);
	    dub = packDouble2DArray(ein.PropagatedCovariance);
	    if (dub != null)
		builder = builder.addAllPropagatedCovariance(dub);
	    dub = packDouble2DArray(ein.InnovationCovariance);
	    if (dub != null)
		builder = builder.addAllInnovationCovariance(dub);
	    dub = packDouble2DArray(ein.EstimatedCovariance);
	    if (dub != null)
		builder = builder.addAllEstimatedCovariance(dub);
	    if (ein.PreFit != null)
	    {
		for (Map.Entry<String, double[]> kv : ein.PreFit.entrySet())
		{
		    builder = builder.putPreFit(kv.getKey(), Messages.DoubleArray.newBuilder().addAllArray(
						    DoubleStream.of(kv.getValue()).boxed().collect(Collectors.toList())).build());
		}
	    }
	    if (ein.PostFit != null)
	    {
		for (Map.Entry<String, double[]> kv : ein.PostFit.entrySet())
		{
		    builder = builder.putPostFit(kv.getKey(), Messages.DoubleArray.newBuilder().addAllArray(
						     DoubleStream.of(kv.getValue()).boxed().collect(Collectors.toList())).build());
		}
	    }
	    output.add(builder.build());
	}

	return(output);
    }

    public static List<Messages.DoubleArray> packDouble2DArray(double[][] in)
    {
	if (in == null || in.length == 0)
	    return(null);
	List<Messages.DoubleArray> out = new ArrayList<Messages.DoubleArray>(in.length);
	for (int i = 0; i < in.length; i++)
	{
	    Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder()
		.addAllArray(DoubleStream.of(in[i]).boxed().collect(Collectors.toList()));
	    out.add(builder.build());
	}

	return(out);
    }
    
    public static int[][] unpackInteger2DArray(List<Messages.IntegerArray> in)
    {
	if (in.size() == 0)
	    return(null);
	int[][] out = new int[in.size()][];
	for (int i = 0; i < out.length; i++)
	    out[i] = in.get(i).getArrayList().stream().mapToInt(Integer::intValue).toArray();
	return(out);
    }

    public static String getStackTrace(Throwable exc)
    {
	StringWriter sw = new StringWriter();
	exc.printStackTrace(new PrintWriter(sw));
	return(sw.toString());
    }
}
