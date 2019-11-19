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
								      v.getPositionVelocityBiasList().stream().mapToDouble(Double::doubleValue).toArray(),
								      v.getBiasEstimation()));
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
		output.rawMeas[i].time = min.getTime();
	    if (min.getStation().length() > 0)
		output.rawMeas[i].station = min.getStation();
	    if (min.getAzimuth() != 0.0)
		output.rawMeas[i].azimuth = min.getAzimuth();
	    if (min.getElevation() != 0.0)
		output.rawMeas[i].elevation = min.getElevation();
	    if (min.getRange() != 0.0)
		output.rawMeas[i].range = min.getRange();
	    if (min.getRangeRate() != 0.0)
		output.rawMeas[i].rangeRate = min.getRangeRate();
	    if (min.getRightAscension() != 0.0)
		output.rawMeas[i].rightAscension = min.getRightAscension();
	    if (min.getDeclination() != 0.0)
		output.rawMeas[i].declination = min.getDeclination();
	    if (min.getPositionCount() > 0)
		output.rawMeas[i].position = min.getPositionList().stream().mapToDouble(Double::doubleValue).toArray();
	    if (min.getPositionVelocityCount() > 0)
		output.rawMeas[i].positionVelocity = min.getPositionVelocityList().stream().mapToDouble(Double::doubleValue).toArray();
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
	    Messages.Measurement.Builder builder = Messages.Measurement.newBuilder().setTime(min.time);
	    if (min.station != null && min.station.length() > 0)
		builder = builder.setStation(min.station);
	    if (min.azimuth != 0.0)
		builder = builder.setAzimuth(min.azimuth);
	    if (min.elevation != 0.0)
		builder = builder.setElevation(min.elevation);
	    if (min.range != 0.0)
		builder = builder.setRange(min.range);
	    if (min.rangeRate != 0.0)
		builder = builder.setRangeRate(min.rangeRate);
	    if (min.rightAscension != 0.0)
		builder = builder.setRightAscension(min.rightAscension);
	    if (min.declination != 0.0)
		builder = builder.setDeclination(min.declination);
	    if (min.position != null)
		builder = builder.addAllPosition(DoubleStream.of(min.position).boxed().collect(Collectors.toList()));
	    if (min.positionVelocity != null)
		builder = builder.addAllPositionVelocity(DoubleStream.of(min.positionVelocity).boxed().collect(Collectors.toList()));
	    if (min.trueState != null)
	    {
		builder = builder.addAllTrueStateCartesian(DoubleStream.of(min.trueState.cartesian).boxed().collect(Collectors.toList()));
		builder = builder.setTrueStateSma(min.trueState.keplerian.sma);
		builder = builder.setTrueStateEcc(min.trueState.keplerian.ecc);
		builder = builder.setTrueStateInc(min.trueState.keplerian.inc);
		builder = builder.setTrueStateRaan(min.trueState.keplerian.raan);
		builder = builder.setTrueStateArgp(min.trueState.keplerian.argP);
		builder = builder.setTrueStateMeanAnom(min.trueState.keplerian.meanAnom);
		builder = builder.setTrueStateEx(min.trueState.equinoctial.ex);
		builder = builder.setTrueStateEy(min.trueState.equinoctial.ey);
		builder = builder.setTrueStateHx(min.trueState.equinoctial.hx);
		builder = builder.setTrueStateHy(min.trueState.equinoctial.hy);
		builder = builder.setTrueStateLm(min.trueState.equinoctial.lm);
	    }
	    if (min.atmDensity != 0.0)
		builder = builder.setAtmosphericDensity(min.atmDensity);
	    if (min.accGravity != null)
		builder = builder.addAllAccelerationGravity(
		    DoubleStream.of(min.accGravity).boxed().collect(Collectors.toList()));
	    if (min.accDrag != null)
		builder = builder.addAllAccelerationDrag(
		    DoubleStream.of(min.accDrag).boxed().collect(Collectors.toList()));
	    if (min.accOceanTides != null)
		builder = builder.addAllAccelerationOceanTides(
		    DoubleStream.of(min.accOceanTides).boxed().collect(Collectors.toList()));
	    if (min.accSolidTides != null)
		builder = builder.addAllAccelerationSolidTides(
		    DoubleStream.of(min.accSolidTides).boxed().collect(Collectors.toList()));
	    if (min.accThirdBodies != null)
		builder = builder.addAllAccelerationThirdBodies(
		    DoubleStream.of(min.accThirdBodies).boxed().collect(Collectors.toList()));
	    if (min.accRadiationPressure != null)
		builder = builder.addAllAccelerationRadiationPressure(
		    DoubleStream.of(min.accRadiationPressure).boxed().collect(Collectors.toList()));
	    if (min.accThrust != null)
		builder = builder.addAllAccelerationThrust(
		    DoubleStream.of(min.accThrust).boxed().collect(Collectors.toList()));
	    if (min.stationState != null)
		builder = builder.addAllStationState(
		    DoubleStream.of(min.stationState).boxed().collect(Collectors.toList()));
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
		.setTime(ein.time)
		.addAllEstimatedState(DoubleStream.of(ein.estimatedState).boxed().collect(Collectors.toList()));
	    if (ein.station != null && ein.station.length() > 0)
		builder = builder.setStation(ein.station);
	    dub = packDouble2DArray(ein.propagatedCovariance);
	    if (dub != null)
		builder = builder.addAllPropagatedCovariance(dub);
	    dub = packDouble2DArray(ein.innovationCovariance);
	    if (dub != null)
		builder = builder.addAllInnovationCovariance(dub);
	    dub = packDouble2DArray(ein.estimatedCovariance);
	    if (dub != null)
		builder = builder.addAllEstimatedCovariance(dub);
	    if (ein.preFit != null)
	    {
		for (Map.Entry<String, double[]> kv : ein.preFit.entrySet())
		{
		    builder = builder.putPreFit(kv.getKey(), Messages.DoubleArray.newBuilder().addAllArray(
						    DoubleStream.of(kv.getValue()).boxed().collect(Collectors.toList())).build());
		}
	    }
	    if (ein.postFit != null)
	    {
		for (Map.Entry<String, double[]> kv : ein.postFit.entrySet())
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
