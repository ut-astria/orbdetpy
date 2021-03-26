/*
 * Tools.java - RPC utility functions.
 * Copyright (C) 2019-2021 University of Texas
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
import java.util.stream.Collectors;
import java.util.stream.DoubleStream;
import java.io.PrintWriter;
import java.io.StringWriter;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import org.astria.Estimation;
import org.astria.Measurements;
import org.astria.Settings;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.time.AbsoluteDate;

public final class Tools
{
    public static Settings buildSettingsFromRequest(Messages.Settings req)
    {
	Settings cfg = new Settings();
	Settings.EstimationType[] estmTypes = Settings.EstimationType.values();

	cfg.rsoMass = req.getRsoMass();
	cfg.rsoArea = req.getRsoArea();
	if (req.getRsoFacetsCount() > 0)
	{
	    cfg.rsoFacets = new Settings.Facet[req.getRsoFacetsCount()];
	    for (int i = 0; i < cfg.rsoFacets.length; i++)
		cfg.rsoFacets[i] = new Settings.Facet(req.getRsoFacets(i).getNormalList().stream().mapToDouble(Double::doubleValue).toArray(),
						      req.getRsoFacets(i).getArea());
	}
	if (req.getRsoSolarArrayAxisCount() > 0)
	    cfg.rsoSolarArrayAxis = req.getRsoSolarArrayAxisList().stream().mapToDouble(Double::doubleValue).toArray();
	cfg.rsoSolarArrayArea = req.getRsoSolarArrayArea();
	cfg.rsoAttitudeProvider = Settings.AttitudeType.values()[req.getRsoAttitudeProvider()];
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

	cfg.dragModel = Settings.DragModel.values()[req.getDragModel()];
	cfg.dragCoefficient = new Settings.Parameter("Cd", req.getDragCoefficient().getMin(), req.getDragCoefficient().getMax(),
						     req.getDragCoefficient().getValue(), estmTypes[req.getDragCoefficient().getEstimation()]);
	cfg.dragMSISEFlags = unpackInteger2DArray(req.getDragMSISEFlagsList());
	cfg.dragExpRho0 = req.getDragExpRho0();
	cfg.dragExpH0 = req.getDragExpH0();
	cfg.dragExpHscale = req.getDragExpHscale();

	cfg.rpSun = req.getRpSun();
	cfg.rpCoeffReflection = new Settings.Parameter("Cr", req.getRpCoeffReflection().getMin(), req.getRpCoeffReflection().getMax(),
						       req.getRpCoeffReflection().getValue(), estmTypes[req.getRpCoeffReflection().getEstimation()]);
	cfg.rpCoeffAbsorption = req.getRpCoeffAbsorption();

	if (req.getManeuversCount() > 0)
	{
	    Settings.ManeuverType[] mnvrVals = Settings.ManeuverType.values();
	    Settings.ManeuverTrigger[] trigVals = Settings.ManeuverTrigger.values();
	    cfg.cfgManeuvers = new Settings.Maneuver[req.getManeuversCount()];
	    for (int i = 0; i < cfg.cfgManeuvers.length; i++)
		cfg.cfgManeuvers[i] = new Settings.Maneuver(AbsoluteDate.J2000_EPOCH.shiftedBy(req.getManeuvers(i).getTime()),
							    trigVals[req.getManeuvers(i).getTriggerEvent()],
							    req.getManeuvers(i).getTriggerParamsList().stream().mapToDouble(Double::doubleValue).toArray(),
							    mnvrVals[req.getManeuvers(i).getManeuverType()],
							    req.getManeuvers(i).getManeuverParamsList().stream().mapToDouble(Double::doubleValue).toArray());
	}

	cfg.propStep = req.getPropStep();
	if (req.getPropStart() != 0.0 || req.getPropEnd() != 0.0)
	{
	    cfg.propStart = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getPropStart());
	    cfg.propEnd = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getPropEnd());
	}
	if (req.getPropInitialStateCount() > 0)
	    cfg.propInitialState = req.getPropInitialStateList().stream().mapToDouble(Double::doubleValue).toArray();
	if (req.getPropInitialTLECount() > 0)
	    cfg.propInitialTLE = req.getPropInitialTLEList().toArray(new String[0]);
	if (req.getPropInertialFrame().length() > 0)
	    cfg.propInertialFrame = FramesFactory.getFrame(Predefined.valueOf(req.getPropInertialFrame()));

	cfg.integMinTimeStep = req.getIntegMinTimeStep();
	cfg.integMaxTimeStep = req.getIntegMaxTimeStep();
	cfg.integAbsTolerance = req.getIntegAbsTolerance();
	cfg.integRelTolerance = req.getIntegRelTolerance();
	cfg.simMeasurements = req.getSimMeasurements();

	if (req.getStationsCount() > 0)
	{
	    cfg.cfgStations = new HashMap<String, Settings.Station>(req.getStationsCount());
	    for (Map.Entry<String, Messages.Station> kv : req.getStationsMap().entrySet())
	    {
		Messages.Station v = kv.getValue();
		cfg.cfgStations.put(kv.getKey(),
				    new Settings.Station(v.getLatitude(), v.getLongitude(), v.getAltitude(),
							 v.getBiasList().stream().mapToDouble(Double::doubleValue).toArray(),
							 estmTypes[v.getBiasEstimation()], v.getFovAzimuth(), v.getFovElevation(), v.getFovAperture()));
	    }
	}

	if (req.getMeasurementsCount() > 0)
	{
	    Settings.MeasurementType[] values = Settings.MeasurementType.values();
	    cfg.cfgMeasurements = new HashMap<Settings.MeasurementType, Settings.Measurement>(req.getMeasurementsCount());
	    for (Map.Entry<Integer, Messages.MeasurementSetting> kv: req.getMeasurementsMap().entrySet())
	    {
		Messages.MeasurementSetting v = kv.getValue();
		cfg.cfgMeasurements.put(values[kv.getKey()],
					new Settings.Measurement(v.getTwoWay(), v.getErrorList().stream().mapToDouble(Double::doubleValue).toArray()));
	    }
	}

	if (req.getGeoZoneLatLonCount() > 0)
	    cfg.geoZoneLatLon = req.getGeoZoneLatLonList().stream().mapToDouble(Double::doubleValue).toArray();

	cfg.estmFilter = Settings.Filter.values()[req.getEstmFilter()];
	if (req.getEstmCovarianceCount() > 0)
	    cfg.estmCovariance = req.getEstmCovarianceList().stream().mapToDouble(Double::doubleValue).toArray();
	if (req.getEstmProcessNoiseCount() > 0)
	    cfg.estmProcessNoise = req.getEstmProcessNoiseList().stream().mapToDouble(Double::doubleValue).toArray();
	cfg.estmDMCCorrTime = req.getEstmDMCCorrTime();
	cfg.estmDMCSigmaPert = req.getEstmDMCSigmaPert();
	cfg.estmDMCAcceleration = new Settings.Parameter("DMC", req.getEstmDMCAcceleration().getMin(), req.getEstmDMCAcceleration().getMax(),
							 req.getEstmDMCAcceleration().getValue(), estmTypes[req.getEstmDMCAcceleration().getEstimation()]);
	cfg.estmOutlierSigma = req.getEstmOutlierSigma();
	cfg.estmOutlierWarmup = req.getEstmOutlierWarmup();
	cfg.estmSmootherIterations = req.getEstmSmootherIterations();
	cfg.estmEnablePDAF = req.getEstmEnablePDAF();
	cfg.estmDetectionProbability = req.getEstmDetectionProbability();
	cfg.estmGatingProbability = req.getEstmGatingProbability();
	cfg.estmGatingThreshold = req.getEstmGatingThreshold();
	cfg.outputFlags = req.getOutputFlags();
	return(cfg.build());
    }

    public static Measurements buildMeasurementsFromRequest(List<Messages.Measurement> req, Settings config)
    {
	Measurements output = new Measurements();
	output.array = new Measurements.Measurement[req.size()];
	for (int i = 0; i < output.array.length; i++)
	{
	    Messages.Measurement min = req.get(i);
	    output.array[i] = new Measurements.Measurement();
	    output.array[i].time = AbsoluteDate.J2000_EPOCH.shiftedBy(min.getTime());
	    output.array[i].station = min.getStation();
	    output.array[i].values = min.getValuesList().stream().mapToDouble(Double::doubleValue).toArray();
	    output.array[i].angleRates = min.getAngleRatesList().stream().mapToDouble(Double::doubleValue).toArray();
	}
	return(output.build(config));
    }

    public static List<Messages.Measurement> buildResponseFromMeasurements(ArrayList<Measurements.Measurement> mlist)
    {
	ArrayList<Messages.Measurement> output = new ArrayList<Messages.Measurement>(mlist.size());
	for (Measurements.Measurement min: mlist)
	{
	    Messages.Measurement.Builder builder = Messages.Measurement.newBuilder()
		.setTime(min.time.durationFrom(AbsoluteDate.J2000_EPOCH));
	    if (min.station != null && min.station.length() > 0)
		builder = builder.setStation(min.station);
	    if (min.values != null)
		builder = builder.addAllValues(DoubleStream.of(min.values).boxed().collect(Collectors.toList()));
	    if (min.angleRates != null)
		builder = builder.addAllAngleRates(DoubleStream.of(min.angleRates).boxed().collect(Collectors.toList()));
	    if (min.trueState != null)
		builder = builder.addAllTrueState(DoubleStream.of(min.trueState).boxed().collect(Collectors.toList()));
	    output.add(builder.build());
	}

	return(output);
    }

    public static List<Messages.EstimationOutput> buildResponseFromOrbitDetermination(ArrayList<Estimation.EstimationOutput> elist)
    {
	ArrayList<Messages.EstimationOutput> output = new ArrayList<Messages.EstimationOutput>(elist.size());
	for (Estimation.EstimationOutput ein: elist)
	{
	    Messages.EstimationOutput.Builder builder = Messages.EstimationOutput.newBuilder()
		.setTime(ein.time.durationFrom(AbsoluteDate.J2000_EPOCH))
		.addAllEstimatedState(DoubleStream.of(ein.estimatedState).boxed().collect(Collectors.toList()));
	    if (ein.station != null && ein.station.length() > 0)
		builder = builder.setStation(ein.station);
	    if (ein.propagatedCovariance != null)
		builder = builder.addAllPropagatedCovariance(DoubleStream.of(ein.propagatedCovariance).boxed().collect(Collectors.toList()));
	    if (ein.innovationCovariance != null)
		builder = builder.addAllInnovationCovariance(DoubleStream.of(ein.innovationCovariance).boxed().collect(Collectors.toList()));
	    if (ein.estimatedCovariance != null)
		builder = builder.addAllEstimatedCovariance(DoubleStream.of(ein.estimatedCovariance).boxed().collect(Collectors.toList()));
	    if (ein.preFit != null)
		builder = builder.addAllPreFit(DoubleStream.of(ein.preFit).boxed().collect(Collectors.toList()));
	    if (ein.postFit != null)
		builder = builder.addAllPostFit(DoubleStream.of(ein.postFit).boxed().collect(Collectors.toList()));
	    if (ein.clutterProbability != null)
		builder = builder.setClutterProbability(ein.clutterProbability);
	    output.add(builder.build());
	}

	return(output);
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
