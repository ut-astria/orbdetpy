/*
 * Measurements.java - Functions to parse OD measurement files.
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

import java.util.Arrays;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.Position;
import org.orekit.estimation.measurements.PV;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.measurements.RangeRate;
import org.orekit.estimation.measurements.modifiers.AngularRadioRefractionModifier;
import org.orekit.estimation.measurements.modifiers.Bias;
import org.orekit.estimation.measurements.modifiers.OutlierFilter;
import org.orekit.models.earth.EarthITU453AtmosphereRefraction;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.PVCoordinates;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class Measurements
{
    public static class Measurement
    {
	public AbsoluteDate time;
	public String station;
	public double[] values;
	public double[] angleRates;
	public double[] trueState;
	public ObservedMeasurement<?>[] helpers;

	public Measurement()
	{
	}

	public Measurement(Measurement src)
	{
	    this.time = src.time;
	    this.station = src.station;
	    this.values = src.values;
	    this.angleRates = src.angleRates;
	    this.trueState = src.trueState;
	    this.helpers = src.helpers;
	}

	public Measurement(TimeStampedPVCoordinates pv, double[] extras)
	{
	    this.time = pv.getDate();
	    double[] p = pv.getPosition().toArray();
	    double[] v = pv.getVelocity().toArray();
	    if (extras != null && extras.length > 0)
		this.trueState = new double[6 + extras.length];
	    else
		this.trueState = new double[6];
	    for (int i = 0; i < this.trueState.length; i++)
	    {
		if (i < 3)
		    this.trueState[i] = p[i];
		else if (i < 6)
		    this.trueState[i] = v[i - 3];
		else
		    this.trueState[i] = extras[i - 6];
	    }
	}

	public void addHelper(ObservedMeasurement<?> om)
	{
	    if (helpers == null)
		helpers = new ObservedMeasurement<?>[]{om};
	    else
	    {
		helpers = Arrays.copyOf(helpers, helpers.length + 1);
		helpers[helpers.length - 1] = om;
	    }
	}
    }

    public Measurement[] array;

    public Measurements build(Settings odCfg)
    {
	if (array.length == 0)
	    throw(new RuntimeException("No measurements provided"));
	Settings.Measurement cazim = odCfg.cfgMeasurements.get(Settings.MeasurementType.AZIMUTH);
	Settings.Measurement celev = odCfg.cfgMeasurements.get(Settings.MeasurementType.ELEVATION);
	Settings.Measurement crigh = odCfg.cfgMeasurements.get(Settings.MeasurementType.RIGHT_ASCENSION);
	Settings.Measurement cdecl = odCfg.cfgMeasurements.get(Settings.MeasurementType.DECLINATION);
	Settings.Measurement crang = odCfg.cfgMeasurements.get(Settings.MeasurementType.RANGE);
	Settings.Measurement crrat = odCfg.cfgMeasurements.get(Settings.MeasurementType.RANGE_RATE);
	Settings.Measurement cpos = odCfg.cfgMeasurements.get(Settings.MeasurementType.POSITION);
	Settings.Measurement cposvel = odCfg.cfgMeasurements.get(Settings.MeasurementType.POSITION_VELOCITY);
	if (cazim == null && celev == null && crigh == null && cdecl == null && crang == null && crrat == null &&
	    cpos == null && cposvel == null)
	    throw(new RuntimeException("Measurement types not defined"));

	boolean addBias = odCfg.estmFilter == Settings.Filter.EXTENDED_KALMAN;
	boolean addOutlier = addBias && odCfg.estmOutlierSigma > 0.0 && odCfg.estmOutlierWarmup > 0;
	ObservableSatellite satellite = new ObservableSatellite(0);
	double[] oneOnes = {1.0};
	double[] twoOnes = {1.0, 1.0};
	double[] oneNegInf = {Double.NEGATIVE_INFINITY};
	double[] onePosInf = {Double.POSITIVE_INFINITY};
	double[] twoNegInf = {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
	double[] twoPosInf = {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
	String[] biasAzEl = {Settings.MeasurementType.AZIMUTH.name(), Settings.MeasurementType.ELEVATION.name()};
	String[] biasRaDec = {Settings.MeasurementType.RIGHT_ASCENSION.name(), Settings.MeasurementType.DECLINATION.name()};
	String[] biasRange = {Settings.MeasurementType.RANGE.name()};
	String[] biasRangeRate = {Settings.MeasurementType.RANGE_RATE.name()};

	for (Measurement m: array)
	{
	    if (m.values.length == 0)
		continue;
	    GroundStation gs = null;
	    Settings.Station jsn = null;
	    if (m.station.length() > 0)
	    {
		gs = odCfg.stations.get(m.station);
		jsn = odCfg.cfgStations.get(m.station);
	    }

	    if (cazim != null && celev != null)
	    {
		AngularAzEl obs = new AngularAzEl(gs, m.time, new double[]{m.values[0], m.values[1]},
						  new double[]{cazim.error[0], celev.error[0]}, twoOnes, satellite);
		m.addHelper(obs);
		obs.addModifier(new AngularRadioRefractionModifier(new EarthITU453AtmosphereRefraction(jsn.altitude)));
		if (addOutlier)
		    obs.addModifier(new OutlierFilter<AngularAzEl>(odCfg.estmOutlierWarmup, odCfg.estmOutlierSigma));
		if (addBias && jsn.bias != null && jsn.bias.length > 0)
		    obs.addModifier(new Bias<AngularAzEl>(biasAzEl, new double[]{jsn.bias[0], jsn.bias[1]}, twoOnes, twoNegInf, twoPosInf));
	    }

	    if (crigh != null && cdecl != null)
	    {
		AngularRaDec obs = new AngularRaDec(gs, odCfg.propInertialFrame, m.time, new double[]{m.values[0], m.values[1]},
						    new double[]{crigh.error[0], cdecl.error[0]}, twoOnes, satellite);
		m.addHelper(obs);
		if (addOutlier)
		    obs.addModifier(new OutlierFilter<AngularRaDec>(odCfg.estmOutlierWarmup, odCfg.estmOutlierSigma));
		if (addBias && jsn.bias != null && jsn.bias.length > 0)
		    obs.addModifier(new Bias<AngularRaDec>(biasRaDec, new double[]{jsn.bias[0], jsn.bias[1]}, twoOnes, twoNegInf, twoPosInf));
	    }

	    if (crang != null)
	    {
		Range obs = new Range(gs, crang.twoWay, m.time, m.values[0], crang.error[0], 1.0, satellite);
		m.addHelper(obs);
		if (addOutlier)
		    obs.addModifier(new OutlierFilter<Range>(odCfg.estmOutlierWarmup, odCfg.estmOutlierSigma));
		if (addBias && jsn.bias != null && jsn.bias.length > 0)
		    obs.addModifier(new Bias<Range>(biasRange, new double[]{jsn.bias[0]}, oneOnes, oneNegInf, onePosInf));
	    }

	    if (crrat != null)
	    {
		RangeRate obs = new RangeRate(gs, m.time, m.values[m.values.length-1], crrat.error[0], 1.0, crrat.twoWay, satellite);
		m.addHelper(obs);
		if (addOutlier)
		    obs.addModifier(new OutlierFilter<RangeRate>(odCfg.estmOutlierWarmup, odCfg.estmOutlierSigma));
		if (addBias && jsn.bias != null && jsn.bias.length > 0)
		    obs.addModifier(new Bias<RangeRate>(biasRangeRate, new double[]{jsn.bias[jsn.bias.length-1]}, oneOnes, oneNegInf, onePosInf));
	    }

	    if (cpos != null)
	    {
		Position obs = new Position(m.time, new Vector3D(m.values[0], m.values[1], m.values[2]), cpos.error, 1.0, satellite);
		m.addHelper(obs);
		if (addOutlier)
		    obs.addModifier(new OutlierFilter<Position>(odCfg.estmOutlierWarmup, odCfg.estmOutlierSigma));
	    }

	    if (cposvel != null)
	    {
		PV obs = new PV(m.time, new Vector3D(m.values[0], m.values[1], m.values[2]),
				new Vector3D(m.values[3], m.values[4], m.values[5]), cposvel.error, 1.0, satellite);
		m.addHelper(obs);
		if (addOutlier)
		    obs.addModifier(new OutlierFilter<PV>(odCfg.estmOutlierWarmup, odCfg.estmOutlierSigma));
	    }
	}
	return(this);
    }
}
