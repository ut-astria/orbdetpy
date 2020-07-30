/*
 * Measurements.java - Functions to parse OD measurement files.
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
import java.util.HashMap;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.estimation.measurements.AbstractMeasurement;
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
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.frames.Transform;
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
	}

	public Measurement(TimeStampedPVCoordinates pv)
	{
	    this.time = pv.getDate();
	    Vector3D p = pv.getPosition();
	    Vector3D v = pv.getVelocity();
	    this.trueState = new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()};
	}
    }

    public Measurement[] rawMeas;
    public ArrayList<ObservedMeasurement<?>> measObjs;

    public Measurements build(Settings odCfg)
    {
	buildMeasurementObjects(odCfg);
	return(this);
    }

    private void buildMeasurementObjects(Settings odCfg)
    {
	ArrayList<Measurement> tempRaw = new ArrayList<Measurement>(rawMeas.length);
	for (Measurement m: rawMeas)
	{
	    if (m.station.length() > 0 || m.values.length >= 3)
		tempRaw.add(m);
	}
	rawMeas = tempRaw.toArray(new Measurement[0]);

	measObjs = new ArrayList<ObservedMeasurement<?>>(rawMeas.length);
	final Settings.Measurement cazim = odCfg.cfgMeasurements.get(Settings.MeasurementType.AZIMUTH);
	final Settings.Measurement celev = odCfg.cfgMeasurements.get(Settings.MeasurementType.ELEVATION);
	final Settings.Measurement crigh = odCfg.cfgMeasurements.get(Settings.MeasurementType.RIGHT_ASCENSION);
	final Settings.Measurement cdecl = odCfg.cfgMeasurements.get(Settings.MeasurementType.DECLINATION);
	final Settings.Measurement crang = odCfg.cfgMeasurements.get(Settings.MeasurementType.RANGE);
	final Settings.Measurement crrat = odCfg.cfgMeasurements.get(Settings.MeasurementType.RANGE_RATE);
	final Settings.Measurement cpos = odCfg.cfgMeasurements.get(Settings.MeasurementType.POSITION);
	final Settings.Measurement cposvel = odCfg.cfgMeasurements.get(Settings.MeasurementType.POSITION_VELOCITY);
	final OutlierFilter outlier = new OutlierFilter(odCfg.estmOutlierWarmup, odCfg.estmOutlierSigma);
	final boolean addBias = odCfg.estmFilter == Settings.Filter.EXTENDED_KALMAN;
	final boolean addOutlier = addBias && odCfg.estmOutlierSigma > 0.0 && odCfg.estmOutlierWarmup > 0;
	final ObservableSatellite satellite = new ObservableSatellite(0);
	final double[] oneOnes = new double[]{1.0};
	final double[] twoOnes = new double[]{1.0, 1.0};
	final double[] oneNegInf = new double[]{Double.NEGATIVE_INFINITY};
	final double[] onePosInf = new double[]{Double.POSITIVE_INFINITY};
	final double[] twoNegInf = new double[]{Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
	final double[] twoPosInf = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
	final String[] biasAzEl = new String[]{Settings.MeasurementType.AZIMUTH.name(), Settings.MeasurementType.ELEVATION.name()};
	final String[] biasRaDec = new String[]{Settings.MeasurementType.RIGHT_ASCENSION.name(), Settings.MeasurementType.DECLINATION.name()};
	final String[] biasRange = new String[]{Settings.MeasurementType.RANGE.name()};
	final String[] biasRangeRate = new String[]{Settings.MeasurementType.RANGE_RATE.name()};

	for (Measurement m: rawMeas)
	{
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
		obs.addModifier(new AngularRadioRefractionModifier(new EarthITU453AtmosphereRefraction(jsn.altitude)));
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && jsn.bias != null && jsn.bias.length > 0)
		    obs.addModifier(new Bias<AngularAzEl>(biasAzEl, new double[]{jsn.bias[0], jsn.bias[1]}, twoOnes, twoNegInf, twoPosInf));
		measObjs.add(obs);
	    }

	    if (crigh != null && cdecl != null)
	    {
		AngularRaDec obs = new AngularRaDec(gs, odCfg.propInertialFrame, m.time, new double[]{m.values[0], m.values[1]},
						    new double[]{crigh.error[0], cdecl.error[0]}, twoOnes, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && jsn.bias != null && jsn.bias.length > 0)
		    obs.addModifier(new Bias<AngularRaDec>(biasRaDec, new double[]{jsn.bias[0], jsn.bias[1]}, twoOnes, twoNegInf, twoPosInf));
		measObjs.add(obs);
	    }

	    if (crang != null)
	    {
		Range obs = new Range(gs, crang.twoWay, m.time, m.values[0], crang.error[0], 1.0, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && jsn.bias != null && jsn.bias.length > 0)
		    obs.addModifier(new Bias<Range>(biasRange, new double[]{jsn.bias[0]}, oneOnes, oneNegInf, onePosInf));
		measObjs.add(obs);
	    }

	    if (crrat != null)
	    {
		RangeRate obs = new RangeRate(gs, m.time, m.values[m.values.length-1], crrat.error[0], 1.0, crrat.twoWay, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && jsn.bias != null && jsn.bias.length > 0)
		    obs.addModifier(new Bias<RangeRate>(biasRangeRate, new double[]{jsn.bias[jsn.bias.length-1]}, oneOnes, oneNegInf, onePosInf));
		measObjs.add(obs);
	    }

	    if (cpos != null)
	    {
		Position obs = new Position(m.time, new Vector3D(m.values[0], m.values[1], m.values[2]), cpos.error, 1.0, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		measObjs.add(obs);
	    }

	    if (cposvel != null)
	    {
		PV obs = new PV(m.time, new Vector3D(m.values[0], m.values[1], m.values[2]),
				new Vector3D(m.values[3], m.values[4], m.values[5]), cposvel.error, 1.0, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		measObjs.add(obs);
	    }
	}
    }
}
