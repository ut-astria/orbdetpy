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
import org.orekit.estimation.measurements.modifiers.Bias;
import org.orekit.estimation.measurements.modifiers.OutlierFilter;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.PVCoordinates;

public final class Measurements
{
    public static class KeplerianElements
    {
	public double sma;
	public double ecc;
	public double inc;
	public double raan;
	public double argP;
	public double meanAnom;

	public KeplerianElements()
	{
	}

	public KeplerianElements(double a, double e, double i, double W, double o, double M)
	{
	    this.sma = a;
	    this.ecc = e;
	    this.inc = i;
	    this.raan = W;
	    this.argP = o;
	    this.meanAnom = M;
	}
    }

    public static class EquinoctialElements
    {
	public double sma;
	public double ex;
	public double ey;
	public double hx;
	public double hy;
	public double lm;

	public EquinoctialElements()
	{
	}

	public EquinoctialElements(double a, double ex, double ey, double hx, double hy, double lm)
	{
	    this.sma = a;
	    this.ex = ex;
	    this.ey = ey;
	    this.hx = hx;
	    this.hy = hy;
	    this.lm = lm;
	}
    }

    public static class State
    {
	public double[] cartesian;
	public KeplerianElements keplerian;
	public EquinoctialElements equinoctial;

	public State()
	{
	}
    }

    public static class Measurement
    {
	public String time;
	public String station;
	public double azimuth;
	public double elevation;
	public double range;
	public double rangeRate;
	public double rightAscension;
	public double declination;
	public double[] position;
	public double[] positionVelocity;
	public double[] angleRates;

	public Measurement()
	{
	}

	public Measurement(Measurement src)
	{
	    this.time = src.time;
	    this.station = src.station;
	    this.azimuth = src.azimuth;
	    this.elevation = src.elevation;
	    this.range = src.range;
	    this.rangeRate = src.rangeRate;
	    this.rightAscension = src.rightAscension;
	    this.declination = src.declination;
	    this.position = src.position;
	    this.positionVelocity = src.positionVelocity;
	    this.angleRates = src.angleRates;
	}
    }

    public static class SimulatedMeasurement extends Measurement
    {
	public State trueState;
	public double atmDensity;
	public double[] accGravity;
	public double[] accDrag;
	public double[] accOceanTides;
	public double[] accSolidTides;
	public double[] accThirdBodies;
	public double[] accRadiationPressure;
	public double[] accThrust;
	public double[] stationState;

	public SimulatedMeasurement()
	{
	}

	public SimulatedMeasurement(SimulatedMeasurement src)
	{
	    super(src);
	    this.trueState = src.trueState;
	    this.atmDensity = src.atmDensity;
	    this.accGravity = src.accGravity;
	    this.accDrag = src.accDrag;
	    this.accOceanTides = src.accOceanTides;
	    this.accSolidTides = src.accSolidTides;
	    this.accThirdBodies = src.accThirdBodies;
	    this.accRadiationPressure = src.accRadiationPressure;
	    this.accThrust = src.accThrust;
	    this.stationState = src.stationState;
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
	ArrayList<Measurement> tempraw = new ArrayList<Measurement>(rawMeas.length);
	for (Measurement m: rawMeas)
	{
	    if (m.station != null || m.position != null || m.positionVelocity != null)
		tempraw.add(m);
	}
	rawMeas = tempraw.toArray(new Measurement[0]);

	measObjs = new ArrayList<ObservedMeasurement<?>>(rawMeas.length);
	final Settings.Measurement cazim = odCfg.cfgMeasurements.get("azimuth");
	final Settings.Measurement celev = odCfg.cfgMeasurements.get("elevation");
	final Settings.Measurement crigh = odCfg.cfgMeasurements.get("rightAscension");
	final Settings.Measurement cdecl = odCfg.cfgMeasurements.get("declination");
	final Settings.Measurement crang = odCfg.cfgMeasurements.get("range");
	final Settings.Measurement crrat = odCfg.cfgMeasurements.get("rangeRate");
	final Settings.Measurement cpos = odCfg.cfgMeasurements.get("position");
	final Settings.Measurement cposvel = odCfg.cfgMeasurements.get("positionVelocity");
	final OutlierFilter outlier = new OutlierFilter(odCfg.estmOutlierWarmup, odCfg.estmOutlierSigma);
	final boolean addBias = odCfg.estmFilter.equalsIgnoreCase("EKF");
	final boolean addOutlier = addBias && odCfg.estmOutlierSigma > 0.0 && odCfg.estmOutlierWarmup > 0;
	final ObservableSatellite satellite = new ObservableSatellite(0);
	final double[] oneOnes = new double[] {1.0};
	final double[] twoOnes = new double[] {1.0, 1.0};
	final double[] oneNegInf = new double[] {Double.NEGATIVE_INFINITY};
	final double[] onePosInf = new double[] {Double.POSITIVE_INFINITY};
	final double[] twoNegInf = new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY};
	final double[] twoPosInf = new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
	final String[] biasAzEl = new String[] {"Az", "El"};
	final String[] biasRaDec = new String[] {"RA", "Dec"};
	final String[] biasRange = new String[] {"range"};
	final String[] biasRangeRate = new String[] {"rangeRate"};

	for (Measurement m: rawMeas)
	{
	    GroundStation gs = null;
	    Settings.Station jsn = null;
	    AbsoluteDate time = DataManager.parseDateTime(m.time);
	    if (m.station != null)
	    {
		gs = odCfg.stations.get(m.station);
		jsn = odCfg.cfgStations.get(m.station);
	    }

	    if (m.azimuth != 0.0 && cazim != null && celev != null)
	    {
		AngularAzEl obs = new AngularAzEl(gs, time, new double[]{m.azimuth, m.elevation},
						  new double[] {cazim.error[0], celev.error[0]}, twoOnes, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && (jsn.azimuthBias != 0.0 || jsn.elevationBias != 0.0))
		    obs.addModifier(new Bias<AngularAzEl>(biasAzEl, new double[] {jsn.azimuthBias, jsn.elevationBias},
							  twoOnes, twoNegInf, twoPosInf));
		measObjs.add(obs);
	    }

	    if (m.rightAscension != 0.0 && crigh != null && cdecl != null)
	    {
		AngularRaDec obs = new AngularRaDec(gs, FramesFactory.getFrame(Predefined.EME2000), time,
						    new double[] {m.rightAscension, m.declination},
						    new double[]{crigh.error[0], cdecl.error[0]}, twoOnes, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && (jsn.rightAscensionBias != 0.0 || jsn.declinationBias != 0.0))
		    obs.addModifier(new Bias<AngularRaDec>(biasRaDec, new double[] {jsn.rightAscensionBias, jsn.declinationBias},
							   twoOnes, twoNegInf, twoPosInf));
		measObjs.add(obs);
	    }

	    if (m.range != 0.0 && crang != null)
	    {
		Range obs = new Range(gs, crang.twoWay, time, m.range, crang.error[0], 1.0, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && jsn.rangeBias != 0.0)
		    obs.addModifier(new Bias<Range>(biasRange, new double[] {jsn.rangeBias}, oneOnes, oneNegInf, onePosInf));
		measObjs.add(obs);
	    }

	    if (m.rangeRate != 0.0 && crrat != null)
	    {
		RangeRate obs = new RangeRate(gs, time, m.rangeRate, crrat.error[0], 1.0, crrat.twoWay, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && jsn.rangeRateBias != 0.0)
		    obs.addModifier(new Bias<RangeRate>(biasRangeRate, new double[] {jsn.rangeRateBias}, oneOnes, oneNegInf, onePosInf));
		measObjs.add(obs);
	    }

	    if (m.position != null && cpos != null)
	    {
		double[] X = m.position;
		Position obs = new Position(time, new Vector3D(X[0], X[1], X[2]), cpos.error, 1.0, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && jsn != null && jsn.positionBias != null)
		    obs.addModifier(new Bias<Position>(
					new String[] {"x", "y", "z"}, jsn.positionBias,	new double[] {1.0, 1.0, 1.0},
					new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
					new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}));
		measObjs.add(obs);
	    }

	    if (m.positionVelocity != null && cposvel != null)
	    {
		double[] X = m.positionVelocity;
		PV obs = new PV(time, new Vector3D(X[0], X[1], X[2]), new Vector3D(X[3], X[4], X[5]), cposvel.error, 1.0, satellite);
		if (addOutlier)
		    obs.addModifier(outlier);
		if (addBias && jsn != null && jsn.positionVelocityBias != null)
		    obs.addModifier(new Bias<PV>(
					new String[] {"x", "y", "z", "Vx", "Vy", "Vz"}, jsn.positionVelocityBias,
					new double[] {1.0, 1.0, 1.0, 1.0, 1.0, 1.0},
					new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY,
						      Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
					new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY,
						      Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}));
		measObjs.add(obs);
	    }
	}
    }
}
