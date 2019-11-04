/*
 * Measurements.java - Functions to parse OD measurement files.
 * Copyright (C) 2018-2019 University of Texas
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
import org.orekit.frames.Frame;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.PVCoordinates;

public final class Measurements
{
    public static class KeplerianElements
    {
	public double SMA;
	public double Ecc;
	public double Inc;
	public double RAAN;
	public double ArgP;
	public double MeanAnom;

	public KeplerianElements(double a, double e, double i, double W, double o, double M)
	{
	    this.SMA = a;
	    this.Ecc = e;
	    this.Inc = i;
	    this.RAAN = W;
	    this.ArgP = o;
	    this.MeanAnom = M;
	}
    }

    public static class EquinoctialElements
    {
	public double SMA;
	public double Ex;
	public double Ey;
	public double Hx;
	public double Hy;
	public double Lm;

	public EquinoctialElements(double a, double ex, double ey, double hx, double hy, double lm)
	{
	    this.SMA = a;
	    this.Ex = ex;
	    this.Ey = ey;
	    this.Hx = hx;
	    this.Hy = hy;
	    this.Lm = lm;
	}
    }

    public static class State
    {
	public double[] Cartesian;
	public KeplerianElements Kepler;
	public EquinoctialElements Equinoctial;
    }

    public static class Measurement
    {
	public String Time = null;
	public String Station = null;
	public Double Azimuth = null;
	public Double Elevation = null;
	public Double Range = null;
	public Double RangeRate = null;
	public Double RightAscension = null;
	public Double Declination = null;
	public Double[] Position = null;
	public Double[] PositionVelocity = null;

	public Measurement()
	{
	}

	public Measurement(Measurement src)
	{
	    this.Time = src.Time;
	    this.Station = src.Station;
	    this.Azimuth = src.Azimuth;
	    this.Elevation = src.Elevation;
	    this.Range = src.Range;
	    this.RangeRate = src.RangeRate;
	    this.RightAscension = src.RightAscension;
	    this.Declination = src.Declination;
	    this.Position = src.Position;
	    this.PositionVelocity = src.PositionVelocity;
	}
    }

    public static class SimulatedMeasurement extends Measurement
    {
	public State TrueState;
	public Double AtmDensity;
	public double[] AccGravity;
	public double[] AccDrag;
	public double[] AccOceanTides;
	public double[] AccSolidTides;
	public double[] AccThirdBodies;
	public double[] AccRadiationPressure;
	public double[] AccThrust;
	public double[] StationState;

	public SimulatedMeasurement()
	{
	}

	public SimulatedMeasurement(SimulatedMeasurement src)
	{
	    super(src);
	    this.TrueState = src.TrueState;
	    this.AtmDensity = src.AtmDensity;
	    this.AccGravity = src.AccGravity;
	    this.AccDrag = src.AccDrag;
	    this.AccOceanTides = src.AccOceanTides;
	    this.AccSolidTides = src.AccSolidTides;
	    this.AccThirdBodies = src.AccThirdBodies;
	    this.AccRadiationPressure = src.AccRadiationPressure;
	    this.AccThrust = src.AccThrust;
	    this.StationState = src.StationState;
	}
    }

    public Measurement[] rawMeas = null;
    public ArrayList<ObservedMeasurement<?>> measObjs = null;

    public Measurements build(Settings odCfg)
    {
	buildMeasurementObjects(odCfg);
	return(this);
    }

    private void buildMeasurementObjects(Settings odCfg)
    {
	Settings.Measurement cazim = odCfg.cfgMeasurements.get("Azimuth");
	Settings.Measurement celev = odCfg.cfgMeasurements.get("Elevation");
	Settings.Measurement crigh = odCfg.cfgMeasurements.get("RightAscension");
	Settings.Measurement cdecl = odCfg.cfgMeasurements.get("Declination");
	Settings.Measurement crang = odCfg.cfgMeasurements.get("Range");
	Settings.Measurement crrat = odCfg.cfgMeasurements.get("RangeRate");
	Settings.Measurement cpos = odCfg.cfgMeasurements.get("Position");
	Settings.Measurement cposvel = odCfg.cfgMeasurements.get("PositionVelocity");

	ArrayList<Measurement> tempraw = new ArrayList<Measurement>(rawMeas.length);
	for (Measurement m: rawMeas)
	{
	    if (m.Station != null || m.Position != null || m.PositionVelocity != null)
		tempraw.add(m);
	}
	rawMeas = tempraw.toArray(new Measurement[0]);

	measObjs = new ArrayList<ObservedMeasurement<?>>();
	for (Measurement m: rawMeas)
	{
	    GroundStation gs = null;
	    Settings.Station jsn = null;
	    if (m.Station != null)
	    {
		gs = odCfg.stations.get(m.Station);
		jsn = odCfg.cfgStations.get(m.Station);
	    }

	    AbsoluteDate time = new AbsoluteDate(DateTimeComponents.parseDateTime(m.Time), DataManager.getTimeScale("UTC"));
	    if (m.Azimuth != null && cazim != null && celev != null)
	    {
		AngularAzEl obs = new AngularAzEl(gs, time, new double[]{m.Azimuth, m.Elevation},
						  new double[]{cazim.error[0], celev.error[0]},
						  new double[]{1.0, 1.0}, new ObservableSatellite(0));
		if (jsn.azimuthBias != 0.0 || jsn.elevationBias != 0.0)
		    obs.addModifier(new Bias<AngularAzEl>(
					new String[] {"Az", "El"}, new double[] {jsn.azimuthBias, jsn.elevationBias},
					new double[] {1.0, 1.0}, new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
					new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}));
		measObjs.add(obs);
	    }

	    if (m.RightAscension != null && crigh != null && cdecl != null)
	    {
		AngularRaDec obs = new AngularRaDec(gs, DataManager.getFrame("EME2000"), time,
						    new double[]{m.RightAscension, m.Declination},
						    new double[]{crigh.error[0], cdecl.error[0]},
						    new double[]{1.0, 1.0}, new ObservableSatellite(0));
		if (jsn.rightAscensionBias != 0.0 || jsn.declinationBias != 0.0)
		    obs.addModifier(new Bias<AngularRaDec>(
					new String[] {"RA", "Dec"}, new double[] {jsn.rightAscensionBias, jsn.declinationBias},
					new double[] {1.0, 1.0}, new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
					new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}));
		measObjs.add(obs);
	    }

	    if (m.Range != null && crang != null)
	    {
		Range obs = new Range(gs, crang.twoWay, time, m.Range, crang.error[0], 1.0, new ObservableSatellite(0));
		if (jsn.rangeBias != 0.0)
		    obs.addModifier(new Bias<Range>(
					new String[] {"Range"}, new double[] {jsn.rangeBias}, new double[] {1.0},
					new double[] {Double.NEGATIVE_INFINITY}, new double[] {Double.POSITIVE_INFINITY}));
		measObjs.add(obs);
	    }

	    if (m.RangeRate != null && crrat != null)
	    {
		RangeRate obs = new RangeRate(gs, time, m.RangeRate, crrat.error[0], 1.0, crrat.twoWay, new ObservableSatellite(0));
		if (jsn.rangeRateBias != 0.0)
		    obs.addModifier(new Bias<RangeRate>(
					new String[] {"RangeRate"}, new double[] {jsn.rangeRateBias}, new double[] {1.0},
					new double[] {Double.NEGATIVE_INFINITY}, new double[] {Double.POSITIVE_INFINITY}));
		measObjs.add(obs);
	    }

	    if (m.Position != null && cpos != null)
	    {
		Double[] X = m.Position;
		Position obs = new Position(time, new Vector3D(X[0], X[1], X[2]), cpos.error, 1.0, new ObservableSatellite(0));
		if (jsn != null && jsn.positionBias != null)
		    obs.addModifier(new Bias<Position>(
					new String[] {"x", "y", "z"}, jsn.positionBias,	new double[] {1.0, 1.0, 1.0},
					new double[] {Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY, Double.NEGATIVE_INFINITY},
					new double[] {Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY}));
		measObjs.add(obs);
	    }

	    if (m.PositionVelocity != null && cposvel != null)
	    {
		Double[] X = m.PositionVelocity;
		PV obs = new PV(time, new Vector3D(X[0], X[1], X[2]), new Vector3D(X[3], X[4], X[5]),
				cposvel.error, 1.0, new ObservableSatellite(0));
		if (jsn != null && jsn.positionVelocityBias != null)
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
