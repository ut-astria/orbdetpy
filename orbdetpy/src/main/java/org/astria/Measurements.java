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
import org.orekit.frames.Frame;
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
	    if (m.station != null || m.position != null || m.positionVelocity != null)
		tempraw.add(m);
	}
	rawMeas = tempraw.toArray(new Measurement[0]);

	measObjs = new ArrayList<ObservedMeasurement<?>>(rawMeas.length);
	for (Measurement m: rawMeas)
	{
	    GroundStation gs = null;
	    Settings.Station jsn = null;
	    if (m.station != null)
	    {
		gs = odCfg.stations.get(m.station);
		jsn = odCfg.cfgStations.get(m.station);
	    }

	    AbsoluteDate time = DataManager.parseDateTime(m.time);
	    if (m.azimuth != 0.0 && cazim != null && celev != null)
		measObjs.add(new AngularAzEl(gs, time, new double[]{m.azimuth, m.elevation}, new double[]{cazim.error[0], celev.error[0]},
					     new double[]{1.0, 1.0}, new ObservableSatellite(0)));

	    if (m.rightAscension != 0.0 && crigh != null && cdecl != null)
		measObjs.add(new AngularRaDec(gs, DataManager.getFrame("EME2000"), time, new double[]{m.rightAscension, m.declination},
					      new double[]{crigh.error[0], cdecl.error[0]}, new double[]{1.0, 1.0}, new ObservableSatellite(0)));

	    if (m.range != 0.0 && crang != null)
		measObjs.add(new Range(gs, crang.twoWay, time, m.range, crang.error[0], 1.0, new ObservableSatellite(0)));

	    if (m.rangeRate != 0.0 && crrat != null)
		measObjs.add(new RangeRate(gs, time, m.rangeRate, crrat.error[0], 1.0, crrat.twoWay, new ObservableSatellite(0)));

	    if (m.position != null && cpos != null)
	    {
		double[] X = m.position;
		measObjs.add(new Position(time, new Vector3D(X[0], X[1], X[2]), cpos.error, 1.0, new ObservableSatellite(0)));
	    }

	    if (m.positionVelocity != null && cposvel != null)
	    {
		double[] X = m.positionVelocity;
		measObjs.add(new PV(time, new Vector3D(X[0], X[1], X[2]), new Vector3D(X[3], X[4], X[5]),
				    cposvel.error, 1.0, new ObservableSatellite(0)));
	    }
	}
    }
}
