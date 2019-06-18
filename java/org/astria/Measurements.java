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

import com.google.gson.Gson;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import org.astria.DataManager;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservableSatellite;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.PV;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.measurements.RangeRate;
import org.orekit.frames.Frame;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.PVCoordinates;

public class Measurements
{
    class JSONKepler
    {
	double SMA;
	double Ecc;
	double Inc;
	double RAAN;
	double ArgP;
	double MeanAnom;

	public JSONKepler(double a, double e, double i, double W, double o, double M)
	{
	    this.SMA = a;
	    this.Ecc = e;
	    this.Inc = i;
	    this.RAAN = W;
	    this.ArgP = o;
	    this.MeanAnom = M;
	}
    }

    class JSONEquinoctial
    {
	double SMA;
	double Ex;
	double Ey;
	double Hx;
	double Hy;
	double Lm;

	public JSONEquinoctial(double a, double ex, double ey, double hx, double hy, double lm)
	{
	    this.SMA = a;
	    this.Ex = ex;
	    this.Ey = ey;
	    this.Hx = hx;
	    this.Hy = hy;
	    this.Lm = lm;
	}
    }

    class JSONState
    {
	double[] Cartesian;
	JSONKepler Kepler;
	JSONEquinoctial Equinoctial;
    }

    class JSONMeasurement
    {
	String Time;
	String Station;
	Double Azimuth;
	Double Elevation;
	Double Range;
	Double RangeRate;
	Double RightAscension;
	Double Declination;
	Double[] PositionVelocity;

	public JSONMeasurement()
	{
	}

	public JSONMeasurement(JSONMeasurement src)
	{
	    Time = src.Time;
	    Station = src.Station;
	    Azimuth = src.Azimuth;
	    Elevation = src.Elevation;
	    Range = src.Range;
	    RangeRate = src.RangeRate;
	    RightAscension = src.RightAscension;
	    Declination = src.Declination;
	    PositionVelocity = src.PositionVelocity;
	}
    }

    class JSONSimulatedMeasurement extends JSONMeasurement
    {
	JSONState TrueState;
	Double AtmDensity;
	double[] AccGravity;
	double[] AccDrag;
	double[] AccOceanTides;
	double[] AccSolidTides;
	double[] AccThirdBodies;
	double[] AccRadiationPressure;
	double[] AccThrust;

	public JSONSimulatedMeasurement()
	{
	}

	public JSONSimulatedMeasurement(JSONSimulatedMeasurement src)
	{
	    super(src);
	    TrueState = src.TrueState;
	    AtmDensity = src.AtmDensity;
	    AccGravity = src.AccGravity;
	    AccDrag = src.AccDrag;
	    AccOceanTides = src.AccOceanTides;
	    AccSolidTides = src.AccSolidTides;
	    AccThirdBodies = src.AccThirdBodies;
	    AccRadiationPressure = src.AccRadiationPressure;
	    AccThrust = src.AccThrust;
	}
    }

    JSONMeasurement[] rawmeas;
    ArrayList<ObservedMeasurement<?>> measobjs;

    public static Measurements loadJSON(Settings odcfg, String json)
    {
	Measurements m = new Measurements();
	m.rawmeas = new Gson().fromJson(json, JSONMeasurement[].class);
	m.getMeasurementObjects(odcfg.stations, odcfg.Measurements, odcfg.propframe);
	return(m);
    }

    private void getMeasurementObjects(HashMap<String, GroundStation> gsta,
				       Map<String, Settings.JSONMeasurement> mcfg, Frame propframe)
    {
	measobjs = new ArrayList<ObservedMeasurement<?>>();

	for (JSONMeasurement m: rawmeas)
	{
	    GroundStation gs = null;
	    if (m.Station != null)
		gs = gsta.get(m.Station);

	    AbsoluteDate time = new AbsoluteDate(DateTimeComponents.parseDateTime(m.Time),
						 DataManager.utcscale);

	    if (m.Azimuth != null)
	    {
                measobjs.add(new AngularAzEl(gs, time, new double[]{m.Azimuth, m.Elevation},
					     new double[]{mcfg.get("Azimuth").Error[0],
							  mcfg.get("Elevation").Error[0]},
					     new double[]{1.0, 1.0}, new ObservableSatellite(0)));
	    }

	    if (m.RightAscension != null)
	    {
                measobjs.add(new AngularRaDec(gs, DataManager.eme2000, time,
					      new double[]{m.RightAscension, m.Declination},
					      new double[]{mcfg.get("RightAscension").Error[0],
							   mcfg.get("Declination").Error[0]},
					      new double[]{1.0, 1.0}, new ObservableSatellite(0)));
	    }

	    if (m.Range != null)
	    {
		Settings.JSONMeasurement c = mcfg.get("Range");
		measobjs.add(new Range(gs, c.TwoWay, time, m.Range, c.Error[0], 1.0,
				       new ObservableSatellite(0)));
	    }

	    if (m.RangeRate != null)
	    {
		Settings.JSONMeasurement c = mcfg.get("RangeRate");
		measobjs.add(new RangeRate(gs, time, m.RangeRate, c.Error[0], 1.0, c.TwoWay,
					   new ObservableSatellite(0)));
	    }

	    if (m.PositionVelocity != null)
	    {
		Double[] X = m.PositionVelocity;
		Settings.JSONMeasurement c = mcfg.get("PositionVelocity");
		if (c.ReferenceFrame != null)
		{
		    Frame fromframe = DataManager.eme2000;
		    if (c.ReferenceFrame.equals("GCRF"))
			fromframe = DataManager.gcrf;
		    else if (c.ReferenceFrame.equals("ITRF"))
			fromframe = DataManager.itrf;

		    Transform xfm = fromframe.getTransformTo(propframe, time);
		    PVCoordinates frompv = new PVCoordinates(new Vector3D(X[0], X[1], X[2]),
							     new Vector3D(X[3], X[4], X[5]));
		    PVCoordinates topv = xfm.transformPVCoordinates(frompv);

		    Vector3D p = topv.getPosition();
		    Vector3D v = topv.getVelocity();
		    X = new Double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()};
		}

		measobjs.add(new PV(time, new Vector3D(X[0], X[1], X[2]),
				    new Vector3D(X[3], X[4], X[5]),
				    c.Error, 1.0, new ObservableSatellite(0)));
	    }
	}
    }
}
