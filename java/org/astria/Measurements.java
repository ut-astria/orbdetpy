/*
 * Measurements.java - Functions to parse OD measurement files.
 * Copyright (C) 2018 Shiva Iyer <shiva.iyer AT utexas DOT edu>
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
import org.orekit.estimation.measurements.AngularAzEl;
import org.orekit.estimation.measurements.AngularRaDec;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.estimation.measurements.ObservedMeasurement;
import org.orekit.estimation.measurements.Range;
import org.orekit.estimation.measurements.RangeRate;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;

public class Measurements
{
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
	Double[] TrueState;
    }

    JSONMeasurement[] rawmeas;
    ArrayList<ObservedMeasurement<?>> measobjs;

    public static Measurements loadJSON(String json, HashMap<String, GroundStation> gsta,
					Map<String, Settings.JSONMeasurement> mcfg)
    {
	Measurements m = new Measurements();
	m.rawmeas = new Gson().fromJson(json, JSONMeasurement[].class);
	m.getMeasurementObjects(gsta, mcfg);
	return(m);
    }

    private void getMeasurementObjects(HashMap<String, GroundStation> gsta,
				       Map<String, Settings.JSONMeasurement> mcfg)
    {
	measobjs = new ArrayList<ObservedMeasurement<?>>();
	for (JSONMeasurement m: rawmeas)
	{
	    GroundStation gs = gsta.get(m.Station);
	    AbsoluteDate time = new AbsoluteDate(DateTimeComponents.parseDateTime(m.Time),
						 DataManager.utcscale);

	    if (m.Azimuth != null)
	    {
                measobjs.add(new AngularAzEl(gs, time, new double[]{m.Azimuth, m.Elevation},
					     new double[]{mcfg.get("Azimuth").Error,
							  mcfg.get("Elevation").Error},
					     new double[]{1.0, 1.0}));
	    }

	    if (m.RightAscension != null)
	    {
                measobjs.add(new AngularRaDec(gs, DataManager.eme2000, time,
					      new double[]{m.RightAscension, m.Declination},
					      new double[]{mcfg.get("RightAscension").Error,
							   mcfg.get("Declination").Error},
					      new double[]{1.0, 1.0}));
	    }

	    if (m.Range != null)
	    {
		Settings.JSONMeasurement c = mcfg.get("Range");
		measobjs.add(new Range(gs, time, m.Range, c.Error, 1.0, c.TwoWay));
	    }

	    if (m.RangeRate != null)
	    {
		Settings.JSONMeasurement c = mcfg.get("RangeRate");
		measobjs.add(new RangeRate(gs, time, m.RangeRate, c.Error, 1.0, c.TwoWay));
	    }
	}
    }
}
