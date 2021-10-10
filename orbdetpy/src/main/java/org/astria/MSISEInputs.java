/*
 * MSISEInputs.java - Class for reading MSISE space weather data.
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

import java.io.File;
import java.util.HashMap;
import java.util.Scanner;
import org.orekit.models.earth.atmosphere.NRLMSISE00InputParameters;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;

public final class MSISEInputs implements NRLMSISE00InputParameters
{
    private static MSISEInputs singleton;
    private AbsoluteDate minDate;
    private AbsoluteDate maxDate;
    private HashMap<String, double[]> weather;

    private MSISEInputs() throws Exception
    {
	String[] tok = null;
	weather = new HashMap<String, double[]>();
	Scanner scan = new Scanner(new File(DataManager.dataPath, "SpaceWeather.dat"));
	while (scan.hasNextLine())
	{
	    tok = scan.nextLine().split(",");
	    double[] vals = new double[tok.length];
	    for (int i = 0; i < tok.length; i++)
	    {
		if (tok[i].length() > 0)
		    vals[i] = Double.parseDouble(tok[i]);
	    }
	    weather.put(String.format("%s%s%s", tok[0], tok[1], tok[2]), vals);
	    if (minDate == null)
		minDate = new AbsoluteDate(Integer.parseInt(tok[0]), Integer.parseInt(tok[1]), Integer.parseInt(tok[2]), TimeScalesFactory.getUTC());
	}

	scan.close();
	if (tok != null)
	    maxDate = new AbsoluteDate(Integer.parseInt(tok[0]), Integer.parseInt(tok[1]), Integer.parseInt(tok[2]), TimeScalesFactory.getUTC());
    }

    public static synchronized MSISEInputs getInstance()
    {
	if (singleton == null)
	{
	    try
	    {
		singleton = new MSISEInputs();
	    }
	    catch (Exception exc)
	    {
		throw(new RuntimeException(exc));
	    }
	}
	return(singleton);
    }

    @Override public AbsoluteDate getMinDate()
    {
	return(minDate);
    }

    @Override public AbsoluteDate getMaxDate()
    {
	return(maxDate);
    }

    @Override public double getDailyFlux(AbsoluteDate date)
    {
	String p = date.shiftedBy(-86400.0).toString();
	String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
	return(weather.get(k)[26]);
    }

    @Override public double getAverageFlux(AbsoluteDate date)
    {
	String p = date.toString();
	String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
	return(weather.get(k)[28]);
    }

    @Override public double[] getAp(AbsoluteDate date)
    {
	double[] apValues = new double[7];
	for (int i = 0; i < 7; i++)
	{
	    if (i == 0)
	    {
		String p = date.toString();
		String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
		apValues[0] = weather.get(k)[22];
	    }
	    else if (i <= 4)
	    {
		String p = date.shiftedBy(-10800.0*(i-1)).toString();
		String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
		apValues[i] = weather.get(k)[Integer.parseInt(p.substring(11, 13))/3 + 14];
	    }
	    else
	    {
		for (int j = 8*i - 36; j <= 8*i - 29; j++)
		{
		    String p = date.shiftedBy(-10800.0*(j-1)).toString();
		    String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
		    apValues[i] += weather.get(k)[Integer.parseInt(p.substring(11, 13))/3 + 14];
		}
		apValues[i] = apValues[i]/8.0;
	    }
	}
	return(apValues);
    }
}
