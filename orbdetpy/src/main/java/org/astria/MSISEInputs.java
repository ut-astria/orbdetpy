/*
 * MSISEInputs.java - Class for reading MSISE space weather data.
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

import java.util.HashMap;
import org.orekit.models.earth.atmosphere.NRLMSISE00InputParameters;
import org.orekit.time.AbsoluteDate;

public final class MSISEInputs implements NRLMSISE00InputParameters
{
    private final AbsoluteDate minDate;
    private final AbsoluteDate maxDate;
    private final HashMap<String, double[]> swData;
    private final double[] apVals;
    
    public MSISEInputs(AbsoluteDate min, AbsoluteDate max, HashMap<String, double[]> sw, int apflag)
    {
	minDate = min;
	maxDate = max;
	swData = sw;
	if (apflag == 1)
	    apVals = new double[1];
	else
	    apVals = new double[7];
    }

    public AbsoluteDate getMinDate()
    {
	return(minDate);
    }

    public AbsoluteDate getMaxDate()
    {
	return(maxDate);
    }

    public double getDailyFlux(AbsoluteDate date)
    {
	String p = date.shiftedBy(-86400.0).toString();
	String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
	return(swData.get(k)[26]);
    }

    public double getAverageFlux(AbsoluteDate date)
    {
	String p = date.toString();
	String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
	return(swData.get(k)[28]);
    }

    public double[] getAp(AbsoluteDate date)
    {
	for (int i = 0; i < apVals.length; i++)
	{
	    if (i == 0)
	    {
		String p = date.toString();
		String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
		apVals[0] = swData.get(k)[22];
	    }
	    else if (i <= 4)
	    {
		String p = date.shiftedBy(-10800.0*(i-1)).toString();
		String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
		apVals[i] = swData.get(k)[Integer.parseInt(p.substring(11, 13))/3+14];
	    }
	    else
	    {
		apVals[i] = 0.0;
		for (int j = 8*i-36; j <= 8*i-29; j++)
		{
		    String p = date.shiftedBy(-10800.0*(j-1)).toString();
		    String k = p.substring(0, 4) + p.substring(5, 7) + p.substring(8, 10);
		    apVals[i] += swData.get(k)[Integer.parseInt(p.substring(11, 13))/3+14];
		}
		apVals[i] = apVals[i]/8;
	    }
	}

	return(apVals);
    }
}
