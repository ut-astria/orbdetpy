/*
 * DataManager.java - Functions for handling data files.
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

import java.io.File;
import java.util.HashMap;
import java.util.Scanner;
import org.orekit.data.DirectoryCrawler;
import org.orekit.data.DataProvidersManager;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.ITRFVersion;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.time.UT1Scale;
import org.orekit.utils.IERSConventions;

public class DataManager
{
    public static class MSISEData
    {
	public AbsoluteDate mindate;
	public AbsoluteDate maxdate;
	public HashMap<String, double[]> data;
    }

    public static String datapath;

    public static Frame gcrf;
    public static Frame itrf;
    public static Frame eme2000;
    public static TimeScale utcscale;
    public static UT1Scale ut1scale;

    public static MSISEData msisedata;

    public static void initialize(String path) throws Exception
    {
	DataManager.datapath = path;
	DataProvidersManager.getInstance().addProvider(
	    new DirectoryCrawler(new File(path)));

	gcrf = FramesFactory.getGCRF();
	itrf = FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, false);
	eme2000 = FramesFactory.getEME2000();
	utcscale = TimeScalesFactory.getUTC();
	ut1scale = TimeScalesFactory.getUT1(IERSConventions.IERS_2010, false);

	loadMSISEData();
    }

    private static void loadMSISEData() throws Exception
    {
	msisedata = new MSISEData();
	msisedata.data = new HashMap<String, double[]>();

	String[] toks = null;
	Scanner scan = new Scanner(new File(datapath, "SpaceWeather.dat"));
	while (scan.hasNextLine())
	{
	    toks = scan.nextLine().split(",");
	    for (int i = 0; i < toks.length; i++)
		toks[i] = toks[i].trim();

	    if (msisedata.mindate == null)
		msisedata.mindate = new AbsoluteDate(Integer.parseInt(toks[0]),
						     Integer.parseInt(toks[1]),
						     Integer.parseInt(toks[2]),
						     utcscale);

	    double[] vals = new double[toks.length];
	    for (int i = 0; i < toks.length; i++)
	    {
		if (toks[i].length() > 0)
		    vals[i] = Double.parseDouble(toks[i]);
		else
		    vals[i] = 0.0;
	    }

	    msisedata.data.put(String.format("%s%s%s", toks[0], toks[1],
					     toks[2]), vals);
	}

	scan.close();
	if (toks != null)
	    msisedata.maxdate = new AbsoluteDate(Integer.parseInt(toks[0]),
						 Integer.parseInt(toks[1]),
						 Integer.parseInt(toks[2]),
						 utcscale);
    }
}
