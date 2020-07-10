/*
 * DataManager.java - Functions for handling data files.
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

import java.io.File;
import java.util.HashMap;
import java.util.Scanner;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import org.hipparchus.util.FastMath;
import org.orekit.data.DirectoryCrawler;
import org.orekit.data.DataProvidersManager;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.TimeScalesFactory;

public final class DataManager
{
    public static class MSISEData
    {
	public AbsoluteDate minDate;
	public AbsoluteDate maxDate;
	public HashMap<String, double[]> data;
    }

    private static String dataPath;
    protected static ExecutorService threadPool;
    protected static MSISEData msiseData;

    private DataManager()
    {
    }

    public static void initialize(String path) throws Exception
    {
	DataManager.dataPath = path;
	DataProvidersManager.getInstance().addProvider(new DirectoryCrawler(new File(path)));
	threadPool = Executors.newFixedThreadPool(FastMath.min(Runtime.getRuntime().availableProcessors(), 1024));
	loadMSISEData();
    }

    private static void loadMSISEData() throws Exception
    {
	msiseData = new MSISEData();
	msiseData.data = new HashMap<String, double[]>();

	String[] toks = null;
	Scanner scan = new Scanner(new File(dataPath, "SpaceWeather.dat"));
	while (scan.hasNextLine())
	{
	    toks = scan.nextLine().split(",");
	    for (int i = 0; i < toks.length; i++)
		toks[i] = toks[i].trim();

	    if (msiseData.minDate == null)
		msiseData.minDate = new AbsoluteDate(Integer.parseInt(toks[0]), Integer.parseInt(toks[1]),
						     Integer.parseInt(toks[2]), TimeScalesFactory.getUTC());

	    double[] vals = new double[toks.length];
	    for (int i = 0; i < toks.length; i++)
	    {
		if (toks[i].length() > 0)
		    vals[i] = Double.parseDouble(toks[i]);
		else
		    vals[i] = 0.0;
	    }

	    msiseData.data.put(String.format("%s%s%s", toks[0], toks[1], toks[2]), vals);
	}

	scan.close();
	if (toks != null)
	    msiseData.maxDate = new AbsoluteDate(Integer.parseInt(toks[0]), Integer.parseInt(toks[1]),
						 Integer.parseInt(toks[2]), TimeScalesFactory.getUTC());
    }
}
