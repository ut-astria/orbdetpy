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

    public static void initialize(String dataPath) throws Exception
    {
	DataManager.dataPath = dataPath;
	DataManager.threadPool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
	DataProvidersManager.getInstance().addProvider(new DirectoryCrawler(new File(dataPath)));
	loadMSISEData();
    }

    public static void shutdown()
    {
	if (threadPool != null)
	    threadPool.shutdown();
    }

    private static void loadMSISEData() throws Exception
    {
	String[] toks = null;
	DataManager.msiseData = new MSISEData();
	DataManager.msiseData.data = new HashMap<String, double[]>();
	Scanner scan = new Scanner(new File(dataPath, "SpaceWeather.dat"));
	while (scan.hasNextLine())
	{
	    toks = scan.nextLine().split(",");
	    double[] vals = new double[toks.length];
	    for (int i = 0; i < toks.length; i++)
		toks[i] = toks[i].trim();

	    if (DataManager.msiseData.minDate == null)
		DataManager.msiseData.minDate = new AbsoluteDate(Integer.parseInt(toks[0]), Integer.parseInt(toks[1]),
								 Integer.parseInt(toks[2]), TimeScalesFactory.getUTC());
	    for (int i = 0; i < toks.length; i++)
	    {
		if (toks[i].length() > 0)
		    vals[i] = Double.parseDouble(toks[i]);
	    }
	    DataManager.msiseData.data.put(String.format("%s%s%s", toks[0], toks[1], toks[2]), vals);
	}

	scan.close();
	if (toks != null)
	    DataManager.msiseData.maxDate = new AbsoluteDate(Integer.parseInt(toks[0]), Integer.parseInt(toks[1]),
							     Integer.parseInt(toks[2]), TimeScalesFactory.getUTC());
    }
}
