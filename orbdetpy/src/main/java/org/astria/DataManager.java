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
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.HashMap;
import java.util.Scanner;
import java.util.TreeMap;
import java.util.concurrent.Executors;
import java.util.concurrent.ExecutorService;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataContext;
import org.orekit.data.DirectoryCrawler;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;

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
    protected static TreeMap<Double, String> wamFileMap;
    protected static OneAxisEllipsoid earthShape;
    protected static HashMap<Integer, AbsoluteDate> epochs;

    private DataManager()
    {
    }

    public static synchronized void initialize(String dataPath) throws Exception
    {
	DataManager.dataPath = dataPath;
	DataManager.threadPool = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());
	DataContext.getDefault().getDataProvidersManager().addProvider(new DirectoryCrawler(new File(dataPath)));
	DataManager.earthShape = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING,
						      FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP));
	loadEpochs();
	loadMSISEData();
	loadWAMData();
    }

    public static synchronized void shutdown()
    {
	if (threadPool != null)
	{
	    threadPool.shutdown();
	    threadPool = null;
	}
    }

    private static void loadEpochs()
    {
	epochs = new HashMap<Integer, AbsoluteDate>();
	epochs.put(0, AbsoluteDate.BEIDOU_EPOCH);
	epochs.put(1, AbsoluteDate.CCSDS_EPOCH);
	epochs.put(2, AbsoluteDate.FIFTIES_EPOCH);
	epochs.put(3, AbsoluteDate.GALILEO_EPOCH);
	epochs.put(4, AbsoluteDate.GLONASS_EPOCH);
	epochs.put(5, AbsoluteDate.GPS_EPOCH);
	epochs.put(6, AbsoluteDate.IRNSS_EPOCH);
	epochs.put(7, AbsoluteDate.J2000_EPOCH);
	epochs.put(8, AbsoluteDate.JAVA_EPOCH);
	epochs.put(9, AbsoluteDate.JULIAN_EPOCH);
	epochs.put(10, AbsoluteDate.MODIFIED_JULIAN_EPOCH);
	epochs.put(11, AbsoluteDate.QZSS_EPOCH);
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

    private static void loadWAMData() throws Exception
    {
	Path wamPath = Paths.get(dataPath, "WAM");
	DataManager.wamFileMap = new TreeMap<Double, String>();
	if (!wamPath.toFile().isDirectory())
	    return;
	Files.find(wamPath, Integer.MAX_VALUE, (p, a) -> a.isRegularFile() && p.toString().endsWith(".nc"))
	    .forEach(cdf -> {
		    String s = cdf.getFileName().toString().substring(4, 19).replace("_", "");
		    String utc = String.format("%s-%s-%sT%s:%s:%sZ", s.substring(0, 4), s.substring(4, 6), s.substring(6, 8),
					       s.substring(8, 10), s.substring(10, 12), s.substring(12, 14));
		    DataManager.wamFileMap.put(new AbsoluteDate(DateTimeComponents.parseDateTime(utc), TimeScalesFactory.getUTC())
					       .durationFrom(AbsoluteDate.J2000_EPOCH), cdf.toString());
		});
    }

    public static AbsoluteDate getEpoch(int epoch)
    {
	return(epochs.get(epoch));
    }
}
