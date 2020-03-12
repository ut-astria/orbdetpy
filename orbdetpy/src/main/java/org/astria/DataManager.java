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
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.ITRFVersion;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateComponents;
import org.orekit.time.DateTimeComponents;
import org.orekit.time.TimeComponents;
import org.orekit.time.TimeScale;
import org.orekit.time.TimeScalesFactory;
import org.orekit.time.UT1Scale;
import org.orekit.utils.IERSConventions;

public final class DataManager
{
    public static class MSISEData
    {
	public AbsoluteDate minDate;
	public AbsoluteDate maxDate;
	public HashMap<String, double[]> data;
    }

    protected static ExecutorService threadPool;
    private static String dataPath;

    private static HashMap<String, TimeScale> timeScales;
    private static HashMap<String, Frame> refFrames;
    protected static MSISEData msiseData;

    private DataManager()
    {
    }

    public static void initialize(String path) throws Exception
    {
	threadPool = Executors.newFixedThreadPool(FastMath.min(Runtime.getRuntime().availableProcessors(), 1024));

	DataManager.dataPath = path;
	DataProvidersManager.getInstance().addProvider(new DirectoryCrawler(new File(path)));

	timeScales = new HashMap<String, TimeScale>();
	timeScales.put("TT", TimeScalesFactory.getTT());
	timeScales.put("UTC", TimeScalesFactory.getUTC());
	timeScales.put("UT1", TimeScalesFactory.getUT1(IERSConventions.IERS_2010, false));

	refFrames = new HashMap<String, Frame>();
	refFrames.put("EME2000", FramesFactory.getEME2000());
	refFrames.put("GCRF", FramesFactory.getGCRF());
	refFrames.put("ICRF", FramesFactory.getICRF());
	refFrames.put("ITRF", FramesFactory.getITRF(ITRFVersion.ITRF_2014, IERSConventions.IERS_2010, false));
	refFrames.put("MOD", FramesFactory.getMOD(IERSConventions.IERS_2010));
	refFrames.put("TEME", FramesFactory.getTEME());
	refFrames.put("TOD", FramesFactory.getTOD(IERSConventions.IERS_2010, false));

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
						     Integer.parseInt(toks[2]), getTimeScale("UTC"));

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
						 Integer.parseInt(toks[2]), getTimeScale("UTC"));
    }

    public static TimeScale getTimeScale(String name)
    {
	if (name == null || !timeScales.containsKey(name))
	    return(timeScales.get("UTC"));
	return(timeScales.get(name));
    }

    public static Frame getFrame(String name)
    {
	if (name != null)
 	    return(refFrames.get(name));
	else
	    return(refFrames.get("EME2000"));
    }

    public static String getUTCString(AbsoluteDate time)
    {
	DateTimeComponents dtc = time.getComponents(getTimeScale("UTC"));
	DateComponents dc = dtc.getDate();
	TimeComponents tc = dtc.getTime();

	StringBuilder sbsec = new StringBuilder(9);
	if (tc.getSecond() < 10.0)
	    sbsec.append("0");
	sbsec.append(tc.getSecond());
	if (sbsec.length() > 9)
	    sbsec.setLength(9);
	else
	{
	    while (sbsec.length() < 9)
		sbsec.append("0");
	}

	StringBuilder sb = new StringBuilder(27).append(dc.getYear()).append("-");
	if (dc.getMonth() < 10)
	    sb.append("0");
	sb.append(dc.getMonth()).append("-");
	if (dc.getDay() < 10)
	    sb.append("0");
	sb.append(dc.getDay()).append("T");
	if (tc.getHour() < 10)
	    sb.append("0");
	sb.append(tc.getHour()).append(":");
	if (tc.getMinute() < 10)
	    sb.append("0");
	sb.append(tc.getMinute()).append(":").append(sbsec).append("Z");
	return(sb.toString());
    }

    public static AbsoluteDate parseDateTime(String time)
    {
	return(new AbsoluteDate(DateTimeComponents.parseDateTime(time), getTimeScale("UTC")));
    }
}
