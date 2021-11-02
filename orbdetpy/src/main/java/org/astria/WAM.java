/*
 * WAM.java - Implementation of NOAA's WAM-IPE atmospheric model.
 * Copyright (C) 2020-2021 University of Texas
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

import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.regex.Pattern;
import org.hipparchus.CalculusFieldElement;
import org.hipparchus.analysis.interpolation.TricubicInterpolatingFunction;
import org.hipparchus.analysis.interpolation.TricubicInterpolator;
import org.hipparchus.geometry.euclidean.threed.FieldVector3D;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.frames.Frame;
import org.orekit.models.earth.atmosphere.Atmosphere;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.time.FieldAbsoluteDate;
import org.orekit.time.TimeScalesFactory;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFiles;
import ucar.nc2.Variable;

public final class WAM implements Atmosphere
{
    private final class CacheEntry
    {
	public float[] latitude, longitude;
	public float[][][] altitude, density, temperature, densityN2, densityO, densityO2;
    }

    private final class CacheMap<String, CacheEntry> extends LinkedHashMap<String, CacheEntry>
    {
	@Override protected boolean removeEldestEntry(Map.Entry eldest)
	{
	    return(size() > 16);
	}
    }

    private static WAM singleton;
    private boolean overrideTime;
    private TreeMap<Double, String> metaData;
    private final TricubicInterpolator interpolator = new TricubicInterpolator();
    private final CacheMap<String, CacheEntry> dataCache = new CacheMap<String, CacheEntry>();
    private static final double[] ATOMIC_MASS = {28.0*1.660539E-27, 16.0*1.660539E-27, 32.0*1.660539E-27}; // N2, O, O

    private WAM() throws Exception
    {
	Path wamPath = Paths.get(DataManager.dataPath, "WAM");
	metaData = new TreeMap<Double, String>();
	if (!wamPath.toFile().isDirectory())
	    return;

	Files.find(wamPath, Integer.MAX_VALUE, (p, a) -> a.isRegularFile() && p.toString().endsWith(".nc"))
	    .forEach(cdf -> {
		    for (String s: cdf.getFileName().toString().split(Pattern.quote(".")))
		    {
			if (!s.startsWith("20") || s.length() != 15)
			    continue;
			s = s.replace("_", "");
			String utc = String.format("%s-%s-%sT%s:%s:%sZ", s.substring(0, 4), s.substring(4, 6), s.substring(6, 8),
						   s.substring(8, 10), s.substring(10, 12), s.substring(12, 14));
			metaData.put(new AbsoluteDate(DateTimeComponents.parseDateTime(utc), TimeScalesFactory.getUTC())
				     .durationFrom(AbsoluteDate.J2000_EPOCH), cdf.toString());
			break;
		    }
		});

	String envVar = System.getenv("ORBDETPY_WAM_OVERRIDE");
	overrideTime = envVar != null && (envVar.equals("1") || envVar.equalsIgnoreCase("true"));
	if (overrideTime && metaData.size() > 0)
	    metaData.put(0.0, metaData.firstEntry().getValue());
    }

    public static synchronized WAM getInstance()
    {
	if (singleton == null)
	{
	    try
	    {
		singleton = new WAM();
	    }
	    catch (Exception exc)
	    {
		throw(new RuntimeException(exc));
	    }
	}
	return(singleton);
    }

    @Override public double getDensity(AbsoluteDate date, Vector3D position, Frame frame)
    {
	double tt = date.durationFrom(AbsoluteDate.J2000_EPOCH);
	Map.Entry<Double, String> finfo = metaData.floorEntry(tt);
	if (finfo == null || (!overrideTime && tt - finfo.getKey() > 86400.0))
	    throw(new RuntimeException("WAM data not found for " + date.toString()));

	GeodeticPoint gp = DataManager.earthShape.transform(position, frame, date);
	double lat = FastMath.toDegrees(gp.getLatitude());
	double lon = FastMath.toDegrees(MathUtils.normalizeAngle(gp.getLongitude(), FastMath.PI));
	double alt = gp.getAltitude();

	CacheEntry entry;
	synchronized (dataCache)
	{
	    entry = dataCache.get(finfo.getValue());
	    if (entry == null)
	    {
		entry = new CacheEntry();
		dataCache.put(finfo.getValue(), entry);
		try (NetcdfFile data = NetcdfFiles.open(finfo.getValue()))
		{
		    entry.latitude = (float[])data.findVariable("lat").read("1:90").copyTo1DJavaArray();
		    entry.longitude = (float[])data.findVariable("lon").read().copyTo1DJavaArray();
		    entry.altitude = (float[][][])data.findVariable("height").read(":,1:90,:").copyToNDJavaArray();

		    Variable var = data.findVariable("thermosphere_mass_density");
		    if (var == null)
			var = data.findVariable("neutral_density");
		    entry.density = (float[][][])var.read(":,1:90,:").copyToNDJavaArray();

		    var = data.findVariable("temp_neutral");
		    if (var != null)
		    {
			entry.temperature = (float[][][])var.read("149,1:90,:").copyToNDJavaArray();
			entry.densityN2 = (float[][][])data.findVariable("N2_Density").read("149,1:90,:").copyToNDJavaArray();
			entry.densityO = (float[][][])data.findVariable("O_Density").read("149,1:90,:").copyToNDJavaArray();
			entry.densityO2 = (float[][][])data.findVariable("O2_Density").read("149,1:90,:").copyToNDJavaArray();
		    }
		}
		catch (Exception exc)
		{
		    throw(new RuntimeException(exc));
		}
	    }
	}

	int[] xb = angleBounds(entry.latitude, (float)lat);
	double[] gridX = {entry.latitude[xb[0]], entry.latitude[xb[1]]};

	int[] yb = angleBounds(entry.longitude, (float)lon);
	double[] gridY = {entry.longitude[yb[0]], entry.longitude[yb[1]]};
	if (gridY[1] < gridY[0])
	    gridY[1] += 360.0;

	int[] zb = altitudeBounds(entry.altitude, xb[0], yb[0], (float)alt);
	double[] gridZ = {entry.altitude[zb[0]][xb[0]][yb[0]], entry.altitude[zb[1]][xb[0]][yb[0]]};

	double[][][] gridF = new double[gridX.length][gridY.length][gridZ.length];
	for (int i = 0; i < gridX.length; i++)
	    for (int j = 0; j < gridY.length; j++)
		for (int k = 0; k < gridZ.length; k++)
		    gridF[i][j][k] = entry.density[zb[k]][xb[i]][yb[j]];

	TricubicInterpolatingFunction function = interpolator.interpolate(gridX, gridY, gridZ, gridF);
	if (alt >= gridZ[0] && alt <= gridZ[gridZ.length - 1])
	    return(function.value(lat, lon, alt));
	if (alt < gridZ[0])
	    return(function.value(lat, lon, gridZ[0]));

	double[] species = {entry.densityN2[0][xb[0]][yb[0]], entry.densityO[0][xb[0]][yb[0]], entry.densityO2[0][xb[0]][yb[0]]};
	while (gridZ[gridZ.length - 1] < alt)
	{
	    double scale = -FastMath.min(alt - gridZ[gridZ.length - 1], 10E3)*9.80665*FastMath.pow(
		6371008.8/(gridZ[gridZ.length - 1] + 6371008.8), 2)/(1.380649E-23*entry.temperature[0][xb[0]][yb[0]]);
	    for (int i = 0; i < ATOMIC_MASS.length; i++)
		species[i] *= FastMath.exp(scale*ATOMIC_MASS[i]);
	    gridZ[gridZ.length - 1] += 10E3;
	}

	species[0] *= ATOMIC_MASS[0];
	for (int i = 1; i < ATOMIC_MASS.length; i++)
	    species[0] += species[i]*ATOMIC_MASS[i];
	return(species[0]);
    }

    @Override public <T extends CalculusFieldElement<T>> T getDensity(FieldAbsoluteDate<T> date, FieldVector3D<T> position, Frame frame)
    {
	throw(new UnsupportedOperationException("Not implemented; use double getDensity(...)."));
    }

    @Override public Frame getFrame()
    {
        return(DataManager.earthShape.getBodyFrame());
    }

    private int[] angleBounds(float[] array, float value)
    {
	int[] bounds = {0, Arrays.binarySearch(array, value)};
	if (bounds[1] < 0)
	    bounds[1] = -1 - bounds[1];
	if (bounds[1] == 0)
	    bounds[0] = array.length - 1;
	else if (bounds[1] == array.length)
	{
	    bounds[0] = array.length - 1;
	    bounds[1] = 0;
	}
	else
	    bounds[0] = bounds[1] - 1;
	return(bounds);
    }

    private int[] altitudeBounds(float[][][] altitudes, int xb, int yb, float value)
    {
	float[] array = new float[altitudes.length];
	for (int i = 0; i < array.length; i++)
	    array[i] = altitudes[i][xb][yb];

	int[] bounds = {0, Arrays.binarySearch(array, value)};
	if (bounds[1] < 0)
	    bounds[1] = -1 - bounds[1];
	if (bounds[1] == 0)
	    bounds[1] = 1;
	else if (bounds[1] == array.length)
	{
	    bounds[0] = array.length - 2;
	    bounds[1] = array.length - 1;
	}
	else
	    bounds[0] = bounds[1] - 1;
	return(bounds);
    }
}
