/*
 * WAM.java - Implementation of NOAA's WAM-IPE atmospheric model.
 * Copyright (C) 2020 University of Texas
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

import java.util.Arrays;
import java.util.Map;
import org.hipparchus.RealFieldElement;
import org.hipparchus.analysis.interpolation.TricubicInterpolator;
import org.hipparchus.geometry.euclidean.threed.FieldVector3D;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.bodies.BodyShape;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.frames.Frame;
import org.orekit.models.earth.atmosphere.Atmosphere;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.FieldAbsoluteDate;
import ucar.nc2.NetcdfFile;
import ucar.nc2.NetcdfFiles;

public final class WAM implements Atmosphere
{
    private BodyShape earth;
    private String currentFile = "";
    private float[] wamLat, wamLon, altSlice;
    private float[][][] wamAlt, wamRho;
    private TricubicInterpolator interpolator = new TricubicInterpolator();

    public WAM(BodyShape earth)
    {
	this.earth = earth;
    }

    @Override public synchronized double getDensity(AbsoluteDate date, Vector3D position, Frame frame)
    {
	double tt = date.durationFrom(AbsoluteDate.J2000_EPOCH);
	Map.Entry<Double, String> entry = DataManager.wamFileMap.floorEntry(tt);
	if (entry == null || (tt - entry.getKey() > 86400.0))
	    throw(new RuntimeException("WAM data unavailable at " + date.toString()));

	if (!currentFile.equals(entry.getValue()))
	{
	    currentFile = entry.getValue();
	    try (NetcdfFile wamData = NetcdfFiles.open(currentFile))
	    {
		wamLat = (float[])wamData.findVariable("lat").read("1:89").copyTo1DJavaArray();
		wamLon = (float[])wamData.findVariable("lon").read().copyTo1DJavaArray();
		wamAlt = (float[][][])wamData.findVariable("height").read(":,1:89,:").copyToNDJavaArray();
		wamRho = (float[][][])wamData.findVariable("thermosphere_mass_density").read(":,1:89,:").copyToNDJavaArray();
		altSlice = new float[wamAlt.length];
	    }
	    catch (Exception exc)
	    {
		throw(new RuntimeException(exc));
	    }
	}

	GeodeticPoint gp = earth.transform(position, frame, date);
	float lat = (float)FastMath.toDegrees(gp.getLatitude());
	float lon = (float)FastMath.toDegrees(MathUtils.normalizeAngle(gp.getLongitude(), FastMath.PI));
	float alt = (float)gp.getAltitude();

	int[] xb = getBounds(wamLat, lat, true);
	int[] yb = getBounds(wamLon, lon, true);
	for (int i = 0; i < altSlice.length; i++)
	    altSlice[i] = wamAlt[i][xb[0]][yb[0]];
	int[] zb = getBounds(altSlice, alt, false);

	double[][][] gridF = new double[2][2][2]; 
	double[] gridX = {wamLat[xb[0]], wamLat[xb[1]]};
	double[] gridY = {wamLon[yb[0]], wamLon[yb[1]]};
	if (yb[1] == 0)
	    gridY[1] += 360.0;
	double[] gridZ = {wamAlt[zb[0]][xb[0]][yb[0]], wamAlt[zb[1]][xb[0]][yb[0]]};
	for (int i = 0; i < 2; i++)
	    for (int j = 0; j < 2; j++)
		for (int k = 0; k < 2; k++)
		    gridF[i][j][k] = wamRho[zb[k]][xb[i]][yb[j]];
	return(interpolator.interpolate(gridX, gridY, gridZ, gridF).value(lat, lon, alt));
    }

    @Override public synchronized <T extends RealFieldElement<T>> T getDensity(FieldAbsoluteDate<T> date, FieldVector3D<T> position, Frame frame)
    {
	throw(new UnsupportedOperationException("Method is not implemented. Call double getDensity(...)."));
    }

    @Override public Frame getFrame()
    {
        return(earth.getBodyFrame());
    }

    private int[] getBounds(float[] array, float key, boolean cyclic)
    {
	int[] bounds = {0, Arrays.binarySearch(array, key)};
	if (bounds[1] < 0)
	    bounds[1] = -1 - bounds[1];
	if (bounds[1] == 0)
	{
	    if (cyclic)
		bounds[0] = array.length - 1;
	    else
		bounds[1] = 1;
	}
	else if (bounds[1] == array.length)
	{
	    if (cyclic)
	    {
		bounds[0] = array.length - 1;
		bounds[1] = 0;
	    }
	    else
	    {
		bounds[0] = array.length - 2;
		bounds[1] = array.length - 1;
	    }
	}
	else
	    bounds[0] = bounds[1] - 1;
	return(bounds);
    }
}
