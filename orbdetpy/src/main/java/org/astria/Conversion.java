/*
 * Conversion.java - Various utility functions.
 * Copyright (C) 2019 University of Texas
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

import java.util.List;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.frames.Frame;
import org.orekit.frames.TopocentricFrame;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

public final class Conversion
{
    public static double[] transformFrame(String srcFrame, String time, List<Double> pva, String destFrame)
    {
	Frame fromFrame = DataManager.getFrame(srcFrame);
	Frame toFrame = DataManager.getFrame(destFrame);
	Transform xfm = fromFrame.getTransformTo(toFrame, DataManager.parseDateTime(time));

	PVCoordinates topv = null;
	if (pva.size() == 6)
	    topv = xfm.transformPVCoordinates(new PVCoordinates(new Vector3D(pva.get(0), pva.get(1), pva.get(2)),
								new Vector3D(pva.get(3), pva.get(4), pva.get(5))));
	else
	    topv = xfm.transformPVCoordinates(new PVCoordinates(new Vector3D(pva.get(0), pva.get(1), pva.get(2)),
								new Vector3D(pva.get(3), pva.get(4), pva.get(5)),
								new Vector3D(pva.get(6), pva.get(7), pva.get(8))));

	Vector3D p = topv.getPosition();
	Vector3D v = topv.getVelocity();
	if (pva.size() == 6)
	    return(new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()});

	Vector3D a = topv.getAcceleration();
	return(new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ(), a.getX(), a.getY(), a.getZ()});
    }

    public static double[] convertAzElToRaDec(String time, double az, double el, double lat, double lon, double alt, String frame)
    {
	OneAxisEllipsoid oae = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
						    Constants.WGS84_EARTH_FLATTENING, DataManager.getFrame("ITRF"));
	TopocentricFrame fromframe = new TopocentricFrame(oae, new GeodeticPoint(lat, lon, alt), "gs");
	Transform xfm = fromframe.getTransformTo(DataManager.getFrame(frame), DataManager.parseDateTime(time));

	Vector3D tovec = xfm.transformVector(new Vector3D(FastMath.cos(el)*FastMath.sin(az),
							  FastMath.cos(el)*FastMath.cos(az), FastMath.sin(el)));
	double x = tovec.getX();
	double y = tovec.getY();
	double[] radec = new double[]{FastMath.atan2(y, x), FastMath.atan2(tovec.getZ(), FastMath.sqrt(x*x + y*y))};
	return(radec);
    }

    public static double[] convertRaDecToAzEl(String frame, String time, double ra, double dec, double lat, double lon, double alt)
    {
	OneAxisEllipsoid oae = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
						    Constants.WGS84_EARTH_FLATTENING, DataManager.getFrame("ITRF"));
	TopocentricFrame toframe = new TopocentricFrame(oae, new GeodeticPoint(lat, lon, alt), "gs");
	Transform xfm = DataManager.getFrame(frame).getTransformTo(toframe, DataManager.parseDateTime(time));

	Vector3D tovec = xfm.transformVector(new Vector3D(FastMath.cos(dec)*FastMath.cos(ra),
							  FastMath.cos(dec)*FastMath.sin(ra), FastMath.sin(dec)));
	double x = tovec.getX();
	double y = tovec.getY();
	double[] azel = new double[]{FastMath.atan2(x, y), FastMath.atan2(tovec.getZ(), FastMath.sqrt(x*x + y*y))};
	return(azel);
    }
}
