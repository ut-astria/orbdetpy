/*
 * Conversion.java - Various conversion functions.
 * Copyright (C) 2019-2020 University of Texas
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
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.frames.TopocentricFrame;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateComponents;
import org.orekit.time.DateTimeComponents;
import org.orekit.time.TimeComponents;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

public final class Conversion
{
    public static double[] transformFrame(Predefined srcFrame, AbsoluteDate time, List<Double> pva, Predefined destFrame)
    {
	Frame fromFrame = FramesFactory.getFrame(srcFrame);
	Frame toFrame = FramesFactory.getFrame(destFrame);
	Transform xfm = fromFrame.getTransformTo(toFrame, time);

	if (pva.size() == 3)
	{
	    Vector3D p = xfm.transformPosition(new Vector3D(pva.get(0), pva.get(1), pva.get(2)));
	    return(new double[]{p.getX(), p.getY(), p.getZ()});
	}
	else if (pva.size() == 6)
	{
	    PVCoordinates toPv = xfm.transformPVCoordinates(new PVCoordinates(new Vector3D(pva.get(0), pva.get(1), pva.get(2)),
									      new Vector3D(pva.get(3), pva.get(4), pva.get(5))));
	    Vector3D p = toPv.getPosition();
	    Vector3D v = toPv.getVelocity();
	    return(new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()});
	}
	else
	{
	    PVCoordinates toPv = xfm.transformPVCoordinates(new PVCoordinates(new Vector3D(pva.get(0), pva.get(1), pva.get(2)),
									      new Vector3D(pva.get(3), pva.get(4), pva.get(5)),
									      new Vector3D(pva.get(6), pva.get(7), pva.get(8))));
	    Vector3D p = toPv.getPosition();
	    Vector3D v = toPv.getVelocity();
	    Vector3D a = toPv.getAcceleration();
	    return(new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ(), a.getX(), a.getY(), a.getZ()});
	}
    }

    public static double[] convertAzElToRaDec(AbsoluteDate time, double az, double el, double lat,
					      double lon, double alt, Predefined frame)
    {
	OneAxisEllipsoid oae = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING,
						    FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP));
	TopocentricFrame fromFrame = new TopocentricFrame(oae, new GeodeticPoint(lat, lon, alt), "gs");
	Transform xfm = fromFrame.getTransformTo(FramesFactory.getFrame(frame), time);

	Vector3D toVec = xfm.transformVector(new Vector3D(FastMath.cos(el)*FastMath.sin(az),
							  FastMath.cos(el)*FastMath.cos(az), FastMath.sin(el)));
	double x = toVec.getX();
	double y = toVec.getY();
	return(new double[]{FastMath.atan2(y, x), FastMath.atan2(toVec.getZ(), FastMath.sqrt(x*x + y*y))});
    }

    public static double[] convertRaDecToAzEl(Predefined frame, AbsoluteDate time, double ra, double dec,
					      double lat, double lon, double alt)
    {
	OneAxisEllipsoid oae = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING,
						    FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP));
	TopocentricFrame toframe = new TopocentricFrame(oae, new GeodeticPoint(lat, lon, alt), "gs");
	Transform xfm = FramesFactory.getFrame(frame).getTransformTo(toframe, time);

	Vector3D toVec = xfm.transformVector(new Vector3D(FastMath.cos(dec)*FastMath.cos(ra),
							  FastMath.cos(dec)*FastMath.sin(ra), FastMath.sin(dec)));
	double x = toVec.getX();
	double y = toVec.getY();
	return(new double[]{FastMath.atan2(x, y), FastMath.atan2(toVec.getZ(), FastMath.sqrt(x*x + y*y))});
    }

    public static String getUTCString(double j2000Offset)
    {
	DateTimeComponents dtc = AbsoluteDate.J2000_EPOCH.shiftedBy(j2000Offset).getComponents(TimeScalesFactory.getUTC());
	DateComponents dc = dtc.getDate();
	TimeComponents tc = dtc.getTime();

	StringBuilder sbSec = new StringBuilder(9);
	if (tc.getSecond() < 10.0)
	    sbSec.append("0");
	sbSec.append(tc.getSecond());
	if (sbSec.length() > 9)
	    sbSec.setLength(9);
	else
	{
	    while (sbSec.length() < 9)
		sbSec.append("0");
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
	sb.append(tc.getMinute()).append(":").append(sbSec).append("Z");
	return(sb.toString());
    }

    public static double getJ2000EpochOffset(String utcTime)
    {
	return(new AbsoluteDate(DateTimeComponents.parseDateTime(utcTime), TimeScalesFactory.getUTC())
	       .durationFrom(AbsoluteDate.J2000_EPOCH));
    }
}
