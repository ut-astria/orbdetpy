/*
 * Conversion.java - Various conversion functions.
 * Copyright (C) 2019-2021 University of Texas
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
import org.hipparchus.util.MathUtils;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.frames.TopocentricFrame;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.PVCoordinates;

public final class Conversion
{
    private Conversion()
    {
    }

    public static double[] transformFrame(Predefined srcFrame, AbsoluteDate time, List<Double> pva, Predefined destFrame)
    {
	Transform xfm = FramesFactory.getFrame(srcFrame).getTransformTo(FramesFactory.getFrame(destFrame), time);
	if (pva.size() == 3)
	    return(xfm.transformPosition(new Vector3D(pva.get(0), pva.get(1), pva.get(2))).toArray());
	else if (pva.size() == 6)
	{
	    PVCoordinates toPv = xfm.transformPVCoordinates(new PVCoordinates(new Vector3D(pva.get(0), pva.get(1), pva.get(2)),
									      new Vector3D(pva.get(3), pva.get(4), pva.get(5))));
	    double[] p = toPv.getPosition().toArray();
	    double[] v = toPv.getVelocity().toArray();
	    return(new double[]{p[0], p[1], p[2], v[0], v[1], v[2]});
	}
	else
	{
	    PVCoordinates toPv = xfm.transformPVCoordinates(new PVCoordinates(new Vector3D(pva.get(0), pva.get(1), pva.get(2)),
									      new Vector3D(pva.get(3), pva.get(4), pva.get(5)),
									      new Vector3D(pva.get(6), pva.get(7), pva.get(8))));
	    double[] p = toPv.getPosition().toArray();
	    double[] v = toPv.getVelocity().toArray();
	    double[] a = toPv.getAcceleration().toArray();
	    return(new double[]{p[0], p[1], p[2], v[0], v[1], v[2], a[0], a[1], a[2]});
	}
    }

    public static double[] convertAzElToRaDec(AbsoluteDate time, double az, double el, double lat,
					      double lon, double alt, Predefined frame)
    {
	TopocentricFrame fromFrame = new TopocentricFrame(DataManager.earthShape, new GeodeticPoint(lat, lon, alt), "gs");
	Transform xfm = fromFrame.getTransformTo(FramesFactory.getFrame(frame), time);
	Vector3D toVec = xfm.transformVector(new Vector3D(FastMath.cos(el)*FastMath.sin(az), FastMath.cos(el)*FastMath.cos(az), FastMath.sin(el)));
	double[] xyz = toVec.toArray();
	return(new double[]{MathUtils.normalizeAngle(FastMath.atan2(xyz[1], xyz[0]), FastMath.PI),
			    FastMath.atan2(xyz[2], FastMath.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]))});
    }

    public static double[] convertRaDecToAzEl(Predefined frame, AbsoluteDate time, double ra, double dec,
					      double lat, double lon, double alt)
    {
	TopocentricFrame toframe = new TopocentricFrame(DataManager.earthShape, new GeodeticPoint(lat, lon, alt), "gs");
	Transform xfm = FramesFactory.getFrame(frame).getTransformTo(toframe, time);
	Vector3D toVec = xfm.transformVector(new Vector3D(FastMath.cos(dec)*FastMath.cos(ra), FastMath.cos(dec)*FastMath.sin(ra), FastMath.sin(dec)));
	double[] xyz = toVec.toArray();
	return(new double[]{MathUtils.normalizeAngle(FastMath.atan2(xyz[0], xyz[1]), FastMath.PI),
			    FastMath.atan2(xyz[2], FastMath.sqrt(xyz[0]*xyz[0]+xyz[1]*xyz[1]))});
    }
}
