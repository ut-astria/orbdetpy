/*
 * Utilities.java - Various utility functions.
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

import org.astria.DataManager;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.estimation.iod.IodGooding;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.frames.Frame;
import org.orekit.frames.TopocentricFrame;
import org.orekit.frames.Transform;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

public class Utilities
{
    public static double[] transformFrame(String srcframe, String time, double[] pv, String destframe)
    {
	Frame fromframe, toframe;
	if (srcframe.equals("GCRF"))
	    fromframe = DataManager.gcrf;
	else if (srcframe.equals("ITRF"))
	    fromframe = DataManager.itrf;
	else
	    fromframe = DataManager.eme2000;

	if (destframe.equals("GCRF"))
	    toframe = DataManager.gcrf;
	else if (destframe.equals("ITRF"))
	    toframe = DataManager.itrf;
	else
	    toframe = DataManager.eme2000;

	Transform xfm = fromframe.getTransformTo(toframe, new AbsoluteDate(DateTimeComponents.parseDateTime(time),
									   DataManager.utcscale));

	PVCoordinates topv = xfm.transformPVCoordinates(new PVCoordinates(new Vector3D(pv[0], pv[1], pv[2]),
									  new Vector3D(pv[3], pv[4], pv[5])));
	Vector3D p = topv.getPosition();
	Vector3D v = topv.getVelocity();

	return(new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()});
    }

    public static KeplerianOrbit iodGooding(double[] gslat, double[] gslon, double[] gsalt, String[] tmstr,
					    double[] azi, double[] ele, double rho1init, double rho3init)
    {
	Vector3D[] los = new Vector3D[3];
	Vector3D[] gspos = new Vector3D[3];
	AbsoluteDate[] time = new AbsoluteDate[3];
	OneAxisEllipsoid oae = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
						    Constants.WGS84_EARTH_FLATTENING, DataManager.itrf);

	for (int i = 0; i < 3; i++)
	{
	    los[i] = new Vector3D(Math.sin(azi[i])*Math.cos(ele[i]), Math.cos(azi[i])*Math.cos(ele[i]), Math.sin(ele[i]));
	    time[i] = new AbsoluteDate(DateTimeComponents.parseDateTime(tmstr[i]), DataManager.utcscale);

	    GroundStation sta = new GroundStation(
		new TopocentricFrame(oae, new GeodeticPoint(gslat[i], gslon[i], gsalt[i]), Integer.toString(i)));
	    sta.getPrimeMeridianOffsetDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
	    sta.getPolarOffsetXDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
	    sta.getPolarOffsetYDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
	    gspos[i] = sta.getBaseFrame().getPVCoordinates(time[i], DataManager.eme2000).getPosition();
	}

	return(new IodGooding(DataManager.eme2000, Constants.EGM96_EARTH_MU).estimate(
		   gspos[0], gspos[1], gspos[2], los[0], time[0], los[1], time[1], los[2], time[2], rho1init, rho3init));
    }
}
