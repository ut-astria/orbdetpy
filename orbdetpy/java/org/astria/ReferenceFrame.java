/*
 * ReferenceFrame.java - Functions for reference frame transformations.
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
import org.orekit.frames.Frame;
import org.orekit.frames.Transform;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.utils.PVCoordinates;

public class ReferenceFrame
{
    public static double[] transform(String srcframe, String time, double[] pv, String destframe)
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

	Transform xfm = fromframe.getTransformTo(toframe, new AbsoluteDate(
						     DateTimeComponents.parseDateTime(time),
						     DataManager.utcscale));

	PVCoordinates topv = xfm.transformPVCoordinates(new PVCoordinates(
							    new Vector3D(pv[0], pv[1], pv[2]),
							    new Vector3D(pv[3], pv[4], pv[5])));
	Vector3D p = topv.getPosition();
	Vector3D v = topv.getVelocity();

	return(new double[]{p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()});
    }
}
