/*
 * EventHandling.java - Functions to handle orbit events.
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

import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.frames.Transform;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.events.EventDetector;
import org.orekit.propagation.events.handlers.EventHandler;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

public class EventHandling<T extends EventDetector> implements EventHandler<T>
{
    protected Settings.JSONManeuver maneuver;

    public EventHandling(Settings.JSONManeuver man)
    {
	maneuver = man;
    }

    @Override public EventHandler.Action eventOccurred(SpacecraftState state, T det, boolean incr)
    {
	return(EventHandler.Action.RESET_STATE);
    }

    @Override public SpacecraftState resetState(T det, SpacecraftState old)
    {
	KeplerianOrbit kep = new KeplerianOrbit(old.getOrbit());
	double a = kep.getA();
	double e = kep.getE();
	double i = kep.getI();
	double O = kep.getRightAscensionOfAscendingNode();
	double w = kep.getPerigeeArgument();
	double theta = kep.getTrueAnomaly();

	Orbit neworb = null;
	if (maneuver.ManeuverType.equals("NorthSouthStationing") || maneuver.ManeuverType.equals("EastWestStationing"))
	{
	    KeplerianOrbit kep1, kep2;
	    double oldmp, mp = 0.0, val1, val2, lb = 0.0, ub;
	    OneAxisEllipsoid earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
							  Constants.WGS84_EARTH_FLATTENING, DataManager.itrf);
	    if (maneuver.ManeuverType.equals("NorthSouthStationing"))
		ub = FastMath.PI;
	    else
		ub = 2.0*FastMath.PI;

	    for (int ii = 0; ii < 50; ii++)
	    {
		oldmp = mp;
		mp = 0.5*(lb + ub);
		if (maneuver.ManeuverType.equals("NorthSouthStationing"))
		{
		    kep1 = new KeplerianOrbit(a, e, lb, w, O, theta, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
		    kep2 = new KeplerianOrbit(a, e, mp, w, O, theta, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
		}
		else
		{
		    kep1 = new KeplerianOrbit(a, e, i, w, O, lb, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
		    kep2 = new KeplerianOrbit(a, e, i, w, O, mp, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
		}

		GeodeticPoint geo1 = earth.transform(kep1.getPVCoordinates().getPosition(), old.getFrame(), old.getDate());
		GeodeticPoint geo2 = earth.transform(kep2.getPVCoordinates().getPosition(), old.getFrame(), old.getDate());
		if (maneuver.ManeuverType.equals("NorthSouthStationing"))
		{
		    val1 = MathUtils.normalizeAngle(geo1.getLatitude(), 0.0) - MathUtils.normalizeAngle(maneuver.Params[0], 0.0);
		    val2 = MathUtils.normalizeAngle(geo2.getLatitude(), 0.0) - MathUtils.normalizeAngle(maneuver.Params[0], 0.0);
		}
		else
		{
		    val1 = MathUtils.normalizeAngle(geo1.getLongitude(), 0.0) - MathUtils.normalizeAngle(maneuver.Params[0], 0.0);
		    val2 = MathUtils.normalizeAngle(geo2.getLongitude(), 0.0) - MathUtils.normalizeAngle(maneuver.Params[0], 0.0);
		}

		if (FastMath.abs(oldmp - mp) < 1.0E-6)
		{
		    neworb = kep2;
		    break;
		}
		if (val1*val2 < 0.0)
		    ub = mp;
		else
		    lb = mp;
	    }
	}
	else
	{
	    if (maneuver.ManeuverType.equals("SemiMajorAxisChange"))
		a = maneuver.Params[0];
	    if (maneuver.ManeuverType.equals("PerigeeChange"))
		a = maneuver.Params[0]/(1 - e);
	    if (maneuver.ManeuverType.equals("EccentricityChange"))
		e = maneuver.Params[0];
	    if (maneuver.ManeuverType.equals("InclinationChange"))
		i = maneuver.Params[0];
	    if (maneuver.ManeuverType.equals("RAANChange"))
		O = maneuver.Params[0];
	    if (maneuver.ManeuverType.equals("ArgPerigeeChange"))
		w = maneuver.Params[0];
	    neworb = new KeplerianOrbit(a, e, i, w, O, theta, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
	}

	if (neworb == null)
	    return(old);
	return(new SpacecraftState(neworb, old.getAttitude(), old.getMass(), old.getAdditionalStates()));
    }
}
