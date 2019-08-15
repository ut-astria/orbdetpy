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
	if (maneuver.TriggerEvent.equals("LongitudeCrossing") && maneuver.ManeuverType.equals("StopPropagation"))
	    return(EventHandler.Action.STOP);
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
	    double lb, ub, target = MathUtils.normalizeAngle(maneuver.ManeuverParams[0], 0.0);
	    OneAxisEllipsoid earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
							  Constants.WGS84_EARTH_FLATTENING, DataManager.itrf);
	    GeodeticPoint geo = earth.transform(kep.getPVCoordinates().getPosition(), old.getFrame(), old.getDate());

	    if (maneuver.ManeuverType.equals("NorthSouthStationing"))
	    {
		lb = 0.0;
		ub = FastMath.PI;
		if (maneuver.ManeuverParams.length > 1)
		    target = MathUtils.normalizeAngle(geo.getLatitude() + maneuver.ManeuverParams[0], 0.0);
	    }
	    else
	    {
		lb = theta;
		ub = lb + 2.0*FastMath.PI;
		if (maneuver.ManeuverParams.length > 1)
		    target = MathUtils.normalizeAngle(geo.getLongitude() + maneuver.ManeuverParams[0], 0.0);
	    }

	    double del, mindel = 1.0E9, minval = lb;
	    for (double val = lb; val <= ub; val += (ub - lb)/3600.0)
	    {
		if (maneuver.ManeuverType.equals("NorthSouthStationing"))
		{
		    kep = new KeplerianOrbit(a, e, val, w, O, theta, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
		    geo = earth.transform(kep.getPVCoordinates().getPosition(), old.getFrame(), old.getDate());
		    del = MathUtils.normalizeAngle(geo.getLatitude(), 0.0) - target;
		}
		else
		{
		    kep = new KeplerianOrbit(a, e, i, w, O, val, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
		    geo = earth.transform(kep.getPVCoordinates().getPosition(), old.getFrame(), old.getDate());
		    del = MathUtils.normalizeAngle(geo.getLongitude(), 0.0) - target;
		}

		if (FastMath.abs(del) < FastMath.abs(mindel))
		{
		    mindel = del;
		    minval = val;
		}
	    }

	    if (maneuver.ManeuverType.equals("NorthSouthStationing"))
		neworb = new KeplerianOrbit(a, e, minval, w, O, theta, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
	    else
		neworb = new KeplerianOrbit(a, e, i, w, O, minval, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
	}
	else
	{
	    if (maneuver.ManeuverParams.length == 1)
	    {
		if (maneuver.ManeuverType.equals("SemiMajorAxisChange"))
		    a = maneuver.ManeuverParams[0];
		if (maneuver.ManeuverType.equals("PerigeeChange"))
		    a = maneuver.ManeuverParams[0]/(1 - e);
		if (maneuver.ManeuverType.equals("EccentricityChange"))
		    e = maneuver.ManeuverParams[0];
		if (maneuver.ManeuverType.equals("InclinationChange"))
		    i = maneuver.ManeuverParams[0];
		if (maneuver.ManeuverType.equals("RAANChange"))
		    O = maneuver.ManeuverParams[0];
		if (maneuver.ManeuverType.equals("ArgPerigeeChange"))
		    w = maneuver.ManeuverParams[0];
	    }
	    else
	    {
		if (maneuver.ManeuverType.equals("SemiMajorAxisChange"))
		    a += maneuver.ManeuverParams[0];
		if (maneuver.ManeuverType.equals("PerigeeChange"))
		    a += maneuver.ManeuverParams[0]/(1 - e);
		if (maneuver.ManeuverType.equals("EccentricityChange"))
		    e += maneuver.ManeuverParams[0];
		if (maneuver.ManeuverType.equals("InclinationChange"))
		    i += maneuver.ManeuverParams[0];
		if (maneuver.ManeuverType.equals("RAANChange"))
		    O += maneuver.ManeuverParams[0];
		if (maneuver.ManeuverType.equals("ArgPerigeeChange"))
		    w += maneuver.ManeuverParams[0];
	    }

	    neworb = new KeplerianOrbit(a, e, i, w, O, theta, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
	}

	if (neworb == null)
	    return(old);
	return(new SpacecraftState(neworb, old.getAttitude(), old.getMass(), old.getAdditionalStates()));
    }
}
