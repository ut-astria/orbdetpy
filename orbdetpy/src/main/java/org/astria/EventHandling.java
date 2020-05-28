/*
 * EventHandling.java - Functions to handle orbit events.
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

import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.events.Action;
import org.hipparchus.util.MathUtils;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
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

public final class EventHandling<T extends EventDetector> implements EventHandler<T>
{
    private final String mnvrType;
    private final double delta;

    public EventHandling(String mnvrType, double delta)
    {
	this.mnvrType = mnvrType;
	this.delta = delta;
    }

    @Override public Action eventOccurred(SpacecraftState state, T det, boolean incr)
    {
	if (mnvrType != null && mnvrType.equalsIgnoreCase("StopPropagation"))
	    return(Action.STOP);
	return(Action.RESET_STATE);
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
	if (mnvrType.equalsIgnoreCase("NorthSouthStationing") || mnvrType.equalsIgnoreCase("EastWestStationing"))
	{
	    PVCoordinates pvc = old.getOrbit().getPVCoordinates();
	    OneAxisEllipsoid earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING,
							  FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP));
	    GeodeticPoint geo = earth.transform(pvc.getPosition(), old.getFrame(), old.getDate());

	    if (mnvrType.equalsIgnoreCase("NorthSouthStationing"))
		geo = new GeodeticPoint(geo.getLatitude() + delta, geo.getLongitude(), geo.getAltitude());
	    else
		geo = new GeodeticPoint(geo.getLatitude(), geo.getLongitude() + delta, geo.getAltitude());

	    Transform xfm = earth.getFrame().getTransformTo(old.getFrame(), old.getDate());
	    Vector3D newpos = xfm.transformPosition(earth.transform(geo));
	    Rotation rot = new Rotation(pvc.getPosition(), newpos);
	    neworb = new CartesianOrbit(new PVCoordinates(newpos, rot.applyTo(pvc.getVelocity())), old.getFrame(),
					old.getDate(), old.getMu());
	}
	else
	{
	    if (mnvrType.equalsIgnoreCase("SemiMajorAxisChange"))
		a += delta;
	    else if (mnvrType.equalsIgnoreCase("PerigeeChange"))
		a += delta/(1 - e);
	    else if (mnvrType.equalsIgnoreCase("EccentricityChange"))
		e += delta;
	    else if (mnvrType.equalsIgnoreCase("InclinationChange"))
		i += delta;
	    else if (mnvrType.equalsIgnoreCase("RAANChange"))
		O += delta;
	    else if (mnvrType.equalsIgnoreCase("ArgPerigeeChange"))
		w += delta;
	    else
		throw(new RuntimeException("Invalid maneuver type"));

	    neworb = new KeplerianOrbit(a, e, i, w, O, theta, PositionAngle.TRUE, old.getFrame(),
					old.getDate(), old.getMu());
	}

	return(new SpacecraftState(neworb, old.getAttitude(), old.getMass(), old.getAdditionalStates()));
    }
}
