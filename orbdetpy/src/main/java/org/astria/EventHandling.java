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

import org.hipparchus.geometry.euclidean.threed.Rotation;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.ode.events.Action;
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

public final class EventHandling<T extends EventDetector> implements EventHandler<T>
{
    private final String trigEvent;
    private final String mnvrType;
    private final double target;
    private int steps;

    public EventHandling(String trigEvent, String mnvrType, double target, int steps)
    {
	this.trigEvent = trigEvent;
	this.mnvrType = mnvrType;
	this.target = target;
	this.steps = steps;
    }

    @Override public Action eventOccurred(SpacecraftState state, T det, boolean incr)
    {
	if (trigEvent.equals("LongitudeCrossing") && mnvrType.equals("StopPropagation"))
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
	if (mnvrType.equals("NorthSouthStationing") || mnvrType.equals("EastWestStationing"))
	{
	    PVCoordinates pvc = old.getOrbit().getPVCoordinates();
	    OneAxisEllipsoid earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
							  Constants.WGS84_EARTH_FLATTENING, DataManager.getFrame("ITRF"));
	    GeodeticPoint geo = earth.transform(pvc.getPosition(), old.getFrame(), old.getDate());

	    if (mnvrType.equals("NorthSouthStationing"))
		geo = new GeodeticPoint(geo.getLatitude() + (target - geo.getLatitude())/steps,
					geo.getLongitude(), geo.getAltitude());
	    else
		geo = new GeodeticPoint(geo.getLatitude(), geo.getLongitude() +
					(target - geo.getLongitude())/steps, geo.getAltitude());

	    Transform xfm = earth.getFrame().getTransformTo(old.getFrame(), old.getDate());
	    Vector3D newpos = xfm.transformPosition(earth.transform(geo));
	    Rotation rot = new Rotation(pvc.getPosition(), newpos);
	    neworb = new CartesianOrbit(new PVCoordinates(newpos, rot.applyTo(pvc.getVelocity())), old.getFrame(),
					old.getDate(), old.getMu());
	}
	else
	{
	    if (mnvrType.equals("SemiMajorAxisChange"))
		a += (target - a)/steps;
	    if (mnvrType.equals("PerigeeChange"))
		a += (target - a/(1 - e))/steps;
	    if (mnvrType.equals("EccentricityChange"))
		e += (target - e)/steps;
	    if (mnvrType.equals("InclinationChange"))
		i += (target - i)/steps;
	    if (mnvrType.equals("RAANChange"))
		O += (target - O)/steps;
	    if (mnvrType.equals("ArgPerigeeChange"))
		w += (target - w)/steps;
	    neworb = new KeplerianOrbit(a, e, i, w, O, theta, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
	}

	steps--;
	if (neworb == null)
	    return(old);
	return(new SpacecraftState(neworb, old.getAttitude(), old.getMass(), old.getAdditionalStates()));
    }
}
