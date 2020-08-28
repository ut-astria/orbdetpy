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
    public final static String GEO_ZONE_NAME = "Geographic Region";

    public final Settings.ManeuverType maneuverType;
    public final double delta;
    public final String stationName;
    public Boolean isVisible;

    public EventHandling(Settings.ManeuverType maneuverType, double delta, String stationName, Boolean isVisible)
    {
	this.maneuverType = maneuverType;
	this.delta = delta;
	this.stationName = stationName;
	this.isVisible = isVisible;
    }

    @Override public Action eventOccurred(SpacecraftState state, T detector, boolean increasing)
    {
	if (maneuverType == Settings.ManeuverType.UNDEFINED && stationName != null)
	{
	    String name = detector.getClass().getSimpleName();
	    if (name.equals("ElevationDetector"))
	    {
		isVisible = increasing;
		if (increasing)
		    return(Action.STOP);
	    }
	    else if (name.equals("GroundFieldOfViewDetector") || name.equals("GeographicZoneDetector"))
	    {
		isVisible = !increasing;
		if (!increasing)
		    return(Action.STOP);
	    }
	}

	if (maneuverType != Settings.ManeuverType.UNDEFINED)
	{
	    if (maneuverType == Settings.ManeuverType.STOP_PROPAGATION)
		return(Action.STOP);
	    if (delta != 0.0)
		return(Action.RESET_STATE);
	}

	return(Action.CONTINUE);
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
	if (maneuverType == Settings.ManeuverType.NORTH_SOUTH_STATIONING || maneuverType == Settings.ManeuverType.EAST_WEST_STATIONING)
	{
	    PVCoordinates pvc = old.getOrbit().getPVCoordinates();
	    OneAxisEllipsoid earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING,
							  FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP));
	    GeodeticPoint geo = earth.transform(pvc.getPosition(), old.getFrame(), old.getDate());

	    if (maneuverType == Settings.ManeuverType.NORTH_SOUTH_STATIONING)
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
	    switch (maneuverType)
	    {
	    case SEMI_MAJOR_AXIS_CHANGE:
		a += delta;
		break;
	    case PERIGEE_CHANGE:
		a += delta/(1 - e);
		break;
	    case ECCENTRICITY_CHANGE:
		e += delta;
		break;
	    case INCLINATION_CHANGE:
		i += delta;
		break;
	    case RAAN_CHANGE:
		O += delta;
		break;
	    case ARG_PERIGEE_CHANGE:
		w += delta;
		break;
	    default:
		throw(new RuntimeException("Invalid maneuver type"));
	    }

	    neworb = new KeplerianOrbit(a, e, i, w, O, theta, PositionAngle.TRUE, old.getFrame(),
					old.getDate(), old.getMu());
	}

	return(new SpacecraftState(neworb, old.getAttitude(), old.getMass(), old.getAdditionalStates()));
    }
}
