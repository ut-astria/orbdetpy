/*
 * EventHandling.java - Functions to handle orbit events.
 * Copyright (C) 2019-2022 University of Texas
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
import org.orekit.frames.Transform;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.Orbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.events.EclipseDetector;
import org.orekit.propagation.events.ElevationDetector;
import org.orekit.propagation.events.EventDetector;
import org.orekit.propagation.events.GeographicZoneDetector;
import org.orekit.propagation.events.GroundFieldOfViewDetector;
import org.orekit.propagation.events.NodeDetector;
import org.orekit.propagation.events.handlers.EventHandler;
import org.orekit.utils.PVCoordinates;

public final class EventHandling<T extends EventDetector> implements EventHandler<T>
{
    public final static String GEO_ZONE_NAME = "GeoRegion";
    public final static String PENUMBRA = "Penumbra";
    public final static String UMBRA = "Umbra";
    public final static String ASC_NODE = "Ascend";
    public final static String DES_NODE = "Descend";

    private final Settings.ManeuverType maneuverType;
    private final double delta;
    public final String observer;
    public Boolean detected;

    public EventHandling(Settings.ManeuverType maneuverType, double delta, String observer, Boolean detected)
    {
        this.maneuverType = maneuverType;
        this.delta = delta;
        this.observer = observer;
        this.detected = detected;
    }

    @Override public Action eventOccurred(SpacecraftState state, T detector, boolean increasing)
    {
        if (maneuverType == Settings.ManeuverType.UNDEFINED)
        {
            if (detector instanceof ElevationDetector)
            {
                detected = increasing;
                if (increasing)
                    return(Action.STOP);
            }
            else if (detector instanceof GroundFieldOfViewDetector || detector instanceof GeographicZoneDetector)
            {
                detected = !increasing;
                if (!increasing)
                    return(Action.STOP);
            }
            else if (detector instanceof EclipseDetector)
                detected = !increasing;
            else if (detector instanceof NodeDetector)
            {
                if ((observer.equals(ASC_NODE) && increasing) || (observer.equals(DES_NODE) && !increasing))
                {
                    detected = true;
                    return(Action.STOP);
                }
            }
        }
        else if (maneuverType == Settings.ManeuverType.STOP_PROPAGATION)
            return(Action.STOP);
        else if (delta != 0.0)
            return(Action.RESET_STATE);
        return(Action.CONTINUE);
    }

    @Override public SpacecraftState resetState(T det, SpacecraftState old)
    {
        Orbit neworb = null;
        KeplerianOrbit kep = new KeplerianOrbit(old.getOrbit());
        double a = kep.getA();
        double e = kep.getE();
        double i = kep.getI();
        double O = kep.getRightAscensionOfAscendingNode();
        double w = kep.getPerigeeArgument();
        double theta = kep.getTrueAnomaly();

        if (maneuverType == Settings.ManeuverType.NORTH_SOUTH_STATIONING || maneuverType == Settings.ManeuverType.EAST_WEST_STATIONING)
        {
            PVCoordinates pvc = old.getOrbit().getPVCoordinates();
            GeodeticPoint geo = DataManager.earthShape.transform(pvc.getPosition(), old.getFrame(), old.getDate());
            if (maneuverType == Settings.ManeuverType.NORTH_SOUTH_STATIONING)
                geo = new GeodeticPoint(geo.getLatitude() + delta, geo.getLongitude(), geo.getAltitude());
            else
                geo = new GeodeticPoint(geo.getLatitude(), geo.getLongitude() + delta, geo.getAltitude());
            Transform xfm = DataManager.earthShape.getFrame().getTransformTo(old.getFrame(), old.getDate());
            Vector3D newpos = xfm.transformPosition(DataManager.earthShape.transform(geo));
            Rotation rot = new Rotation(pvc.getPosition(), newpos);
            neworb = new CartesianOrbit(new PVCoordinates(newpos, rot.applyTo(pvc.getVelocity())), old.getFrame(), old.getDate(), old.getMu());
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
            neworb = new KeplerianOrbit(a, e, i, w, O, theta, PositionAngle.TRUE, old.getFrame(), old.getDate(), old.getMu());
        }
        return(new SpacecraftState(neworb, old.getAttitude(), old.getMass(), old.getAdditionalStates()));
    }
}
