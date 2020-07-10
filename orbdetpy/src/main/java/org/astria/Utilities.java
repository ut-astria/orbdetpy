/*
 * Utilities.java - Various utility functions.
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

import java.util.ArrayList;
import java.util.List;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.estimation.iod.IodGooding;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.files.ccsds.TDMFile;
import org.orekit.files.ccsds.TDMParser;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.frames.TopocentricFrame;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.Ephemeris;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class Utilities
{
    public static KeplerianOrbit iodGooding(double[] gslat, double[] gslon, double[] gsalt, Predefined frame, double[] tm,
					    double[] ra, double[] dec, double rho1init, double rho3init)
    {
	Vector3D[] los = new Vector3D[3];
	Vector3D[] gspos = new Vector3D[3];
	AbsoluteDate[] time = new AbsoluteDate[3];
	OneAxisEllipsoid oae = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING,
						    FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP));

	for (int i = 0; i < 3; i++)
	{
	    los[i] = new Vector3D(FastMath.cos(dec[i])*FastMath.cos(ra[i]),
				  FastMath.cos(dec[i])*FastMath.sin(ra[i]), FastMath.sin(dec[i]));
	    time[i] = AbsoluteDate.J2000_EPOCH.shiftedBy(tm[i]);

	    GroundStation sta = new GroundStation(
		new TopocentricFrame(oae, new GeodeticPoint(gslat[i], gslon[i], gsalt[i]), Integer.toString(i)));
	    gspos[i] = sta.getBaseFrame().getPVCoordinates(time[i], FramesFactory.getFrame(frame)).getPosition();
	}

	IodGooding good = new IodGooding(FramesFactory.getFrame(frame), Constants.EGM96_EARTH_MU);
	KeplerianOrbit orb = good.estimate(gspos[0], gspos[1], gspos[2], los[0], time[0], los[1],
					   time[1], los[2], time[2], rho1init, rho3init);
	return(orb);
    }

    public static ArrayList<ArrayList<Measurements.Measurement>> importTDM(String fileName, TDMParser.TDMFileFormat fileFormat)
    {
	Measurements.Measurement obj = null;
	ArrayList<ArrayList<Measurements.Measurement>> output = new ArrayList<ArrayList<Measurements.Measurement>>();

	TDMFile tdm = new TDMParser().withFileFormat(fileFormat).parse(fileName);
	for (TDMFile.ObservationsBlock blk: tdm.getObservationsBlocks())
	{
	    int i = 0;
	    String atype = blk.getMetaData().getAngleType();
	    ArrayList<Measurements.Measurement> mall = new ArrayList<Measurements.Measurement>();
	    for (TDMFile.Observation obs: blk.getObservations())
	    {
		String keyw = obs.getKeyword();
		if (!(keyw.equals("RANGE") || keyw.equals("DOPPLER_INSTANTANEOUS") || keyw.equals("ANGLE_1") || keyw.equals("ANGLE_2")))
		    continue;
		if (i == 0)
		{
		    obj = new Measurements.Measurement();
		    obj.values = new double[2];
		}

		if (atype == null)
		{
		    if (keyw.equals("RANGE"))
			obj.values[0] = obs.getMeasurement()*1000.0;
		    else if (keyw.equals("DOPPLER_INSTANTANEOUS"))
			obj.values[1] = obs.getMeasurement()*1000.0;
		}
		else if (atype.equals("RADEC"))
		{
		    if (keyw.equals("ANGLE_1"))
			obj.values[0] = obs.getMeasurement()*FastMath.PI/180.0;
		    else if (keyw.equals("ANGLE_2"))
			obj.values[1] = obs.getMeasurement()*FastMath.PI/180.0;
		}
		else if (atype.equals("AZEL"))
		{
		    if (keyw.equals("ANGLE_1"))
			obj.values[0] = obs.getMeasurement()*FastMath.PI/180.0;
		    else if (keyw.equals("ANGLE_2"))
			obj.values[1] = obs.getMeasurement()*FastMath.PI/180.0;
		}

		if (++i == 2)
		{
		    i = 0;
		    obj.time = obs.getEpoch();
		    mall.add(obj);
		}
	    }
	    output.add(mall);
	}

	return(output);
    }

    public static ArrayList<Measurements.Measurement> interpolateEphemeris(Predefined sourceFrame, List<AbsoluteDate> time, List<Double[]> ephem,
									   int numPoints, Predefined destFrame, AbsoluteDate interpStart,
									   AbsoluteDate interpEnd, double stepSize)
    {
	Frame fromFrame = FramesFactory.getFrame(sourceFrame);
	ArrayList<SpacecraftState> states = new ArrayList<SpacecraftState>(time.size());
	for (int i = 0; i < time.size(); i++)
	{
	    Double[] pv = ephem.get(i);
	    states.add(new SpacecraftState(new CartesianOrbit(new TimeStampedPVCoordinates(time.get(i), new Vector3D(pv[0], pv[1], pv[2]),
								  new Vector3D(pv[3], pv[4], pv[5])), fromFrame, Constants.EGM96_EARTH_MU)));
	}

	AbsoluteDate tm = interpStart;
	Frame toFrame = FramesFactory.getFrame(destFrame);
	ArrayList<Measurements.Measurement> output = new ArrayList<Measurements.Measurement>();
	Ephemeris interpolator = new Ephemeris(states, numPoints);

	while (true)
	{
	    output.add(new Measurements.Measurement(interpolator.getPVCoordinates(tm, toFrame)));
	    double deltat = interpEnd.durationFrom(tm);
	    if (deltat <= 0.0)
		break;
	    tm = new AbsoluteDate(tm, FastMath.min(deltat, stepSize));
	}

	return(output);
    }
}
