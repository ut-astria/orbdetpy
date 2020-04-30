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
import org.orekit.files.ccsds.TDMParser.TDMFileFormat;
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
    public static class InterpolationOutput
    {
	public String time;
	public double[] state;

	public InterpolationOutput(AbsoluteDate time, TimeStampedPVCoordinates pv)
	{
	    Vector3D p = pv.getPosition();
	    Vector3D v = pv.getVelocity();
	    this.time = DataManager.getUTCString(time);
	    this.state = new double[] {p.getX(), p.getY(), p.getZ(), v.getX(), v.getY(), v.getZ()};
	}
    }

    public static KeplerianOrbit iodGooding(double[] gslat, double[] gslon, double[] gsalt, Predefined frame, String[] tmstr,
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
	    time[i] = DataManager.parseDateTime(tmstr[i]);

	    GroundStation sta = new GroundStation(
		new TopocentricFrame(oae, new GeodeticPoint(gslat[i], gslon[i], gsalt[i]), Integer.toString(i)));
	    gspos[i] = sta.getBaseFrame().getPVCoordinates(time[i], FramesFactory.getFrame(frame)).getPosition();
	}

	IodGooding good = new IodGooding(FramesFactory.getFrame(frame), Constants.EGM96_EARTH_MU);
	KeplerianOrbit orb = good.estimate(gspos[0], gspos[1], gspos[2], los[0], time[0], los[1],
					   time[1], los[2], time[2], rho1init, rho3init);
	return(orb);
    }

    public static ArrayList<ArrayList<Measurements.SimulatedMeasurement>> importTDM(String fileName, String fileFormat)
    {
	Measurements.SimulatedMeasurement obj = null;
	ArrayList<ArrayList<Measurements.SimulatedMeasurement>> output =
	    new ArrayList<ArrayList<Measurements.SimulatedMeasurement>>();

	TDMFile tdm = new TDMParser().withFileFormat(TDMFileFormat.valueOf(fileFormat)).parse(fileName);
	for (TDMFile.ObservationsBlock blk : tdm.getObservationsBlocks())
	{
	    int i = 0;
	    String atype = blk.getMetaData().getAngleType();
	    ArrayList<Measurements.SimulatedMeasurement> mall = new ArrayList<Measurements.SimulatedMeasurement>();
	    for (TDMFile.Observation obs : blk.getObservations())
	    {
		String keyw = obs.getKeyword();
		if (!(keyw.equalsIgnoreCase("RANGE") || keyw.equalsIgnoreCase("DOPPLER_INSTANTANEOUS") ||
		      keyw.equalsIgnoreCase("ANGLE_1") || keyw.equalsIgnoreCase("ANGLE_2")))
		    continue;
		if (i == 0)
		    obj = new Measurements.SimulatedMeasurement();

		if (atype == null)
		{
		    if (keyw.equalsIgnoreCase("RANGE"))
			obj.range = obs.getMeasurement()*1000.0;
		    else if (keyw.equalsIgnoreCase("DOPPLER_INSTANTANEOUS"))
			obj.rangeRate = obs.getMeasurement()*1000.0;
		}
		else if (atype.equalsIgnoreCase("RADEC"))
		{
		    if (keyw.equalsIgnoreCase("ANGLE_1"))
			obj.rightAscension = obs.getMeasurement()*FastMath.PI/180.0;
		    else if (keyw.equalsIgnoreCase("ANGLE_2"))
			obj.declination = obs.getMeasurement()*FastMath.PI/180.0;
		}
		else if (atype.equalsIgnoreCase("AZEL"))
		{
		    if (keyw.equalsIgnoreCase("ANGLE_1"))
			obj.azimuth = obs.getMeasurement()*FastMath.PI/180.0;
		    else if (keyw.equalsIgnoreCase("ANGLE_2"))
			obj.elevation = obs.getMeasurement()*FastMath.PI/180.0;
		}

		if (++i == 2)
		{
		    i = 0;
		    obj.time = DataManager.getUTCString(obs.getEpoch());
		    mall.add(obj);
		}
	    }
	    output.add(mall);
	}

	return(output);
    }

    public static ArrayList<InterpolationOutput> interpolateEphemeris(Predefined sourceFrame, List<String> time, List<Double[]> ephem,
								      int numPoints, Predefined destFrame, String interpStart,
								      String interpEnd, double stepSize)
    {
	Frame fromFrame = FramesFactory.getFrame(sourceFrame);
	ArrayList<SpacecraftState> states = new ArrayList<SpacecraftState>(time.size());
	for (int i = 0; i < time.size(); i++)
	{
	    Double[] pv = ephem.get(i);
	    states.add(new SpacecraftState(new CartesianOrbit(new TimeStampedPVCoordinates(
								  DataManager.parseDateTime(time.get(i)), new Vector3D(pv[0], pv[1], pv[2]),
								  new Vector3D(pv[3], pv[4], pv[5])), fromFrame, Constants.EGM96_EARTH_MU)));
	}

	Frame toFrame = FramesFactory.getFrame(destFrame);
	ArrayList<InterpolationOutput> output = new ArrayList<InterpolationOutput>();
	Ephemeris interpolator = new Ephemeris(states, numPoints);
	AbsoluteDate tm = DataManager.parseDateTime(interpStart);
	AbsoluteDate tmFinal = DataManager.parseDateTime(interpEnd);

	while (true)
	{
	    output.add(new InterpolationOutput(tm, interpolator.getPVCoordinates(tm, toFrame)));
	    double deltat = tmFinal.durationFrom(tm);
	    if (deltat <= 0.0)
		break;
	    tm = new AbsoluteDate(tm, FastMath.min(deltat, stepSize));
	}

	return(output);
    }
}
