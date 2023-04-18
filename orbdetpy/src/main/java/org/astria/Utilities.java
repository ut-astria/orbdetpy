/*
 * Utilities.java - Various utility functions.
 * Copyright (C) 2019-2023 University of Texas
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
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.data.DataSource;
import org.orekit.estimation.iod.IodGooding;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.files.ccsds.ndm.ParserBuilder;
import org.orekit.files.ccsds.ndm.tdm.AngleType;
import org.orekit.files.ccsds.ndm.tdm.Observation;
import org.orekit.files.ccsds.ndm.tdm.ObservationsBlock;
import org.orekit.files.ccsds.ndm.tdm.ObservationType;
import org.orekit.files.ccsds.ndm.tdm.Tdm;
import org.orekit.files.ccsds.ndm.tdm.TdmMetadata;
import org.orekit.files.ccsds.ndm.tdm.TdmParser;
import org.orekit.files.ccsds.section.Segment;
import org.orekit.files.ccsds.utils.FileFormat;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.frames.TopocentricFrame;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;

public final class Utilities {
	private Utilities() {
	}

	public static KeplerianOrbit iodGooding(double[] gslat, double[] gslon, double[] gsalt, Predefined frame,
			double[] tm, double[] ra, double[] dec, double rho1init, double rho3init) {
		Vector3D[] los = new Vector3D[3];
		Vector3D[] gspos = new Vector3D[3];
		AbsoluteDate[] time = new AbsoluteDate[3];
		for (int i = 0; i < 3; i++) {
			los[i] = new Vector3D(FastMath.cos(dec[i]) * FastMath.cos(ra[i]),
					FastMath.cos(dec[i]) * FastMath.sin(ra[i]), FastMath.sin(dec[i]));
			time[i] = AbsoluteDate.J2000_EPOCH.shiftedBy(tm[i]);
			GroundStation sta = new GroundStation(new TopocentricFrame(DataManager.earthShape,
					new GeodeticPoint(gslat[i], gslon[i], gsalt[i]), Integer.toString(i)));
			gspos[i] = sta.getBaseFrame().getPVCoordinates(time[i], FramesFactory.getFrame(frame)).getPosition();
		}

		IodGooding good = new IodGooding(Constants.EGM96_EARTH_MU);
		KeplerianOrbit orb = good.estimate(FramesFactory.getFrame(frame), gspos[0], gspos[1], gspos[2], los[0], time[0],
				los[1], time[1], los[2], time[2], rho1init, rho3init);
		return (orb);
	}

	public static ArrayList<ArrayList<Measurements.Measurement>> importTDM(String fileName, FileFormat fileFormat) {
		TdmParser parser = new ParserBuilder().buildTdmParser();
		parser.reset(fileFormat);

		Measurements.Measurement obj = null;
		ArrayList<ArrayList<Measurements.Measurement>> output = new ArrayList<ArrayList<Measurements.Measurement>>();

		for (Segment<TdmMetadata, ObservationsBlock> blk : parser.parseMessage(new DataSource(fileName))
				.getSegments()) {
			int i = 0;
			AngleType angleType = blk.getMetadata().getAngleType();
			ArrayList<Measurements.Measurement> allMeas = new ArrayList<Measurements.Measurement>();

			for (Observation obs : blk.getData().getObservations()) {
				ObservationType key = obs.getType();
				if (key != ObservationType.RANGE && key != ObservationType.DOPPLER_INSTANTANEOUS
						&& key != ObservationType.ANGLE_1 && key != ObservationType.ANGLE_2)
					continue;
				if (i == 0) {
					obj = new Measurements.Measurement();
					obj.values = new double[2];
				}

				if (angleType == null) {
					if (key == ObservationType.RANGE)
						obj.values[0] = obs.getMeasurement();
					else if (key == ObservationType.DOPPLER_INSTANTANEOUS)
						obj.values[1] = obs.getMeasurement();
				} else if (angleType == AngleType.RADEC) {
					if (key == ObservationType.ANGLE_1)
						obj.values[0] = obs.getMeasurement();
					else if (key == ObservationType.ANGLE_2)
						obj.values[1] = obs.getMeasurement();
				} else if (angleType == AngleType.AZEL) {
					if (key == ObservationType.ANGLE_1)
						obj.values[0] = obs.getMeasurement();
					else if (key == ObservationType.ANGLE_2)
						obj.values[1] = obs.getMeasurement();
				}

				if (++i == 2) {
					i = 0;
					obj.time = obs.getEpoch();
					allMeas.add(obj);
				}
			}
			output.add(allMeas);
		}

		return (output);
	}
}
