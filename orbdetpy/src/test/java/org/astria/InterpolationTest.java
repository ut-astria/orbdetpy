/*
 * InterpolationTest.java - Compare interpolation functions.
 * Copyright (C) 2023 University of Texas
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

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import org.astria.DataManager;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.data.DataContext;
import org.orekit.data.DataSource;
import org.orekit.files.ccsds.ndm.ParserBuilder;
import org.orekit.files.ccsds.ndm.odm.oem.Oem;
import org.orekit.files.ccsds.ndm.odm.oem.OemParser;
import org.orekit.files.ccsds.ndm.odm.oem.OemSatelliteEphemeris;
import org.orekit.frames.FramesFactory;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.BoundedPropagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.Ephemeris;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.CartesianDerivativesFilter;
import org.orekit.utils.Constants;
import org.orekit.utils.TimeStampedPVCoordinates;

public class InterpolationTest {
	public static void main(String[] args) {
		try {
			DataManager.initialize("./orekit-data");

			OemParser parser = new ParserBuilder(DataContext.getDefault()).buildOemParser();
			Oem oem = parser.parse(new DataSource(new File(args[0])));
			OemSatelliteEphemeris sat = oem.getSatellites().entrySet().iterator().next().getValue();
			BoundedPropagator propagator = sat.getPropagator();

			AbsoluteDate time = new AbsoluteDate(DateTimeComponents.parseDateTime(args[1]), TimeScalesFactory.getUTC());

			TimeStampedPVCoordinates pv = propagator.getPVCoordinates(time, FramesFactory.getEME2000());
			Vector3D pos = pv.getPosition();
			Vector3D vel = pv.getVelocity();
			System.out.println(
					"Oem interpolator " + Arrays.toString(pos.toArray()) + " " + Arrays.toString(vel.toArray()));

			List<TimeStampedPVCoordinates> oemStates = sat.getSegments().get(0).getCoordinates();
			AbsoluteDate[] oemTimes = new AbsoluteDate[oemStates.size()];

			for (int i = 0; i < oemStates.size(); i++) {
				oemTimes[i] = oemStates.get(i).getDate();
			}

			int index = Arrays.binarySearch(oemTimes, time);
			if (index < 0) {
				index = -index - 2;
			}

			pv = TimeStampedPVCoordinates.interpolate(time, CartesianDerivativesFilter.USE_PV,
					oemStates.subList(index - 2, index + 3));
			pos = pv.getPosition();
			vel = pv.getVelocity();
			System.out.println("TimeStampedPVCoordinates interpolator " + Arrays.toString(pos.toArray()) + " "
					+ Arrays.toString(vel.toArray()));

			ArrayList<SpacecraftState> ss = new ArrayList<SpacecraftState>(5);

			for (int i = index - 2; i <= index + 2; i++) {
				ss.add(new SpacecraftState(
						new CartesianOrbit(oemStates.get(i), FramesFactory.getEME2000(), Constants.EGM96_EARTH_MU)));
			}

			Ephemeris ephemeris = new Ephemeris(ss, 5);
			pv = ephemeris.getPVCoordinates(time, FramesFactory.getEME2000());
			pos = pv.getPosition();
			vel = pv.getVelocity();
			System.out.println(
					"Ephemeris interpolator " + Arrays.toString(pos.toArray()) + " " + Arrays.toString(vel.toArray()));
		} catch (Exception exc) {
			exc.printStackTrace(System.out);
		}

		System.exit(0);
	}
}
