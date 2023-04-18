/*
 * UtilitiesService.java - Utilities service handler.
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

package org.astria.rpc;

import java.util.ArrayList;
import java.util.List;
import java.util.Map;
import io.grpc.Status;
import io.grpc.StatusRuntimeException;
import io.grpc.stub.StreamObserver;
import org.astria.Measurements;
import org.astria.MSISEInputs;
import org.astria.Settings;
import org.astria.Utilities;
import org.astria.WAM;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.orekit.bodies.CelestialBodyFactory;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.data.DataSource;
import org.orekit.files.ccsds.utils.FileFormat;
import org.orekit.files.sp3.SP3;
import org.orekit.files.sp3.SP3Parser;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.models.earth.atmosphere.Atmosphere;
import org.orekit.models.earth.atmosphere.NRLMSISE00;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.propagation.BoundedPropagator;
import org.orekit.propagation.SpacecraftState;
import org.orekit.propagation.analytical.Ephemeris;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class UtilitiesService extends UtilitiesGrpc.UtilitiesImplBase
{
    @Override public void importSP3(Messages.InterpolateEphemerisInput req, StreamObserver<Messages.Measurement2DArray> resp)
    {
        try
        {
            SP3 parser = new SP3Parser().parse(new DataSource(req.getSourceFrame()));
            Frame toFrame = FramesFactory.getFrame(Predefined.valueOf(req.getDestFrame()));
            Messages.Measurement2DArray.Builder outer = Messages.Measurement2DArray.newBuilder();
            for (Map.Entry<String, SP3.SP3Ephemeris> keyVal: parser.getSatellites().entrySet())
            {
                BoundedPropagator prop = keyVal.getValue().getPropagator();
                ArrayList<Measurements.Measurement> mlist = new ArrayList<Measurements.Measurement>(req.getInterpTimeCount());
                for (int i = 0; i < req.getInterpTimeCount(); i++)
                {
                    Measurements.Measurement meas = new Measurements.Measurement(
                        prop.getPVCoordinates(AbsoluteDate.J2000_EPOCH.shiftedBy(req.getInterpTime(i)), toFrame), null);
                    meas.station = keyVal.getKey();
                    mlist.add(meas);
                }
                Messages.MeasurementArray.Builder inner = Messages.MeasurementArray.newBuilder()
                    .addAllArray(Tools.buildResponseFromMeasurements(mlist));
                outer = outer.addArray(inner);
            }
            resp.onNext(outer.build());
            resp.onCompleted();
        }
        catch (Throwable exc)
        {
            resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
        }
    }

    @Override public void importTDM(Messages.ImportTDMInput req, StreamObserver<Messages.Measurement2DArray> resp)
    {
        try
        {
            ArrayList<ArrayList<Measurements.Measurement>> mlist = Utilities.importTDM(req.getFileName(), FileFormat.values()[req.getFileFormat()]);
            Messages.Measurement2DArray.Builder outer = Messages.Measurement2DArray.newBuilder();
            for (ArrayList<Measurements.Measurement> m: mlist)
            {
                Messages.MeasurementArray.Builder inner = Messages.MeasurementArray.newBuilder()
                    .addAllArray(Tools.buildResponseFromMeasurements(m));
                outer = outer.addArray(inner);
            }
            resp.onNext(outer.build());
            resp.onCompleted();
        }
        catch (Throwable exc)
        {
            resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
        }
    }

    @Override public void interpolateEphemeris(Messages.InterpolateEphemerisInput req, StreamObserver<Messages.MeasurementArray> resp)
    {
        try
        {
            Frame fromFrame = FramesFactory.getFrame(Predefined.valueOf(req.getSourceFrame()));
            Frame toFrame = FramesFactory.getFrame(Predefined.valueOf(req.getDestFrame()));
            ArrayList<SpacecraftState> states = new ArrayList<SpacecraftState>(req.getTimeCount());
            for (int i = 0; i < req.getTimeCount(); i++)
            {
                List<Double> pv = req.getEphem(i).getArrayList();
                TimeStampedPVCoordinates tspv = new TimeStampedPVCoordinates(AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(i)),
                                                                             new Vector3D(pv.get(0), pv.get(1), pv.get(2)),
                                                                             new Vector3D(pv.get(3), pv.get(4), pv.get(5)));
                states.add(new SpacecraftState(new CartesianOrbit(tspv, fromFrame, Constants.EGM96_EARTH_MU)));
            }

            Ephemeris interpolator = new Ephemeris(states, req.getNumPoints());
            ArrayList<Measurements.Measurement> output = new ArrayList<Measurements.Measurement>(req.getInterpTimeCount());
            for (int i = 0; i < req.getInterpTimeCount(); i++)
                output.add(new Measurements.Measurement(
                               interpolator.getPVCoordinates(AbsoluteDate.J2000_EPOCH.shiftedBy(req.getInterpTime(i)), toFrame), null));

            Messages.MeasurementArray.Builder builder =
                Messages.MeasurementArray.newBuilder().addAllArray(Tools.buildResponseFromMeasurements(output));
            resp.onNext(builder.build());
            resp.onCompleted();
        }
        catch (Throwable exc)
        {
            resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
        }
    }

    @Override public void getDensity(Messages.TransformFrameInput req, StreamObserver<Messages.DoubleArray> resp)
    {
        try
        {
            Atmosphere atmosphere;
            AbsoluteDate time = null;
            Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder();
            Settings.DragModel drag = Settings.DragModel.values()[Integer.parseInt(req.getDestFrame())];
            OneAxisEllipsoid earth = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING,
                                                          FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP));
            if (drag == Settings.DragModel.MSISE2000)
                atmosphere = new NRLMSISE00(MSISEInputs.getInstance(), CelestialBodyFactory.getSun(), earth);
            else if (drag == Settings.DragModel.WAM)
                atmosphere = WAM.getInstance();
            else
                throw(new RuntimeException("Invalid atmospheric drag model"));

            for (int i = 0; i < req.getPvaCount(); i++)
            {
                if (i < req.getUTCTimeCount())
                    time = new AbsoluteDate(DateTimeComponents.parseDateTime(req.getUTCTime(i)), TimeScalesFactory.getUTC());
                else if (i < req.getTimeCount())
                    time = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(i));

                Messages.DoubleArray lla = req.getPva(i);
                Vector3D pos = earth.transform(new GeodeticPoint(lla.getArray(0), lla.getArray(1), lla.getArray(2)));
                builder = builder.addArray(atmosphere.getDensity(time, pos, earth.getFrame()));
            }
            resp.onNext(builder.build());
            resp.onCompleted();
        }
        catch (Throwable exc)
        {
            resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
        }
    }
}
