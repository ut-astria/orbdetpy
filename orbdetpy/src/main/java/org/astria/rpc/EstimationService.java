/*
 * EstimationService.java - Estimation service handler.
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
import io.grpc.Status;
import io.grpc.StatusRuntimeException;
import io.grpc.stub.StreamObserver;
import org.astria.Estimation;
import org.astria.Measurements;
import org.astria.MultiTargetEstimation;
import org.astria.Settings;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.orekit.bodies.GeodeticPoint;
import org.orekit.bodies.OneAxisEllipsoid;
import org.orekit.estimation.iod.IodLaplace;
import org.orekit.estimation.measurements.GroundStation;
import org.orekit.frames.Frame;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.frames.TopocentricFrame;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.time.AbsoluteDate;
import org.orekit.utils.Constants;
import org.orekit.utils.TimeStampedPVCoordinates;

public final class EstimationService extends EstimationGrpc.EstimationImplBase
{
    @Override public void determineOrbit(Messages.DetermineOrbitInput req, StreamObserver<Messages.EstimationOutputArray> resp)
    {
        try
        {
            Settings odCfg = Tools.buildSettingsFromRequest(req.getConfig());
            Measurements odObs = Tools.buildMeasurementsFromRequest(req.getMeasurementsList(), odCfg);
            ArrayList<Estimation.EstimationOutput> estOut = new Estimation(odCfg, odObs).determineOrbit();

            Messages.EstimationOutputArray.Builder builder = Messages.EstimationOutputArray
                .newBuilder().addAllArray(Tools.buildResponseFromOrbitDetermination(estOut));
            resp.onNext(builder.build());
            resp.onCompleted();
        }
        catch (Throwable exc)
        {
            resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
        }
    }

    @Override public void multiTargetOD(Messages.MultiTargetInput req, StreamObserver<Messages.MultiTargetOutput> resp)
    {
        try
        {
            ArrayList<Settings> cfgList = new ArrayList<Settings>(req.getConfigCount());
            for (int i = 0; i < req.getConfigCount(); i++)
                cfgList.add(Tools.buildSettingsFromRequest(req.getConfig(i)));

            ArrayList<Measurements> obsList = new ArrayList<Measurements>(req.getMeasurementsCount());
            for (int i = 0; i < req.getMeasurementsCount(); i++)
                obsList.add(Tools.buildMeasurementsFromRequest(req.getMeasurements(i).getArrayList(), cfgList.get(i)));

            MultiTargetEstimation.MultiTargetOutput multiOut = new MultiTargetEstimation(cfgList, obsList).multiTargetDetermineOrbit();
            Messages.MultiTargetOutput.Builder outer = Messages.MultiTargetOutput.newBuilder();
            for (ArrayList<Estimation.EstimationOutput> a: multiOut.estOutput)
            {
                Messages.EstimationOutputArray.Builder inner = Messages.EstimationOutputArray
                    .newBuilder().addAllArray(Tools.buildResponseFromOrbitDetermination(a));
                outer = outer.addEstOutput(inner);
            }

            for (ArrayList<Integer> a: multiOut.associatedObs)
            {
                Messages.IntegerArray.Builder inner = Messages.IntegerArray.newBuilder().addAllArray(a);
                outer = outer.addAssociatedObs(inner);
            }

            outer = outer.addAllUnassociatedObs(multiOut.unassociatedObs);
            resp.onNext(outer.build());
            resp.onCompleted();
        }
        catch (Throwable exc)
        {
            resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
        }
    }

    @Override public void iodLaplace(Messages.AnglesInput req, StreamObserver<Messages.DoubleArray> resp)
    {
        try
        {
            Vector3D[] los = new Vector3D[3];
            AbsoluteDate[] time = new AbsoluteDate[3];
            for (int i = 0; i < 3; i++)
            {
                time[i] = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(i));
                los[i] = new Vector3D(FastMath.cos(req.getAngle2(i))*FastMath.cos(req.getAngle1(i)),
                                      FastMath.cos(req.getAngle2(i))*FastMath.sin(req.getAngle1(i)),
                                      FastMath.sin(req.getAngle2(i)));
            }

            OneAxisEllipsoid body = new OneAxisEllipsoid(Constants.WGS84_EARTH_EQUATORIAL_RADIUS, Constants.WGS84_EARTH_FLATTENING,
                                                         FramesFactory.getFrame(Predefined.ITRF_CIO_CONV_2010_ACCURATE_EOP));
            GroundStation gs = new GroundStation(new TopocentricFrame(body, new GeodeticPoint(req.getLatitude(),
                                                                                              req.getLongitude(), req.getAltitude()), "sensor"));
            gs.getPrimeMeridianOffsetDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
            gs.getPolarOffsetXDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);
            gs.getPolarOffsetYDriver().setReferenceDate(AbsoluteDate.J2000_EPOCH);

            Frame frame = FramesFactory.getFrame(Predefined.valueOf(req.getFrame()));
            TimeStampedPVCoordinates obsPv = gs.getBaseFrame().getPVCoordinates(time[1], frame);
            CartesianOrbit estOrbit = new IodLaplace(Constants.EGM96_EARTH_MU).estimate(
                frame, obsPv, time[0], los[0], time[1], los[1], time[2], los[2]);
            TimeStampedPVCoordinates estPv = estOrbit.getPVCoordinates();
            double[] estPos = estPv.getPosition().toArray();
            double[] estVel = estPv.getVelocity().toArray();

            Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder();
            for (int i = 0; i < 6; i++)
            {
                if (i < 3)
                    builder = builder.addArray(estPos[i]);
                else
                    builder = builder.addArray(estVel[i-3]);
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
