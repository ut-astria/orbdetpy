/*
 * ConversionService.java - Conversion service handler.
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

package org.astria.rpc;

import com.google.protobuf.DoubleValue;
import com.google.protobuf.StringValue;
import io.grpc.Status;
import io.grpc.StatusRuntimeException;
import io.grpc.stub.StreamObserver;
import org.astria.Conversion;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.frames.FramesFactory;
import org.orekit.frames.Predefined;
import org.orekit.orbits.CartesianOrbit;
import org.orekit.orbits.KeplerianOrbit;
import org.orekit.orbits.PositionAngle;
import org.orekit.time.AbsoluteDate;
import org.orekit.time.DateTimeComponents;
import org.orekit.time.TimeScalesFactory;
import org.orekit.utils.Constants;
import org.orekit.utils.PVCoordinates;

public final class ConversionService extends ConversionGrpc.ConversionImplBase
{
    @Override public void transformFrame(Messages.TransformFrameInput req, StreamObserver<Messages.DoubleArray> resp)
    {
	try
	{
	    AbsoluteDate time;
	    if (req.getUTCTime().length() == 0)
		time = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime());
	    else
		time = new AbsoluteDate(DateTimeComponents.parseDateTime(req.getUTCTime()), TimeScalesFactory.getUTC());
	    double[] pva = Conversion.transformFrame(Predefined.valueOf(req.getSrcFrame()), time, req.getPvaList(), Predefined.valueOf(req.getDestFrame()));

	    Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder();
	    for (int i = 0; i < pva.length; i++)
		builder = builder.addArray(pva[i]);
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void convertAzElToRaDec(Messages.AnglesInput req, StreamObserver<Messages.DoubleArray> resp)
    {
	try
	{
	    double[] raDec = Conversion.convertAzElToRaDec(AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(0)), req.getAngle1(0), req.getAngle2(0),
							   req.getLatitude(), req.getLongitude(), req.getAltitude(), Predefined.valueOf(req.getFrame()));

	    Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder();
	    for (int i = 0; i < raDec.length; i++)
		builder = builder.addArray(raDec[i]);
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void convertRaDecToAzEl(Messages.AnglesInput req, StreamObserver<Messages.DoubleArray> resp)
    {
	try
	{
	    double[] azEl = Conversion.convertRaDecToAzEl(Predefined.valueOf(req.getFrame()), AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(0)),
							  req.getAngle1(0), req.getAngle2(0), req.getLatitude(), req.getLongitude(), req.getAltitude());

	    Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder();
	    for (int i = 0; i < azEl.length; i++)
		builder = builder.addArray(azEl[i]);
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void convertPosToLLA(Messages.TransformFrameInput req, StreamObserver<Messages.DoubleArray> resp)
    {
	try
	{
	    AbsoluteDate time;
	    if (req.getUTCTime().length() == 0)
		time = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime());
	    else
		time = new AbsoluteDate(DateTimeComponents.parseDateTime(req.getUTCTime()), TimeScalesFactory.getUTC());
	    double[] lla = Conversion.convertPosToLLA(Predefined.valueOf(req.getSrcFrame()), time, req.getPvaList());

	    Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder();
	    for (int i = 0; i < lla.length; i++)
		builder = builder.addArray(lla[i]);
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void convertElemToPv(Messages.TransformFrameInput req, StreamObserver<Messages.DoubleArray> resp)
    {
	try
	{
	    AbsoluteDate time;
	    if (req.getUTCTime().length() == 0)
		time = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime());
	    else
		time = new AbsoluteDate(DateTimeComponents.parseDateTime(req.getUTCTime()), TimeScalesFactory.getUTC());
	    KeplerianOrbit kep = new KeplerianOrbit(req.getPva(0), req.getPva(1), req.getPva(2), req.getPva(4), req.getPva(3), req.getPva(5),
						    PositionAngle.values()[(int)req.getPva(6)], FramesFactory.getFrame(Predefined.valueOf(req.getSrcFrame())),
						    time, Constants.EGM96_EARTH_MU);
	    PVCoordinates pvc = kep.getPVCoordinates();
	    double[] pos = pvc.getPosition().toArray();
	    double[] vel = pvc.getVelocity().toArray();

	    Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder()
		.addArray(pos[0]).addArray(pos[1]).addArray(pos[2]).addArray(vel[0]).addArray(vel[1]).addArray(vel[2]);
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void convertPvToElem(Messages.TransformFrameInput req, StreamObserver<Messages.DoubleArray> resp)
    {
	try
	{
	    AbsoluteDate time;
	    if (req.getUTCTime().length() == 0)
		time = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime());
	    else
		time = new AbsoluteDate(DateTimeComponents.parseDateTime(req.getUTCTime()), TimeScalesFactory.getUTC());
	    PVCoordinates pvc = new PVCoordinates(new Vector3D(req.getPva(0), req.getPva(1), req.getPva(2)),
						  new Vector3D(req.getPva(3), req.getPva(4), req.getPva(5)));
	    KeplerianOrbit orb = new KeplerianOrbit(pvc, FramesFactory.getFrame(Predefined.valueOf(req.getSrcFrame())), time, Constants.EGM96_EARTH_MU);

	    Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder()
		.addArray(orb.getA()).addArray(orb.getE()).addArray(orb.getI())
		.addArray(MathUtils.normalizeAngle(orb.getRightAscensionOfAscendingNode(), FastMath.PI))
		.addArray(MathUtils.normalizeAngle(orb.getPerigeeArgument(), FastMath.PI))
		.addArray(MathUtils.normalizeAngle(orb.getMeanAnomaly(), FastMath.PI))
		.addArray(MathUtils.normalizeAngle(orb.getTrueAnomaly(), FastMath.PI))
		.addArray(MathUtils.normalizeAngle(orb.getEccentricAnomaly(), FastMath.PI));
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void getUTCString(DoubleValue req, StreamObserver<StringValue> resp)
    {
	try
	{
	    StringValue.Builder builder = StringValue.newBuilder().setValue(Conversion.getUTCString(req.getValue()));
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void getJ2000EpochOffset(StringValue req, StreamObserver<DoubleValue> resp)
    {
	try
	{
	    DoubleValue.Builder builder = DoubleValue.newBuilder().setValue(
		new AbsoluteDate(DateTimeComponents.parseDateTime(req.getValue()), TimeScalesFactory.getUTC()).durationFrom(AbsoluteDate.J2000_EPOCH));
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }
}
