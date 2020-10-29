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
import java.util.ArrayList;
import java.util.List;
import org.astria.Conversion;
import org.hipparchus.geometry.euclidean.threed.Vector3D;
import org.hipparchus.util.FastMath;
import org.hipparchus.util.MathUtils;
import org.orekit.frames.Frame;
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
    @Override public void transformFrame(Messages.TransformFrameInput req, StreamObserver<Messages.Double2DArray> resp)
    {
	try
	{
	    AbsoluteDate time;
	    boolean stringTime = req.getUTCTimeCount() > 0;
	    Predefined srcFrame = Predefined.valueOf(req.getSrcFrame());
	    Predefined destFrame = Predefined.valueOf(req.getDestFrame());
	    Messages.Double2DArray.Builder outer = Messages.Double2DArray.newBuilder();

	    for (int i = 0; i < req.getPvaCount(); i++)
	    {
		if (stringTime)
		    time = new AbsoluteDate(DateTimeComponents.parseDateTime(req.getUTCTime(i)), TimeScalesFactory.getUTC());
		else
		    time = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(i));
		double[] pva = Conversion.transformFrame(srcFrame, time, req.getPva(i).getArrayList(), destFrame);

		Messages.DoubleArray.Builder inner = Messages.DoubleArray.newBuilder();
		for (int j = 0; j < pva.length; j++)
		    inner = inner.addArray(pva[j]);
		if (stringTime)
		    inner = inner.addArray(time.durationFrom(AbsoluteDate.J2000_EPOCH));
		outer = outer.addArray(inner.build());
	    }
	    resp.onNext(outer.build());
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

    @Override public void convertPosToLLA(Messages.TransformFrameInput req, StreamObserver<Messages.Double2DArray> resp)
    {
	try
	{
	    AbsoluteDate time;
	    boolean stringTime = req.getUTCTimeCount() > 0;
	    Predefined srcFrame = Predefined.valueOf(req.getSrcFrame());
	    Messages.Double2DArray.Builder outer = Messages.Double2DArray.newBuilder();

	    for (int i = 0; i < req.getPvaCount(); i++)
	    {
		if (stringTime)
		    time = new AbsoluteDate(DateTimeComponents.parseDateTime(req.getUTCTime(i)), TimeScalesFactory.getUTC());
		else
		    time = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(i));
		double[] lla = Conversion.convertPosToLLA(srcFrame, time, req.getPva(i).getArrayList());

		Messages.DoubleArray.Builder inner = Messages.DoubleArray.newBuilder()
		    .addArray(lla[0]).addArray(lla[1]).addArray(lla[2]);
		if (stringTime)
		    inner = inner.addArray(time.durationFrom(AbsoluteDate.J2000_EPOCH));
		outer = outer.addArray(inner.build());
	    }
	    resp.onNext(outer.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void convertElemToPv(Messages.TransformFrameInput req, StreamObserver<Messages.Double2DArray> resp)
    {
	try
	{
	    AbsoluteDate time;
	    boolean stringTime = req.getUTCTimeCount() > 0;
	    Frame srcFrame = FramesFactory.getFrame(Predefined.valueOf(req.getSrcFrame()));
	    Messages.Double2DArray.Builder outer = Messages.Double2DArray.newBuilder();

	    for (int i = 0; i < req.getPvaCount(); i++)
	    {
		if (stringTime)
		    time = new AbsoluteDate(DateTimeComponents.parseDateTime(req.getUTCTime(i)), TimeScalesFactory.getUTC());
		else
		    time = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(i));
		List<Double> pva = req.getPva(i).getArrayList();
		KeplerianOrbit kep = new KeplerianOrbit(pva.get(0), pva.get(1), pva.get(2), pva.get(4), pva.get(3), pva.get(5),
							PositionAngle.values()[pva.get(6).intValue()], srcFrame, time, Constants.EGM96_EARTH_MU);
		PVCoordinates pvc = kep.getPVCoordinates();
		double[] pos = pvc.getPosition().toArray();
		double[] vel = pvc.getVelocity().toArray();

		Messages.DoubleArray.Builder inner = Messages.DoubleArray.newBuilder()
		    .addArray(pos[0]).addArray(pos[1]).addArray(pos[2]).addArray(vel[0]).addArray(vel[1]).addArray(vel[2]);
		if (stringTime)
		    inner = inner.addArray(time.durationFrom(AbsoluteDate.J2000_EPOCH));
		outer = outer.addArray(inner.build());
	    }
	    resp.onNext(outer.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void convertPvToElem(Messages.TransformFrameInput req, StreamObserver<Messages.Double2DArray> resp)
    {
	try
	{
	    AbsoluteDate time;
	    boolean stringTime = req.getUTCTimeCount() > 0;
	    Frame srcFrame = FramesFactory.getFrame(Predefined.valueOf(req.getSrcFrame()));
	    Messages.Double2DArray.Builder outer = Messages.Double2DArray.newBuilder();

	    for (int i = 0; i < req.getPvaCount(); i++)
	    {
		if (stringTime)
		    time = new AbsoluteDate(DateTimeComponents.parseDateTime(req.getUTCTime(i)), TimeScalesFactory.getUTC());
		else
		    time = AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(i));
		List<Double> pva = req.getPva(i).getArrayList();
		PVCoordinates pvc = new PVCoordinates(new Vector3D(pva.get(0), pva.get(1), pva.get(2)), new Vector3D(pva.get(3), pva.get(4), pva.get(5)));
		KeplerianOrbit orb = new KeplerianOrbit(pvc, srcFrame, time, Constants.EGM96_EARTH_MU);

		Messages.DoubleArray.Builder inner = Messages.DoubleArray.newBuilder()
		    .addArray(orb.getA()).addArray(orb.getE()).addArray(orb.getI())
		    .addArray(MathUtils.normalizeAngle(orb.getRightAscensionOfAscendingNode(), FastMath.PI))
		    .addArray(MathUtils.normalizeAngle(orb.getPerigeeArgument(), FastMath.PI))
		    .addArray(MathUtils.normalizeAngle(orb.getMeanAnomaly(), FastMath.PI))
		    .addArray(MathUtils.normalizeAngle(orb.getTrueAnomaly(), FastMath.PI))
		    .addArray(MathUtils.normalizeAngle(orb.getEccentricAnomaly(), FastMath.PI));
		if (stringTime)
		    inner = inner.addArray(time.durationFrom(AbsoluteDate.J2000_EPOCH));
		outer = outer.addArray(inner.build());
	    }
	    resp.onNext(outer.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void getUTCString(Messages.DoubleArray req, StreamObserver<StringValue> resp)
    {
	try
	{
	    boolean truncate = req.getArray(0) != 0.0;
	    ArrayList<String> utc = new ArrayList<String>(req.getArrayCount() - 1);
	    for (int i = 1; i < req.getArrayCount(); i++)
	    {
		if (truncate)
		    utc.add(AbsoluteDate.J2000_EPOCH.shiftedBy(req.getArray(i)).toString(TimeScalesFactory.getUTC()) + "Z");
		else
		    utc.add(AbsoluteDate.J2000_EPOCH.shiftedBy(req.getArray(i)).toStringRfc3339(TimeScalesFactory.getUTC()));
	    }

	    StringValue.Builder builder = StringValue.newBuilder().setValue(String.join(" ", utc));
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }

    @Override public void getJ2000EpochOffset(StringValue req, StreamObserver<Messages.DoubleArray> resp)
    {
	try
	{
	    Messages.DoubleArray.Builder builder = Messages.DoubleArray.newBuilder();
	    for (String time: req.getValue().split(" ", 0))
		builder = builder.addArray(new AbsoluteDate(DateTimeComponents.parseDateTime(time),
							    TimeScalesFactory.getUTC()).durationFrom(AbsoluteDate.J2000_EPOCH));
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }
}
