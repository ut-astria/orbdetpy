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
import org.orekit.frames.Predefined;
import org.orekit.time.AbsoluteDate;

public final class ConversionService extends ConversionGrpc.ConversionImplBase
{
    @Override public void transformFrame(Messages.TransformFrameInput req, StreamObserver<Messages.DoubleArray> resp)
    {
	try
	{
	    double[] pva = Conversion.transformFrame(Predefined.valueOf(req.getSrcFrame()), AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime()),
						     req.getPvaList(), Predefined.valueOf(req.getDestFrame()));
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
	    double[] raDec = Conversion.convertAzElToRaDec(AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(0)), req.getAngle1(0),
							   req.getAngle2(0), req.getLatitude(), req.getLongitude(),
							   req.getAltitude(), Predefined.valueOf(req.getFrame()));
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
							  req.getAngle1(0), req.getAngle2(0), req.getLatitude(),
							  req.getLongitude(), req.getAltitude());
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
	    DoubleValue.Builder builder = DoubleValue.newBuilder().setValue(Conversion.getJ2000EpochOffset(req.getValue()));
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }
}
