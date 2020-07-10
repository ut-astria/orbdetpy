/*
 * UtilitiesService.java - Utilities service handler.
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

import java.util.ArrayList;
import io.grpc.Status;
import io.grpc.StatusRuntimeException;
import io.grpc.stub.StreamObserver;
import org.astria.Measurements;
import org.astria.Utilities;
import org.orekit.files.ccsds.TDMParser.TDMFileFormat;
import org.orekit.frames.Predefined;
import org.orekit.time.AbsoluteDate;

public final class UtilitiesService extends UtilitiesGrpc.UtilitiesImplBase
{
    @Override public void importTDM(Messages.ImportTDMInput req, StreamObserver<Messages.Measurement2DArray> resp)
    {
	try
	{
	    ArrayList<ArrayList<Measurements.Measurement>> mlist = Utilities.importTDM(req.getFileName(), TDMFileFormat.values()[req.getFileFormat()]);

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
	    ArrayList<Double[]> ephem = new ArrayList<Double[]>(req.getEphemCount());
	    for (int i = 0; i < req.getEphemCount(); i++)
		ephem.add(req.getEphem(i).getArrayList().toArray(new Double[0]));

	    ArrayList<AbsoluteDate> times = new ArrayList<AbsoluteDate>(req.getTimeCount());
	    for (int i = 0; i < req.getTimeCount(); i++)
		times.add(AbsoluteDate.J2000_EPOCH.shiftedBy(req.getTime(i)));

	    ArrayList<Measurements.Measurement> interp = Utilities.interpolateEphemeris(
		Predefined.valueOf(req.getSourceFrame()), times, ephem, req.getNumPoints(), Predefined.valueOf(req.getDestFrame()),
		AbsoluteDate.J2000_EPOCH.shiftedBy(req.getInterpStart()),
		AbsoluteDate.J2000_EPOCH.shiftedBy(req.getInterpEnd()), req.getStepSize());

	    Messages.MeasurementArray.Builder builder = Messages.MeasurementArray.newBuilder()
		.addAllArray(Tools.buildResponseFromMeasurements(interp));
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }
}
