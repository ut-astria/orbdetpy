/*
 * UtilitiesService.java - Utilities service handler.
 * Copyright (C) 2019 University of Texas
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

public final class UtilitiesService extends UtilitiesGrpc.UtilitiesImplBase
{
    @Override public void importTDM(Messages.ImportTDMInput req, StreamObserver<Messages.Measurement2DArray> resp)
    {
	try
	{
	    ArrayList<ArrayList<Measurements.SimulatedMeasurement>> mlist =
		Utilities.importTDM(req.getFileName(), req.getFileFormat());

	    Messages.Measurement2DArray.Builder outer = Messages.Measurement2DArray.newBuilder();
	    for (ArrayList<Measurements.SimulatedMeasurement> m : mlist)
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

    @Override public void interpolateEphemeris(Messages.InterpolateEphemerisInput req, StreamObserver<Messages.InterpolationOutputArray> resp)
    {
	try
	{
	    ArrayList<Double[]> ephem = new ArrayList<Double[]>(req.getEphemCount());
	    for (int i = 0; i < req.getTimeCount(); i++)
		ephem.add(req.getEphem(i).getArrayList().toArray(new Double[0]));

	    ArrayList<Utilities.InterpolationOutput> interp = Utilities.interpolateEphemeris(req.getSourceFrame(), req.getTimeList(), ephem, req.getNumPoints(),
											     req.getDestFrame(), req.getInterpStart(),
											     req.getInterpEnd(), req.getStepSize());

	    Messages.InterpolationOutputArray.Builder outer = Messages.InterpolationOutputArray.newBuilder();
	    for (Utilities.InterpolationOutput i: interp)
	    {
		Messages.InterpolationOutput.Builder inner = Messages.InterpolationOutput.newBuilder().setTime(i.time);
		for (int j = 0; j < 6; j++)
		    inner = inner.addState(i.state[j]);
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
}
