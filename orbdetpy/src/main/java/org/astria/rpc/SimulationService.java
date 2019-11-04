/*
 * SimulationService.java - Simulation service handler.
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
import org.astria.Settings;
import org.astria.Simulation;

public final class SimulationService extends SimulationGrpc.SimulationImplBase
{
    @Override public void simulateMeasurements(Messages.Settings req, StreamObserver<Messages.MeasurementArray> resp)
    {
	try
	{
	    Settings simCfg = Tools.buildSettingsFromRequest(req);
	    ArrayList<Measurements.SimulatedMeasurement> mlist = new Simulation(simCfg).simulateMeasurements();

	    Messages.MeasurementArray.Builder builder = Messages.MeasurementArray
		.newBuilder().addAllArray(Tools.buildResponseFromMeasurements(mlist));
	    resp.onNext(builder.build());
	    resp.onCompleted();
	}
	catch (Throwable exc)
	{
	    resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
	}
    }
}
