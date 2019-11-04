/*
 * EstimationService.java - Estimation service handler.
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
import org.astria.Estimation;
import org.astria.Measurements;
import org.astria.Settings;

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
}
