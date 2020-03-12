/*
 * MultiTargetService.java - Multiple target estimation.
 * Copyright (C) 2020 University of Texas
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
import org.astria.MultiTargetEstimation;
import org.astria.Measurements;
import org.astria.Settings;

public final class MultiTargetService extends MultiTargetGrpc.MultiTargetImplBase
{
    @Override public void determineOrbit(Messages.MultiTargetInput req, StreamObserver<Messages.EstimationOutputArray> resp)
    {
	try
	{
	    ArrayList<Settings> cfgList = new ArrayList<Settings>(req.getConfigCount());
	    for (int i = 0; i < req.getConfigCount(); i++)
		cfgList.add(Tools.buildSettingsFromRequest(req.getConfig(i)));

	    ArrayList<Measurements> obsList = new ArrayList<Measurements>(req.getMeasurementsCount());
	    for (int i = 0; i < req.getMeasurementsCount(); i++)
		obsList.add(Tools.buildMeasurementsFromRequest(req.getMeasurements(i).getArrayList(), cfgList.get(i)));

	    ArrayList<Estimation.EstimationOutput> estOut = new MultiTargetEstimation(cfgList, obsList).multiTargetDetermineOrbit();

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
