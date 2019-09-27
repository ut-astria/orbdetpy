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

import io.grpc.stub.StreamObserver;
import org.astria.Estimation;
import org.astria.Settings;

public class EstimationService extends EstimationGrpc.EstimationImplBase
{
    @Override public void determineOrbit(EstimationRequest.DetermineOrbitInput inp,
					 StreamObserver<EstimationRequest.DetermineOrbitOutput> out)
    {
	EstimationRequest.DetermineOrbitOutput.Builder builder = EstimationRequest.DetermineOrbitOutput.newBuilder();
	builder = builder.setEstimationJson(new Estimation(inp.getConfigJson(), inp.getMeasurementsJson()).determineOrbit());
	out.onNext(builder.build());
	out.onCompleted();
    }
}
