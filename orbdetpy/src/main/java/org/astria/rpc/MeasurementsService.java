/*
 * MeasurementsService.java - Measurements service handler.
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
import org.astria.Settings;
import org.astria.Simulation;

public class MeasurementsService extends MeasurementsGrpc.MeasurementsImplBase
{
    @Override public void generateMeasurements(MeasurementsRequest.GenerateMeasurementsInput inp,
					       StreamObserver<MeasurementsRequest.GenerateMeasurementsOutput> out)
    {
	MeasurementsRequest.GenerateMeasurementsOutput.Builder builder = MeasurementsRequest.GenerateMeasurementsOutput.newBuilder();
	builder = builder.setMeasurementsJson(new Simulation(inp.getConfigJson()).simulateMeasurements());
	out.onNext(builder.build());
	out.onCompleted();
    }
}
