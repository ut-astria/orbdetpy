/*
 * PropagationService.java - Propagation service handler.
 * Copyright (C) 2019-2023 University of Texas
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

import io.grpc.Status;
import io.grpc.StatusRuntimeException;
import io.grpc.stub.StreamObserver;
import java.util.ArrayList;
import org.astria.Measurements;
import org.astria.ParallelPropagation;
import org.astria.Settings;

public final class PropagationService extends PropagationGrpc.PropagationImplBase {
	@Override
	public void propagate(Messages.SettingsArray req, StreamObserver<Messages.Measurement2DArray> resp) {
		try {
			ArrayList<Settings> cfgObjs = new ArrayList<Settings>(req.getArrayCount());
			for (int i = 0; i < req.getArrayCount(); i++)
				cfgObjs.add(Tools.buildSettingsFromRequest(req.getArray(i)));

			ArrayList<ArrayList<Measurements.Measurement>> propOut = new ParallelPropagation(cfgObjs).propagate();

			Messages.Measurement2DArray.Builder outer = Messages.Measurement2DArray.newBuilder();
			for (ArrayList<Measurements.Measurement> p : propOut) {
				Messages.MeasurementArray.Builder inner = Messages.MeasurementArray.newBuilder()
						.addAllArray(Tools.buildResponseFromMeasurements(p));
				outer = outer.addArray(inner);
			}
			resp.onNext(outer.build());
			resp.onCompleted();
		} catch (Throwable exc) {
			resp.onError(new StatusRuntimeException(Status.INTERNAL.withDescription(Tools.getStackTrace(exc))));
		}
	}
}
