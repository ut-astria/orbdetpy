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

import io.grpc.stub.StreamObserver;
import org.astria.Utilities;

public class UtilitiesService extends UtilitiesGrpc.UtilitiesImplBase
{
    @Override public void importTDM(UtilitiesRequest.ImportTDMInput inp,
				    StreamObserver<UtilitiesRequest.ImportTDMOutput> out)
    {
	String[] meas = Utilities.importTDM(inp.getFileName(), inp.getFileFormat());
	UtilitiesRequest.ImportTDMOutput.Builder builder = UtilitiesRequest.ImportTDMOutput.newBuilder();
	for (int i = 0; i < meas.length; i++)
	    builder = builder.addMeasurementsJson(meas[i]);
	out.onNext(builder.build());
	out.onCompleted();
    }
}
