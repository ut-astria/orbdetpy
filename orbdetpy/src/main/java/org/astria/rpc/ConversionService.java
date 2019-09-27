/*
 * ConversionService.java - Conversion service handler.
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
import org.astria.Conversion;

public class ConversionService extends ConversionGrpc.ConversionImplBase
{
    @Override public void transformFrame(ConversionRequest.TransformFrameInput inp,
					 StreamObserver<ConversionRequest.TransformFrameOutput> out)
    {
	double[] pva = new double[inp.getPvaCount()];
	for (int i = 0; i < pva.length; i++)
	    pva[i] = inp.getPva(i);
	pva = Conversion.transformFrame(inp.getSrcFrame(), inp.getTime(), pva, inp.getDestFrame());

	ConversionRequest.TransformFrameOutput.Builder builder = ConversionRequest.TransformFrameOutput.newBuilder();
	for (int i = 0; i < pva.length; i++)
	    builder = builder.addPva(pva[i]);
	out.onNext(builder.build());
	out.onCompleted();
    }
}
