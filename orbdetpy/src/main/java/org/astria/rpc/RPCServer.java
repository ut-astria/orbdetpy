/*
 * RPCServer.java - RPC server entry point.
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

import io.grpc.Server;
import io.grpc.ServerBuilder;
import org.astria.DataManager;
import org.hipparchus.util.FastMath;
import java.nio.file.Paths;

public final class RPCServer
{
    private int port;
    private int poolSize;
    private String dataPath;
    private Server server;

    private RPCServer(int port, String dataPath, int poolSize)
    {
	this.port = port;
	this.dataPath = dataPath;
	this.poolSize = poolSize;
    }

    private void start() throws Exception
    {
	Package pack = Package.getPackage("org.astria.rpc");
	if (pack != null && pack.getImplementationTitle() != null && pack.getImplementationVersion() != null)
	    System.out.println(String.format("%s version %s", pack.getImplementationTitle(), pack.getImplementationVersion()));

	DataManager.initialize(dataPath, poolSize);
	server = ServerBuilder.forPort(port).maxInboundMessageSize(Integer.MAX_VALUE).addService(new ConversionService())
	    .addService(new EstimationService()).addService(new PropagationService()).addService(new UtilitiesService()).build().start();
	Runtime.getRuntime().addShutdownHook(new Thread()
	{
	    @Override public void run()
	    {
		System.out.println("Shutting down server...");
		RPCServer.this.stop();
	    }
	});
    }

    private void block() throws Exception
    {
	if (server != null)
	    server.awaitTermination();
    }

    private void stop()
    {
	if (server != null)
	    server.shutdown();
    }

    public static void main(String[] args) throws Exception
    {
	String dataPath = Paths.get(".").toAbsolutePath().toString();
	int port = 50051, poolSize = FastMath.min(Runtime.getRuntime().availableProcessors(), 64);
	if (args.length > 0)
	    port = Integer.parseInt(args[0]);
	if (args.length > 1)
	    dataPath = args[1];
	if (args.length > 2)
	    poolSize = Integer.parseInt(args[2]);

	RPCServer server = new RPCServer(port, dataPath, poolSize);
	server.start();
	server.block();
    }
}
