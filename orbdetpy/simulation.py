# simulation.py - Measurement simulation functions.
# Copyright (C) 2019-2020 University of Texas
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

from orbdetpy import write_output_file
from orbdetpy.rpc import simulation_pb2_grpc
from orbdetpy.rpc.server import RemoteServer
from orbdetpy.rpc.tools import build_settings, convert_measurements

def simulate_measurements(config, output_file = None,
                          async_callback = None, async_extra = None):
    """ Simulates measurement data given a configuration.

    Args:
        config: Simulation configuration (Dictionary, file name, text
                file-like object, or JSON encoded string).
        output_file: If specified, the measurements will be written to
                     the file name or text file-like object given.
        async_callback: Callback function to invoke asynchronously when
                        results become available. Processing is done
                        synchronously by default.
        async_extra: Data to pass to the callback function in addition
                     to the simulation results.

    Returns:
        Simulated measurements in the same format as config (Dictionary 
        or JSON encoded string) or None in asynchronous mode.
    """

    def async_helper(resp):
        try:
            sim_data = convert_measurements(resp.result().array)
            if (output_file):
                write_output_file(output_file, sim_data)

            channel.close()
            if (async_callback):
                async_callback(sim_data, async_extra)
            return(sim_data)
        except Exception as exc:
            if (async_callback):
                async_callback(exc, async_extra)
            else:
                raise
        return(None)

    channel = RemoteServer.channel()
    stub = simulation_pb2_grpc.SimulationStub(channel)

    resp = stub.simulateMeasurements.future(build_settings(config))
    if (async_callback):
        resp.add_done_callback(async_helper)
        return(None)

    return(async_helper(resp))
