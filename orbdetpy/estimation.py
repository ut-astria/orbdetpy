# estimation.py - Orbit estimation functions.
# Copyright (C) 2019 University of Texas
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
from orbdetpy.rpc import messages_pb2, estimation_pb2_grpc
from orbdetpy.rpc.server import RemoteServer
from orbdetpy.rpc.tools import (build_settings, build_measurements,
                                convert_estimation)

def determine_orbit(config, meas, output_file = None):
    """ Performs orbit determination given config and measurements.

    Args:
        config: OD configuration (Dictionary, file name, text
                file-like object, or JSON encoded string).
        meas: List of measurements (List, file name, text
                file-like object, or JSON encoded string).
        output_file: If specified, the orbit fit will be written to
                     the file name or text file-like object given. 

    Returns:
        Orbit determination results in the same format as config 
        (Dictionary or JSON encoded string).
    """

    with RemoteServer.channel() as chan:
        stub = estimation_pb2_grpc.EstimationStub(chan)
        resp = stub.determineOrbit(messages_pb2.DetermineOrbitInput(
            config = build_settings(config),
            measurements = build_measurements(meas)))

    fitdata = convert_estimation(resp.array)
    if (output_file):
        write_output_file(output_file, fitdata)
    return(fitdata)
