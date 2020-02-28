# multi_target.py - Multiple target estimation.
# Copyright (C) 2020 University of Texas
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

from traceback import format_exc
from orbdetpy import write_output_file
from orbdetpy.rpc import messages_pb2, multi_target_pb2_grpc
from orbdetpy.rpc.server import RemoteServer
from orbdetpy.rpc.tools import (build_settings, build_measurements,
                                convert_estimation)

def determine_orbit(config_list, meas_list, output_file = None):
    """ Multiple target orbit determination.

    Args:
        config_list: List of object configurations.
        meas_list: List of list of measurements, one sub-list per object.
        output_file: If specified, the orbit fit will be written to
                     the given file name or text file-like object. 

    Returns:
        Multiple target orbit determination results.
    """

    marray = []
    for m in meas_list:
        ma = messages_pb2.MeasurementArray()
        ma.array.extend(build_measurements(m))
        marray.append(ma)
    
    with RemoteServer.channel() as channel:
        stub = multi_target_pb2_grpc.MultiTargetStub(channel)
        request = stub.determineOrbit(messages_pb2.MultiTargetInput(
            config = [build_settings(c) for c in config_list],
            measurements = marray))

    fit_data = convert_estimation(request.array)
    if (output_file):
        write_output_file(output_file, fit_data)
    return(fit_data)
