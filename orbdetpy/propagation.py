# propagation.py - Orbit propagation functions.
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
from orbdetpy.rpc import messages_pb2,  propagation_pb2_grpc
from orbdetpy.rpc.server import RemoteServer
from orbdetpy.rpc.tools import build_settings, convert_propagation

def propagate_orbits(config_list, output_file = None):
    """ Propagates the given objects in parallel.

    Args:
        config_list: List of configurations (each a dictionary,
                     file name, text file-like object, or
                     JSON encoded string).
        output_file: If specified, the output will be written to
                     the file name or text file-like object given. 

    Returns:
        Propagated state vectors at desired time intervals
        (List of dictionaries).
    """

    with RemoteServer.channel() as chan:
        stub = propagation_pb2_grpc.PropagationStub(chan)
        resp = stub.propagate(messages_pb2.SettingsArray(
            array = [build_settings(p) for p in config_list]))

    propdata = convert_propagation(resp.array)
    if (output_file):
        write_output_file(output_file, propdata)
    return(propdata)
