# utilities.py - Various utilities.
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

import grpc
from orbdetpy.rpc import messages_pb2, utilities_pb2_grpc
from orbdetpy.rpc.server import RemoteServer

def interpolate_ephemeris(source_frame, times, states, num_points,
                          dest_frame, interp_start, interp_end, step_size):
    """ Interpolates the given state vectors.

    Args:
        source_frame: Source reference frame; a constant from the enum Frame
        times: Time strings of state vectors to interpolate
        states: State vectors to interpolate
        num_points: Number of states to use for interpolation
        dest_frame: Destination reference frame; a constant from the enum Frame
        interp_start: Interpolation start time
        interp_end: Interpolation end time
        step_size: Interpolation step size [s]

    Returns:
        Interpolated time strings and state vectors.
    """

    state_list = []
    for s in states:
        da = messages_pb2.DoubleArray()
        da.array.MergeFrom(s)
        state_list.append(da)

    with RemoteServer.channel() as chan:
        stub = utilities_pb2_grpc.UtilitiesStub(chan)
        resp = stub.interpolateEphemeris(messages_pb2.InterpolateEphemerisInput(
            source_frame=source_frame.name, time=times, ephem=state_list,
            num_points=num_points, dest_frame=dest_frame.name, interp_start=interp_start,
            interp_end=interp_end, step_size=step_size))

    return([r.time, list(r.state)] for r in resp.array)
