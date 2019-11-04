# conversion.py - Various conversion functions.
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

from orbdetpy.rpc import messages_pb2, conversion_pb2_grpc
from orbdetpy.rpc.server import RemoteServer

def transform_frame(srcframe, time, pva, destframe):
    """ Transforms a state vector from one frame to another.

    Args:
        srcframe: Source reference frame ("EME2000", "GCRF",
                  "ITRF", "MOD", "TOD", or "TEME").
        time: State vector epoch (ISO-8601 formatted UTC string).
        pva: State vector to transform [m, m, m, m/s, m/s, m/s] or 
             [m, m, m, m/s, m/s, m/s, m/s^2, m/s^2, m/s^2].
        destframe: Destination reference frame ("EME2000", "GCRF",
                   "ITRF", "MOD", "TOD", or "TEME")..

    Returns:
        State vector transformed to the destination frame.
    """

    with RemoteServer.channel() as chan:
        stub = conversion_pb2_grpc.ConversionStub(chan)
        resp = stub.transformFrame(messages_pb2.TransformFrameInput(
            src_frame = srcframe, time = time,
            pva = pva, dest_frame = destframe))

    return(list(resp.array))
