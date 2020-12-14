# propagation.py - Orbit propagation functions.
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

from typing import List
from orbdetpy.rpc.messages_pb2 import Settings, SettingsArray
from orbdetpy.rpc.propagation_pb2_grpc import PropagationStub
from orbdetpy.rpc.server import RemoteServer

def propagate_orbits(cfg_list: List[Settings]):
    """Propagate orbits and optionally simulate measurements.

    Parameters
    ----------
    cfg_list: List of Settings objects.

    Returns
    -------
    Propagated state vectors and simulated measurements.
    """

    resp = _propagation_stub.propagate(SettingsArray(array=[p for p in cfg_list]))
    return(resp.array)

if (__name__ != '__main__'):
    __pdoc__ = {m: False for m in ("Settings", "SettingsArray")}
    _propagation_stub = PropagationStub(RemoteServer.channel())
