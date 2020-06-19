# test_interpolation.py - Test ephemeris interpolation.
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

from orbdetpy import configure, Frame
from orbdetpy.conversion import get_J2000_epoch_offset
from orbdetpy.propagation import propagate_orbits
from orbdetpy.utilities import interpolate_ephemeris

# Propagate for 1 hour at 5 minute intervals with default settings
cfg = configure(prop_start=get_J2000_epoch_offset("2020-03-09T22:00:02.000Z"),
                prop_initial_state=[-152408.166, -958234.785, 6908448.381,
                                    -7545.691, 285.553, -126.766],
                prop_end=get_J2000_epoch_offset("2020-03-09T23:00:02.000Z"),
                prop_step=300.0)

times, states = [], []
for o in propagate_orbits([cfg])[0].array:
    times.append(o.time)
    states.append(o.true_state)

# Interpolate over the same period at 1 minute intervals
interp = interpolate_ephemeris(Frame.EME2000, times, states, len(times),
                               Frame.EME2000, cfg.prop_start, cfg.prop_end, 60.0)
for i in interp:
    print(i.time, i.true_state)
