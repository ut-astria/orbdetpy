# tle_propagation.py - Validate TLE propagation with the test cases
# in the paper "Revisiting Spacetrack Report #3" by Vallado et al.
# Copyright (C) 2021-2023 University of Texas
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

import math
from os import path
from datetime import datetime, timedelta
from orbdetpy import configure, DragModel, Frame
from orbdetpy.conversion import get_J2000_epoch_offset
from orbdetpy.propagation import propagate_orbits

with open(path.join(path.dirname(path.realpath(__file__)), "tle_tests.txt"), "r") as fp:
    test_cases = [[]]
    for l in fp.read().splitlines():
        if (l.strip()):
            test_cases[-1].append(l)
        else:
            test_cases.append([])

config = [configure(prop_step=60.0, prop_inertial_frame=Frame.TEME, gravity_degree=-1, gravity_order=-1, ocean_tides_degree=-1,
                    ocean_tides_order=-1, third_body_sun=False, third_body_moon=False, solid_tides_sun=False,
                    solid_tides_moon=False, drag_model=DragModel.UNDEFINED, rp_sun=False)]
diff_norm = lambda x, y: math.sqrt((x[0]-y[0])**2 + (x[1]-y[1])**2 + (x[2]-y[2])**2)

for test in test_cases:
    epstr = "19" + test[0][18:32] if (test[0][18:20] >= "57") else "20" + test[0][18:32]
    eputc = (datetime.strptime(epstr[:7], "%Y%j") + timedelta(days=float(epstr[7:]))).strftime("%Y-%m-%dT%H:%M:%S.%f")
    eptdt = get_J2000_epoch_offset(eputc)
    config[0].prop_initial_TLE[:] = test[:2]

    for line in test[2:]:
        truth = [float(t)*1000.0 for t in line.split()[:7]]
        truth[0] = truth[0]*60.0/1000.0
        config[0].prop_start = eptdt + truth[0]
        config[0].prop_end = config[0].prop_start
        prop = propagate_orbits(config)[0].array[0]
        assert(diff_norm(truth[1:4], prop.true_state[:3]) < 2.0E-3)     # 2.0 mm position error tolerance
        assert(diff_norm(truth[4:7], prop.true_state[3:6]) < 0.003E-3)  # 0.003 mm/s velocity error tolerance
