# simtest.py - Test spacecraft measurement simulator.
# Copyright (C) 2018 Shiva Iyer <shiva.iyer AT utexas DOT edu>
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

import os
import sys
import time
import json

sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from orbdetpy import init
init()

if (len(sys.argv) < 3):
    print("Usage: python %s config_file output_file" % sys.argv[0])
    exit()

from orbdetpy.simdata import simulate

with open(sys.argv[1], "r") as f:
    config = json.load(f)

print("Simulation start : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
res = simulate(config)
with open(sys.argv[2], "w") as fout:
    json.dump(res, fout, indent = 1)
print("Simulation end   : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
