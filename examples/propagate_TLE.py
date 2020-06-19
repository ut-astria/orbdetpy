# propagate_TLE.py - Propagate many TLEs in parallel.
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

import os
import sys
import argparse
from datetime import datetime, timedelta

day0 = datetime.today()
day1 = day0 + timedelta(days = 1)
timefmt = "%Y-%m-%dT%H:%M:%S.%fZ"

parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("tle-file", help="Input TLE file", type=argparse.FileType("r"))
parser.add_argument("--start-time", help="Propagation start time",
                    default=datetime(day0.year, day0.month, day0.day).strftime(timefmt))
parser.add_argument("--end-time", help="Propagation end time",
                    default = datetime(day1.year, day1.month, day1.day).strftime(timefmt))
parser.add_argument("--step-size", help="Step size [sec.]", type=float, default=900.0)

if (len(sys.argv) == 1):
    parser.print_help()
    exit(1)
args = parser.parse_args()

from orbdetpy import configure
from orbdetpy.conversion import get_J2000_epoch_offset
from orbdetpy.propagation import propagate_orbits

start = get_J2000_epoch_offset(args.start_time)
end = get_J2000_epoch_offset(args.end_time)
config, elements = [], getattr(args, "tle-file").read().splitlines()
for i in range(0, len(elements)-2, 3):
    config.append(configure(prop_start=start, prop_initial_TLE=elements[i+1:i+3],
                            prop_end=end, prop_step=args.step_size))

try:
    for o in propagate_orbits(config):
        for i in o.array:
            print(i.time, i.true_state)
except Exception as exc:
    print(exc)
