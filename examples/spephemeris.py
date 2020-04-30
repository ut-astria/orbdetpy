# spephemeris.py - USSPACECOM SP ephemeris interpolation.
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

import os
import sys
import json
import argparse
from datetime import datetime, timedelta
from orbdetpy import Frame
from orbdetpy.utilities import interpolate_ephemeris

start = datetime.now()
final = start + timedelta(hours = 1)
timefmt = "%Y-%m-%dT%H:%M:%S.%fZ"

parser = argparse.ArgumentParser(formatter_class = argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("ephem-file", help="Input SP ephemeris file", type=argparse.FileType("r"))
parser.add_argument("output-file", help="File for interpolated output", type=argparse.FileType("w"))
parser.add_argument("--num-points", help="Number of interpolation points", type=int, default=24)
parser.add_argument("--start-time", help="Propagation start time", default=start.strftime(timefmt))
parser.add_argument("--end-time", help="Propagation end time", default=final.strftime(timefmt))
parser.add_argument("--step-size", help="Time step size [s]", type=float, default=60.0)

if (len(sys.argv) == 1):
    parser.print_help()
    exit(1)
args = parser.parse_args()

times, states = [], []
tokens = [[19, 34], [35, 50], [51, 66], [67, 82], [83, 98], [99, 114]]
for l in getattr(args, "ephem-file").read().splitlines()[1:]:
    times.append((datetime(year=int(l[1:5]), month=1, day=1, hour=int(l[8:10]),
                          minute=int(l[10:12]), second=int(l[12:14]), microsecond=
                          int(l[15:18])*1000) + timedelta(days = int(l[5:8])-1)).strftime(timefmt))
    states.append([float(l[t[0]:t[1]])*1000.0 for t in tokens])

interp = interpolate_ephemeris(Frame.TEME, times, states, args.num_points,
                               Frame.EME2000, args.start_time, args.end_time, args.step_size)
for i in interp:
    json.dump(i, getattr(args, "output-file"), indent = 1)
