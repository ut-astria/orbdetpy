# interpolate_oem.py - Interpolate CCSDS OEM files.
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

import sys
import argparse
from orbdetpy import Frame
from orbdetpy.conversion import get_J2000_epoch_offset, get_UTC_string
from orbdetpy.utilities import interpolate_ephemeris

parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument("oem-file", help="CCSDS OEM file", type=argparse.FileType("r"))
parser.add_argument("step-size", help="step size [sec]", type=float)
parser.add_argument("degree", help="interpolating polynomial degree", type=int, nargs="?", default=5)
if (len(sys.argv) == 1):
    parser.print_help()
    exit(1)

args = parser.parse_args()
step = getattr(args, "step-size")
degree = getattr(args, "degree")

# Read OEM file, skipping over blank lines and comments
lines = [l for l in getattr(args, "oem-file").read().splitlines() if (len(l) > 0 and not l.startswith("COMMENT"))]

# *Rough* mapping of CCSDS reference frames to their Orekit counterparts 
frame_map = {"EME2000": Frame.EME2000, "GCRF": Frame.GCRF, "ICRF": Frame.ICRF, "ITRF2000": Frame.ITRF_CIO_CONV_2003_ACCURATE_EOP,
             "ITRF-93": Frame.ITRF_CIO_CONV_1996_ACCURATE_EOP, "ITRF-97": Frame.ITRF_CIO_CONV_1996_ACCURATE_EOP,
             "TEME": Frame.TEME, "TOD": Frame.TOD_CONVENTIONS_2010_ACCURATE_EOP}

body, utc, states = False, [], []
ref_frame, time_system, start_utc, final_utc = [None]*4
for l in lines:
    # Read and parse OEM header block entries
    if (l.startswith("REF_FRAME")):
        ref_frame = frame_map.get(l.split("=")[-1].strip())
        continue
    if (l.startswith("TIME_SYSTEM")):
        time_system = l.split("=")[-1].strip()
        continue
    if (l.startswith("START_TIME")):
        start_utc = l.split("=")[-1].strip()
        continue
    if (l.startswith("STOP_TIME")):
        final_utc = l.split("=")[-1].strip()
        continue
    if (l.startswith("META_STOP")):
        if (ref_frame is None or time_system != "UTC" or start_utc is None or final_utc is None):
            print("Error: unsupported OEM file format")
            exit(1)
        body = True
        continue

    # Read and parse ephemeris entries
    if (body):
        toks = l.split()
        utc.append(toks[0])
        states.append([float(t)*1000.0 for t in toks[1:]])
        if (toks[0] == final_utc):
            break

# Bulk convert UTC timestamps to TT offsets 
tt = get_J2000_epoch_offset(utc)

# Interpolate ephemeris
for i in interpolate_ephemeris(ref_frame, tt, states, degree, ref_frame, tt[0], tt[-1], step):
    print(get_UTC_string(i.time), [x/1000.0 for x in i.true_state])
