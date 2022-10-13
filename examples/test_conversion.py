# test_conversion.py - Test conversion functions.
# Copyright (C) 2020-2022 University of Texas
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

from datetime import datetime, timezone
import orbdetpy.conversion as conv
from orbdetpy import Constant, Epoch, Frame, PositionAngle
import numpy as np
from numpy.linalg import norm

# TT <=> UTC
now = datetime.now(timezone.utc).isoformat()
print(f"J2000_EPOCH = {conv.get_UTC_string(0.0)}")
print(f"""J2000_EPOCH => 2000-01-01T12:00:00Z = {conv.get_J2000_epoch_offset("2000-01-01T12:00:00")}s""")
print(f"""J2000_EPOCH => {now} (now) = {conv.get_J2000_epoch_offset(now)}s""")

# Time difference between pairs of epochs
for key, val in vars(Epoch).items():
    if (key.endswith("_EPOCH") and key != "J2000_EPOCH"):
        print(f"J2000_EPOCH => {key} = {conv.get_epoch_difference(Epoch.J2000_EPOCH, val)}s")
print()

time = conv.get_J2000_epoch_offset("1998-08-10T23:10:00")
# Inertial GCRF state vector
gcrf_pv = [-6045E3, -3490E3, 2500E3, -3.457E3, 6.618E3, 2.533E3]

# Frame transformations
itrf_pv = conv.transform_frame(Frame.GCRF, time, gcrf_pv, Frame.ITRF_CIO_CONV_2010_ACCURATE_EOP)
print("GCRF=>ITRF: {:.2f}m, {:.2f}m, {:.2f}m, {:.2f}m/s, {:.2f}m/s, {:.2f}m/s".format(*itrf_pv))
print("ITRF=>GCRF: {:.2f}m, {:.2f}m, {:.2f}m, {:.2f}m/s, {:.2f}m/s, {:.2f}m/s\n".format(
    *conv.transform_frame(Frame.ITRF_CIO_CONV_2010_ACCURATE_EOP, time, itrf_pv, Frame.GCRF)))

print("converting from GCRF to ICRF")
icrf_pv = conv.transform_frame(Frame.GCRF, time, gcrf_pv, Frame.ICRF)
print("GCRF=>ICRF: {:.2f}m, {:.2f}m, {:.2f}m, {:.2f}m/s, {:.2f}m/s, {:.2f}m/s".format(*icrf_pv))
print("ICRF=>GCRF: {:.2f}m, {:.2f}m, {:.2f}m, {:.2f}m/s, {:.2f}m/s, {:.2f}m/s\n".format(
    *conv.transform_frame(Frame.ICRF, time, icrf_pv, Frame.GCRF)))

print("pos norm icrf: ", norm(np.array(icrf_pv[:3])))
print("vel norm icrf: ", norm(np.array(icrf_pv[3:])))
print("pos norm gcrf: ", norm(np.array(gcrf_pv[:3])))
print("vel norm gcrf: ", norm(np.array(gcrf_pv[3:])))
exit()


# Position => Lat/Lon/Alt
lla = conv.pos_to_lla(Frame.GCRF, time, gcrf_pv[:3])
print(f"pos_to_lla: Lat={lla[0]/Constant.DEGREE_TO_RAD:.2f}d, Lon={lla[1]/Constant.DEGREE_TO_RAD:.2f}d, Alt={lla[2]:.2f}m")

# Lat/Lon/Alt => Position
print("lla_to_pos: {:.2f}m, {:.2f}m, {:.2f}m\n".format(*conv.lla_to_pos(time, lla)))

# State vector <=> orbital elements
elem = conv.pv_to_elem(Frame.GCRF, time, gcrf_pv)
print("pv_to_elem: a={:.2f}m, e={:.4f}, i={:.2f}d, O={:.2f}d, w={:.2f}d, M={:.2f}d, theta={:.2f}d, E={:.2f}d".format(
    elem[0], elem[1], *[e/Constant.DEGREE_TO_RAD for e in elem[2:]]))
print("elem_to_pv: {:.2f}m, {:.2f}m, {:.2f}m, {:.2f}m/s, {:.2f}m/s, {:.2f}m/s\n".format(
    *conv.elem_to_pv(Frame.GCRF, 0.0, *elem[:-2], PositionAngle.MEAN)))

# RA/Dec <=> Az/El
lat, lon, alt = 52.5*Constant.DEGREE_TO_RAD, -1.917*Constant.DEGREE_TO_RAD, 0.0
ra, dec = 250.425*Constant.DEGREE_TO_RAD, 36.467*Constant.DEGREE_TO_RAD
az, el = conv.radec_to_azel(Frame.GCRF, time, ra, dec, lat, lon, alt)
print(f"radec_to_azel: Azimuth={az/Constant.DEGREE_TO_RAD:.3f}d, Elevation={el/Constant.DEGREE_TO_RAD:.3f}d")

r, d = conv.azel_to_radec(time, az, el, lat, lon, alt, Frame.GCRF)
print(f"azel_to_radec: RA={r/Constant.DEGREE_TO_RAD:.3f}d, Dec={d/Constant.DEGREE_TO_RAD:.3f}d\n")
