# test_conversion.py - Test conversion functions.
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

from orbdetpy import Constant, Frame, PositionAngle
from orbdetpy.conversion import (azel_to_radec, elem_to_pv, get_J2000_epoch_offset, get_UTC_string,
                                 pos_to_lla, pv_to_elem, radec_to_azel, transform_frame)

# TT <=> UTC
# NOTE: J2000.0 epoch is "2000-01-01T12:00:00 TT" and NOT UTC!
print("UTC at J2000.0 epoch = {}".format(get_UTC_string(0.0)))
print("Time offset of UTC 2000-01-01T12:00:00Z from J2000.0 epoch = {}s\n".format(
    get_J2000_epoch_offset("2000-01-01T12:00:00Z")))

time = get_J2000_epoch_offset("1998-08-10T23:10:00Z")
# Inertial GCRF state vector
gcrf_pv = [-6045E3, -3490E3, 2500E3, -3.457E3, 6.618E3, 2.533E3]

# Frame transformations
itrf_pv = transform_frame(Frame.GCRF, time, gcrf_pv, Frame.ITRF_CIO_CONV_2010_ACCURATE_EOP)
print("GCRF->ITRF: {:.2f}m, {:.2f}m, {:.2f}m, {:.2f}m/s, {:.2f}m/s, {:.2f}m/s".format(*itrf_pv))
print("ITRF->GCRF: {:.2f}m, {:.2f}m, {:.2f}m, {:.2f}m/s, {:.2f}m/s, {:.2f}m/s\n".format(
    *transform_frame(Frame.ITRF_CIO_CONV_2010_ACCURATE_EOP, time, itrf_pv, Frame.GCRF)))

# Position => Lat/Lon/Alt
lla = pos_to_lla(Frame.GCRF, time, gcrf_pv[:3])
print("pos_to_lla: Lat={:.2f}d, Lon={:.2f}d, Alt={:.2f}m\n".format(
    lla[0]/Constant.DEGREE_TO_RAD, lla[1]/Constant.DEGREE_TO_RAD, lla[2]))

# State vector <=> orbital elements
elem = pv_to_elem(Frame.GCRF, time, gcrf_pv)
print("pv_to_elem: a={:.2f}m, e={:.4f}, i={:.2f}d, O={:.2f}d, w={:.2f}d, M={:.2f}d, theta={:.2f}d, E={:.2f}d".format(
    elem[0], elem[1], *[e/Constant.DEGREE_TO_RAD for e in elem[2:]]))
print("elem_to_pv: {:.2f}m, {:.2f}m, {:.2f}m, {:.2f}m/s, {:.2f}m/s, {:.2f}m/s\n".format(
    *elem_to_pv(Frame.GCRF, 0.0, *elem[:-2], PositionAngle.MEAN)))

# RA/Dec <=> Az/El
lat, lon, alt = 52.5*Constant.DEGREE_TO_RAD, -1.917*Constant.DEGREE_TO_RAD, 0.0
ra, dec = 250.425*Constant.DEGREE_TO_RAD, 36.467*Constant.DEGREE_TO_RAD
az, el = radec_to_azel(Frame.GCRF, time, ra, dec, lat, lon, alt)
print("radec_to_azel: Azimuth={:.3f}d, Elevation={:.3f}d".format(az/Constant.DEGREE_TO_RAD, el/Constant.DEGREE_TO_RAD))

r, d = azel_to_radec(time, az, el, lat, lon, alt, Frame.GCRF)
print("azel_to_radec: RA={:.3f}d, Dec={:.3f}d\n".format(r/Constant.DEGREE_TO_RAD, d/Constant.DEGREE_TO_RAD))
