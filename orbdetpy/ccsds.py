# ccsds.py - CCSDS file format conversion functions.
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

import json
import math
from datetime import datetime
from jnius import autoclass

def format_data(cfg, obs):
    for s in cfg["Stations"].values():
        s["Latitude"] *= 180.0/math.pi
        s["Longitude"] *= 180.0/math.pi
        s["Altitude"] /= 1000.0

    mit = cfg["Measurements"].keys()
    for o in obs:
        for i in mit:
            if (i in ["Range", "RangeRate"]):
                o[i] /= 1000.0
            else:
                o[i] *= 180.0/math.pi

def export_TDM(obj_id, cfg_file, obs_file):
    with open(cfg_file, "r") as fh:
        cfg = json.load(fh)
    with open(obs_file, "r") as fh:
        obs = json.load(fh)
    format_data(cfg, obs)

    miter = cfg["Measurements"].keys()
    if ("RightAscension" in miter and "Declination" in miter):
        frame = cfg["Propagation"].get("InertialFrame", "EME2000")
        if (frame == "GCRF"):
            frame = "ICRF"
        obstype = "ANGLE_TYPE = RADEC\nREFERENCE_FRAME = %s" % frame
        obspath = "1,2"
    if ("Azimuth" in miter and "Elevation" in miter):
        obstype = "ANGLE_TYPE = AZEL"
        obspath = "1,2"
    if ("Range" in miter):
        obspath = "2,1,2"
        if ("Azimuth" in miter and "Elevation" in miter):
            obstype = "RANGE_UNITS = km\nANGLE_TYPE = AZEL"
        else:
            obstype = "RANGE_UNITS = km"

    utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    tdm_data = "CCSDS_TDM_VERS = 1.0\nCREATION_DATE = {utc_time}\n" \
               "ORIGINATOR = UT-ASTRIA\n\n""".format(utc_time = utcnow)

    for sname, sinfo in cfg["Stations"].items():
        sensor = "%s (Lat: %f, Lon: %f, Alt: %f km)" % (
            sname, sinfo["Latitude"], sinfo["Longitude"], sinfo["Altitude"])
        block_header = "META_START\nTIME_SYSTEM = UTC\n" \
                       "PARTICIPANT_1 = {part1}\nPARTICIPANT_2 = {sensor}\n" \
                       "MODE = SEQUENTIAL\nPATH = {part_path}\n{obs_opt}\nMETA_STOP\n\n".format(
                           part1=obj_id, part_path=obspath, obs_opt=obstype, sensor=sensor)

        block_ent = "DATA_START\n"
        for o in obs:
            if (o["Station"] != sname):
                continue
            if ("Range" in o):
                block_ent += "RANGE = {Time} {Range}\n".format(**o)
                if ("RangeRate" in o):
                    block_ent += "DOPPLER_INSTANTANEOUS = {Time} {RangeRate}\n".format(**o)
            if ("Azimuth" in o and "Elevation" in o):
                block_ent += "ANGLE_1 = {Time} {Azimuth}\n" \
                                "ANGLE_2 = {Time} {Elevation}\n".format(**o)
            if ("RightAscension" in o and "Declination" in o):
                block_ent += "ANGLE_1 = {Time} {RightAscension}\n" \
                                "ANGLE_2 = {Time} {Declination}\n".format(**o)

        tdm_data += block_header + block_ent + "DATA_STOP\n\n"

    return(tdm_data)

def import_TDM(file_name):
    return(autoclass("org.astria.Utilities").importTDM(file_name))
