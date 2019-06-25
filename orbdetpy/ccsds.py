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
    mit = cfg["Measurements"].keys()
    for o in obs:
        for i in mit:
            if (i in ["Range", "RangeRate"]):
                o[i] /= 1000.0
            else:
                o[i] *= 180.0/math.pi

def export_TDM(obj_id_list, cfg_json_list, obs_json_list):
    utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    tdm_data = "CCSDS_TDM_VERS = 1.0\nCREATION_DATE = {utc_time}\n" \
               "ORIGINATOR = UT-ASTRIA\n\n""".format(utc_time = utcnow)

    for obj_id, cfg_json, obs_json in zip(obj_id_list, cfg_json_list, obs_json_list):
        cfg = json.loads(cfg_json)
        obs = json.loads(obs_json)
        format_data(cfg, obs)

        miter = cfg["Measurements"].keys()
        if ("RightAscension" in miter):
            obstype = "ANGLE_TYPE = RADEC\nREFERENCE_FRAME = EME2000"
            obspath = "1,2"
        elif ("Azimuth" in miter):
            obstype = "ANGLE_TYPE = AZEL"
            obspath = "1,2"
        elif ("Range" in miter):
            obstype = "RANGE_UNITS = km"
            obspath = "2,1,2"
        else:
            return(None)

        block_header = "META_START\nTIME_SYSTEM = UTC\nSTART_TIME = {start_time}\n" \
                       "STOP_TIME = {stop_time}\nPARTICIPANT_1 = {part1}\nPARTICIPANT_2 = {sensor}\n" \
                       "MODE = SEQUENTIAL\nPATH = {part_path}\n{obs_opt}\nMETA_STOP\n\n".format(
                        start_time=obs[0]["Time"], stop_time=obs[-1]["Time"],part1=obj_id,
                        part_path=obspath, obs_opt=obstype, sensor=list(cfg["Stations"].keys())[0])

        block_ent = "DATA_START\n"
        for o in obs:
            if ("RightAscension" in miter):
                block_ent += "ANGLE_1 = {Time} {RightAscension}\n".format(**o)
                block_ent += "ANGLE_2 = {Time} {Declination}\n".format(**o)
            elif ("Azimuth" in miter):
                block_ent += "ANGLE_1 = {Time} {Azimuth}\n".format(**o)
                block_ent += "ANGLE_2 = {Time} {Elevation}\n".format(**o)
            else:
                block_ent += "RANGE = {Time} {Range}\n".format(**o)
                block_ent += "DOPPLER_INSTANTANEOUS = {Time} {RangeRate}\n".format(**o)
        tdm_data += block_header + block_ent + "DATA_STOP\n\n"

    return(tdm_data)

def import_TDM(file_name):
    return(autoclass("org.astria.Utilities").importTDM(file_name))
