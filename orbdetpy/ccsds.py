# ccsds.py - CCSDS file format conversion functions.
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

import json
import math
import grpc
from datetime import datetime
from orbdetpy.rpc import messages_pb2, utilities_pb2_grpc
from orbdetpy.rpc.server import RemoteServer
from orbdetpy.rpc.tools import convert_measurements

def format_data(cfg, obs):
    for s in cfg["cfgStations"].values():
        s["latitude"] *= 180.0/math.pi
        s["longitude"] *= 180.0/math.pi
        s["altitude"] /= 1000.0

    mit = cfg["cfgMeasurements"].keys()
    for o in obs:
        if ("station" in o):
            for i in mit:
                if (i in ["range", "rangeRate"]):
                    o[i] /= 1000.0
                else:
                    o[i] *= 180.0/math.pi

def export_OEM(cfg_file, obs_file, obj_id, obj_name):
    with open(cfg_file, "r") as fh:
        cfg = json.load(fh)
    with open(obs_file, "r") as fh:
        obs = json.load(fh)

    utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    frame = cfg.get("propInertialFrame", "EME2000")
    if (frame == "GCRF"):
        frame = "ICRF"
    cov_data = ""
    oem_data = "CCSDS_OEM_VERS = 2.0\nCREATION_DATE = {utc_time}\n" \
               "ORIGINATOR = UT-ASTRIA\n\nMETA_START\nOBJECT_NAME = {obj_name}\n" \
               "OBJECT_ID = {obj_id}\nCENTER_NAME = EARTH\n" \
               "REF_FRAME = {ref_frame}\nTIME_SYSTEM = UTC\nSTART_TIME = {start_time}\n" \
               "STOP_TIME = {stop_time}\nMETA_STOP\n\n". \
               format(utc_time=utcnow, obj_name=obj_name, obj_id=obj_id, 
                      ref_frame=frame, start_time=obs[0]["time"], stop_time=obs[-1]["time"])

    added = set()
    state_key = "estimatedState" if ("estimatedState" in obs[0]) else "trueStateCartesian"
    for o in obs:
        if (o["time"] in added):
            continue
        added.add(o["time"])

        oem_data += "{epoch} {X[0]} {X[1]} {X[2]} {X[3]} {X[4]} {X[5]}\n". \
                    format(epoch=o["time"], X=[x/1E3 for x in o[state_key][:6]])

        cov = o.get("estimatedCovariance", [])
        if (len(cov) == 0):
            cov = o.get("propagatedCovariance", [])
        if (len(cov) > 0):
            cov_data += "EPOCH = {}\n".format(o["time"])
        for i in range(1, min(len(cov), 6) + 1):
            cov_data += " ".join([str(cov[j][j]/1E6) for j in range(i)]) + "\n"

    if (len(cov_data) > 0):
        return("{}\nCOVARIANCE_START\n{}COVARIANCE_STOP\n".format(oem_data, cov_data))
    return(oem_data)

def export_TDM(cfg_file, obs_file, station_list, obj_id):
    with open(cfg_file, "r") as fh:
        cfg = json.load(fh)
    with open(obs_file, "r") as fh:
        obs = json.load(fh)
    format_data(cfg, obs)

    miter = cfg["cfgMeasurements"].keys()
    if ("rightAscension" in miter and "declination" in miter):
        frame = cfg.get("propInertialFrame", "EME2000")
        if (frame == "GCRF"):
            frame = "ICRF"
        obstype = "ANGLE_TYPE = RADEC\nREFERENCE_FRAME = {}".format(frame)
        obspath = "1,2"
    if ("azimuth" in miter and "elevation" in miter):
        obstype = "ANGLE_TYPE = AZEL"
        obspath = "1,2"
    if ("range" in miter):
        obspath = "2,1,2"
        if ("azimuth" in miter and "elevation" in miter):
            obstype = "RANGE_UNITS = km\nANGLE_TYPE = AZEL"
        else:
            obstype = "RANGE_UNITS = km"

    utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    tdm_data = "CCSDS_TDM_VERS = 1.0\nCREATION_DATE = {utc_time}\n" \
               "ORIGINATOR = UT-ASTRIA\n\n".format(utc_time=utcnow)

    for sname, sinfo in cfg["cfgStations"].items():
        if (station_list and sname not in station_list):
            continue
        sensor = "{} (Lat: {}, Lon: {}, Alt: {} km)".format(
            sname, sinfo["latitude"], sinfo["longitude"], sinfo["altitude"])
        block_header = "META_START\nTIME_SYSTEM = UTC\n" \
                       "PARTICIPANT_1 = {part1}\nPARTICIPANT_2 = {sensor}\n" \
                       "MODE = SEQUENTIAL\nPATH = {part_path}\n{obs_opt}\nMETA_STOP\n\n".format(
                           part1=obj_id, part_path=obspath, obs_opt=obstype, sensor=sensor)

        block_ent = "DATA_START\n"
        for o in obs:
            if (o.get("station") != sname):
                continue
            if ("range" in o):
                block_ent += "RANGE = {time} {range}\n".format(**o)
                if ("rangeRate" in o):
                    block_ent += "DOPPLER_INSTANTANEOUS = {time} {rangeRate}\n".format(**o)
            if ("azimuth" in o and "elevation" in o):
                block_ent += "ANGLE_1 = {time} {azimuth}\n" \
                                "ANGLE_2 = {time} {elevation}\n".format(**o)
            if ("rightAscension" in o and "declination" in o):
                block_ent += "ANGLE_1 = {time} {rightAscension}\n" \
                                "ANGLE_2 = {time} {declination}\n".format(**o)

        tdm_data += block_header + block_ent + "DATA_STOP\n\n"

    return(tdm_data)

def import_TDM(file_name, file_format):
    with RemoteServer.channel() as chan:
        stub = utilities_pb2_grpc.UtilitiesStub(chan)
        resp = stub.importTDM(messages_pb2.ImportTDMInput(
            file_name = file_name, file_format = file_format))

    marr = list(resp.array)
    return([convert_measurements(list(marr[idx].array)) for idx in range(len(marr))])
