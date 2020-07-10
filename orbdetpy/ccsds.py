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

import math
from datetime import datetime
from typing import List, Optional
from orbdetpy import Frame, MeasurementType
from orbdetpy.conversion import get_UTC_string
from orbdetpy.rpc.messages_pb2 import ImportTDMInput, Settings
from orbdetpy.rpc.server import RemoteServer
from orbdetpy.rpc.utilities_pb2_grpc import UtilitiesStub

def export_OEM(cfg: Settings, obs, obj_id: str, obj_name: str)->str:
    """Export ephemerides in CCSDS OEM format.

    Parameters
    ----------
    cfg : Settings object.
    obs : Measurements or estimation results to export.
    obj_id : Object identifier.
    obj_name : Object name.

    Returns
    -------
    Ephemerides in OEM format.
    """

    frame = Frame.ICRF if (cfg.prop_inertial_frame == Frame.GCRF) else cfg.prop_inertial_frame
    utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    oem_data = "CCSDS_OEM_VERS = 2.0\nCREATION_DATE = {utc_time}\n" \
               "ORIGINATOR = UT-ASTRIA\n\nMETA_START\nOBJECT_NAME = {obj_name}\n" \
               "OBJECT_ID = {obj_id}\nCENTER_NAME = EARTH\n" \
               "REF_FRAME = {ref_frame}\nTIME_SYSTEM = UTC\nSTART_TIME = {start_time}\n" \
               "STOP_TIME = {stop_time}\nMETA_STOP\n\n". \
               format(utc_time=utcnow, obj_name=obj_name, obj_id=obj_id, ref_frame=frame,
                      start_time=get_UTC_string(obs[0].time), stop_time=get_UTC_string(obs[-1].time))

    cov_data, added  = "", set()
    state_key = "estimated_state" if (hasattr(obs[0], "estimated_state")) else "true_state"
    for o in obs:
        if (o.time in added):
            continue
        added.add(o.time)

        utc = get_UTC_string(o.time)
        oem_data += ("{epoch} {X[0]} {X[1]} {X[2]} {X[3]} {X[4]} {X[5]}\n".
                     format(epoch=utc, X=[x/1000.0 for x in getattr(o, state_key)[:6]]))

        cov = []
        if (hasattr(o, "estimated_covariance")):
            cov = o.estimated_covariance
        elif (hasattr(o, "propagated_covariance")):
            cov = o.propagated_covariance
        if (len(cov) >= 21):
            cov_data += "EPOCH = {}\n".format(utc)
            cov_data += " ".join([str(c/1E6) for c in cov[:21]]) + "\n"

    if (len(cov_data) > 0):
        return("{}\nCOVARIANCE_START\n{}COVARIANCE_STOP\n".format(oem_data, cov_data))
    return(oem_data)

def export_TDM(cfg: Settings, obs, obj_id: str, station_list: Optional[List[str]]=None)->str:
    """Export tracking data in CCSDS TDM format.

    Parameters
    ----------
    cfg : Settings object.
    obs : Measurements to export.
    obj_id : Object identifier.
    station_list : List of ground stations to include; None to include all.

    Returns
    -------
    Tracking data in TDM format.
    """

    miter = cfg.measurements.keys()
    if (MeasurementType.RIGHT_ASCENSION in miter and MeasurementType.DECLINATION in miter):
        frame = Frame.ICRF if (cfg.prop_inertial_frame == Frame.GCRF) else cfg.prop_inertial_frame
        obstype = "ANGLE_TYPE = RADEC\nREFERENCE_FRAME = {}".format(frame)
        obspath = "1,2"
    if (MeasurementType.AZIMUTH in miter and MeasurementType.ELEVATION in miter):
        obstype = "ANGLE_TYPE = AZEL"
        obspath = "1,2"
    if (MeasurementType.RANGE in miter):
        obspath = "2,1,2"
        if (MeasurementType.AZIMUTH in miter and MeasurementType.ELEVATION in miter):
            obstype = "RANGE_UNITS = km\nANGLE_TYPE = AZEL"
        else:
            obstype = "RANGE_UNITS = km"

    utcnow = datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")
    tdm_data = "CCSDS_TDM_VERS = 1.0\nCREATION_DATE = {utc_time}\n" \
               "ORIGINATOR = UT-ASTRIA\n\n".format(utc_time=utcnow)

    for sname, sinfo in cfg.stations.items():
        if (station_list is not None and sname not in station_list):
            continue
        sensor = "{} (Lat: {}, Lon: {}, Alt: {} km)".format(
            sname, sinfo.latitude*180.0/math.pi, sinfo.longitude*180.0/math.pi, sinfo.altitude/1E3)
        block_header = "META_START\nTIME_SYSTEM = UTC\n" \
                       "PARTICIPANT_1 = {part1}\nPARTICIPANT_2 = {sensor}\n" \
                       "MODE = SEQUENTIAL\nPATH = {part_path}\n{obs_opt}\nMETA_STOP\n\n".format(
                           part1=obj_id, part_path=obspath, obs_opt=obstype, sensor=sensor)

        block_ent = "DATA_START\n"
        for o in obs:
            if (o.station != sname):
                continue
            utc = get_UTC_string(o.time)
            if (MeasurementType.RANGE in miter):
                block_ent += "RANGE = {} {}\n".format(utc, o.values[0]/1E3)
                if ("rangeRate" in o):
                    block_ent += "DOPPLER_INSTANTANEOUS = {} {}\n".format(utc, o.values[1]/1E3)
            if ((MeasurementType.AZIMUTH in miter and MeasurementType.ELEVATION in miter) or
                (MeasurementType.RIGHT_ASCENSION in miter and MeasurementType.DECLINATION in miter)):
                block_ent += "ANGLE_1 = {} {}\nANGLE_2 = {} {}\n".format(
                    utc, o.values[0]*180.0/math.pi, utc, o.values[1]*180.0/math.pi)
        tdm_data += block_header + block_ent + "DATA_STOP\n\n"

    return(tdm_data)

def import_TDM(file_name: str, file_format: int):
    """Import tracking data from CCSDS TDM file.

    Parameters
    ----------
    file_name : Fully qualified TDM file.
    file_format : Constant from TDMFileFormat.

    Returns
    -------
    Measurements object.
    """

    with RemoteServer.channel() as chan:
        resp = UtilitiesStub(chan).importTDM(ImportTDMInput(
            file_name=file_name, file_format=file_format))
    return(resp.array)
