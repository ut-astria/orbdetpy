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

from datetime import datetime
from typing import List, Optional
from orbdetpy import Constant, Frame, MeasurementType
from orbdetpy.conversion import get_UTC_string
from orbdetpy.rpc.messages_pb2 import ImportTDMInput, Settings
from orbdetpy.rpc.server import RemoteServer
from orbdetpy.rpc.utilities_pb2_grpc import UtilitiesStub

_ccsds_stub = UtilitiesStub(RemoteServer.channel())

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
    oem_data = f"""CCSDS_OEM_VERS = 2.0
CREATION_DATE = {datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")}
ORIGINATOR = UTexas-Austin

META_START
OBJECT_NAME = {obj_name}
OBJECT_ID = {obj_id}
CENTER_NAME = EARTH
REF_FRAME = {frame}
TIME_SYSTEM = UTC
START_TIME = {get_UTC_string(obs[0].time)}
STOP_TIME = {get_UTC_string(obs[-1].time)}
META_STOP
"""

    cov_data, added  = "", set()
    state_key = "estimated_state" if (hasattr(obs[0], "estimated_state")) else "true_state"
    for o in obs:
        if (o.time in added):
            continue
        added.add(o.time)
        utc = get_UTC_string(o.time)
        X=[x/1E3 for x in getattr(o, state_key)[:6]]
        oem_data = f"{oem_data}\n{utc} {X[0]} {X[1]} {X[2]} {X[3]} {X[4]} {X[5]}"

        cov = []
        if (hasattr(o, "estimated_covariance")):
            cov = o.estimated_covariance
        elif (hasattr(o, "propagated_covariance")):
            cov = o.propagated_covariance
        if (len(cov) >= 21):
            cov_data += f"EPOCH = {utc}\n"
            cov_data += " ".join([str(c/1E6) for c in cov[:21]]) + "\n"

    if (len(cov_data) > 0):
        return(f"{oem_data}\n\nCOVARIANCE_START\n{cov_data}COVARIANCE_STOP")
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
        obstype = f"ANGLE_TYPE = RADEC\nREFERENCE_FRAME = {frame}"
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

    tdm_data = f"""CCSDS_TDM_VERS = 1.0
CREATION_DATE = {datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%SZ")}
ORIGINATOR = UTexas-Austin
"""
    for sname, sinfo in cfg.stations.items():
        if (station_list is not None and sname not in station_list):
            continue
        sensor = f"{sname} (Lat:{sinfo.latitude/Constant.DEGREE_TO_RAD},Lon:{sinfo.longitude/Constant.DEGREE_TO_RAD},Alt:{sinfo.altitude/1E3}km)"
        block_header = f"""
META_START
TIME_SYSTEM = UTC
PARTICIPANT_1 = {obj_id}
PARTICIPANT_2 = {sensor}
MODE = SEQUENTIAL
PATH = {obspath}
{obstype}
META_STOP

"""
        block_ent = "DATA_START\n"
        for o in obs:
            if (o.station != sname):
                continue
            utc = get_UTC_string(o.time)
            if (MeasurementType.RANGE in miter):
                block_ent = f"{block_ent}RANGE = {utc} {o.values[0]/1E3}\n"
                if ("rangeRate" in o):
                    block_ent = f"{block_ent}DOPPLER_INSTANTANEOUS = {utc} {o.values[1]/1E3}\n"
            if ((MeasurementType.AZIMUTH in miter and MeasurementType.ELEVATION in miter) or
                (MeasurementType.RIGHT_ASCENSION in miter and MeasurementType.DECLINATION in miter)):
                block_ent = f"""{block_ent}ANGLE_1 = {utc} {o.values[0]/Constant.DEGREE_TO_RAD}
ANGLE_2 = {utc} {o.values[1]/Constant.DEGREE_TO_RAD}
"""
        tdm_data = f"{tdm_data}{block_header}{block_ent}DATA_STOP\n"

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

    resp = _ccsds_stub.importTDM(ImportTDMInput(file_name=file_name, file_format=file_format))
    return(resp.array)
