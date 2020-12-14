# conversion.py - Various conversion functions.
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

from numpy import array, cross, vstack
from numpy.linalg import norm
from typing import List, Tuple
from google.protobuf.wrappers_pb2 import StringValue
from orbdetpy.rpc.conversion_pb2_grpc import ConversionStub
from orbdetpy.rpc.messages_pb2 import AnglesInput, DoubleArray, IntegerArray, TransformFrameInput
from orbdetpy.rpc.server import RemoteServer

def transform_frame(src_frame: int, time: float, pva: List[float], dest_frame: int)->List[float]:
    """Transform a state vector from one frame to another.

    Parameters
    ----------
    src_frame : Source reference frame; a constant from Frame.
    time : Offset in TT from J2000 epoch [s]. Give a list for bulk transforms.
    pva : State vector to transform, can be pos or pos+vel or pos+vel+acc. Provide a 
          list of lists for bulk frame transforms.
    dest_frame : Destination reference frame; a constant from Frame.

    Returns
    -------
    State vector transformed to the destination frame.
    """

    if (isinstance(time, float) or isinstance(time, str)):
        single = True
        time, pva = [time], [pva]
    else:
        single = False

    if (isinstance(time[0], float)):
        resp = _conversion_stub.transformFrame(TransformFrameInput(
            src_frame=src_frame, time=time, pva=[DoubleArray(array=x) for x in pva], dest_frame=dest_frame))
    else:
        resp = _conversion_stub.transformFrame(TransformFrameInput(
            src_frame=src_frame, UTC_time=time, pva=[DoubleArray(array=x) for x in pva], dest_frame=dest_frame))
    return(resp.array[0].array if (single) else resp.array)

def get_lvlh_rotation(state: List[float])->array:
    """Construct the inertial->LVLH rotation matrix.

    Parameters
    ----------
    state : Inertial state vector.

    Returns
    -------
    Inertial->LVLH frame rotation matrix.
    """

    r = array(state[:3])
    r /= norm(r)
    v = array(state[3:6])
    v /= norm(v)
    h = cross(r, v)
    h /= norm(h)
    return(vstack((r, cross(h, r), h)))

def azel_to_radec(time: float, az: float, el: float, lat: float, lon: float, alt: float, frame: int)->Tuple[float, float]:
    """Convert Azimuth/Elevation to Right Ascension/Declination.

    Parameters
    ----------
    time : Offset in TT from J2000 epoch [s].
    az : Azimuth [rad].
    el : Elevation [rad].
    lat : Observer WGS-84 latitude [rad].
    lon : Observer WGS-84 longitude [rad].
    alt : Observer height above WGS-84 reference ellipsoid [m].
    frame : Destination reference frame; Frame.GCRF or Frame.EME2000.

    Returns
    -------
    Right Ascension and Declination.
    """

    resp = _conversion_stub.convertAzElToRaDec(AnglesInput(time=[time], angle1=[az], angle2=[el],
                                                           latitude=lat, longitude=lon, altitude=alt, frame=frame))
    return(resp.array)

def radec_to_azel(frame: int, time: float, ra: float, dec: float, lat: float, lon: float, alt: float)->Tuple[float, float]:
    """Convert Right Ascension/Declination to Azimuth/Elevation.

    Parameters
    ----------
    frame : Source reference frame; Frame.GCRF or Frame.EME2000.
    time : Offset in TT from J2000 epoch [s].
    ra : Right Ascension [rad].
    dec : Declination [rad].
    lat : Observer WGS-84 latitude [rad].
    lon : Observer WGS-84 longitude [rad].
    alt : Observer height above WGS-84 reference ellipsoid [m].

    Returns
    -------
    Azimuth and Elevation.
    """

    resp = _conversion_stub.convertRaDecToAzEl(AnglesInput(time=[time], angle1=[ra], angle2=[dec],
                                                           latitude=lat, longitude=lon, altitude=alt, frame=frame))
    return(resp.array)

def pos_to_lla(frame: int, time: float, pos: List[float])->List[float]:
    """Convert an inertial position to WGS-84 lat/lon/alt.

    Parameters
    ----------
    frame : Inertial reference frame; a constant from Frame.
    time : Offset in TT from J2000 epoch [s]. Give a list for bulk conversions.
    pos : Inertial position vector. Provide a list of lists for bulk conversions.

    Returns
    -------
    WGS-84 latitude [rad], longitude [rad], altitude [m].
    """

    if (isinstance(time, float) or isinstance(time, str)):
        single = True
        time, pos = [time], [pos]
    else:
        single = False

    if (isinstance(time[0], float)):
        resp = _conversion_stub.convertPosToLLA(TransformFrameInput(
            src_frame=frame, time=time, pva=[DoubleArray(array=x) for x in pos]))
    else:
        resp = _conversion_stub.convertPosToLLA(TransformFrameInput(
            src_frame=frame, UTC_time=time, pva=[DoubleArray(array=x) for x in pos]))
    return(resp.array[0].array if (single) else resp.array)

def elem_to_pv(frame: int, time: float, sma: float, ecc: float, inc: float,
               raan: float, argp: float, anom: float, anom_type: int)->List[float]:
    """Convert osculating orbital elements to Cartesian state vector.

    Parameters
    ----------
    frame : Inertial reference frame; a constant from Frame.
    time : Offset in TT from J2000 epoch [s]. Give a list for bulk conversions.
    sma : Semi-major axis [m]. Give a list for bulk conversions.
    ecc : Eccentricity. Give a list for bulk conversions.
    inc : Inclination [rad]. Give a list for bulk conversions.
    raan : RA of ascending node [rad]. Give a list for bulk conversions.
    argp : Argument of perigee [rad]. Give a list for bulk conversions.
    anom : Anomaly angle [rad]. Give a list for bulk conversions.
    anom_type : Anomaly angle type; a constant from PositionAngle. Give a list for bulk conversions.

    Returns
    -------
    Cartesian state vector.
    """

    if (isinstance(time, float) or isinstance(time, str)):
        single = True
        time, sma, ecc, inc, raan, argp, anom, anom_type = (
            [time], [sma], [ecc], [inc], [raan], [argp], [anom], [anom_type])
    else:
        single = False

    if (isinstance(time[0], float)):
        resp = _conversion_stub.convertElemToPv(TransformFrameInput(src_frame=frame, time=time, pva=[
            DoubleArray(array=[a,e,i,W,w,n,t]) for a,e,i,W,w,n,t in zip(sma, ecc, inc, raan, argp, anom, anom_type)]))
    else:
        resp = _conversion_stub.convertElemToPv(TransformFrameInput(src_frame=frame, UTC_time=time, pva=[
            DoubleArray(array=[a,e,i,W,w,n,t]) for a,e,i,W,w,n,t in zip(sma, ecc, inc, raan, argp, anom, anom_type)]))
    return(resp.array[0].array if (single) else resp.array)

def pv_to_elem(frame: int, time: float, pv: List[float])->List[float]:
    """Convert Cartesian state vector to osculating orbital elements.

    Parameters
    ----------
    frame : Inertial reference frame; a constant from Frame.
    time : Offset in TT from J2000 epoch [s]. Give a list for bulk conversions.
    pv : Inertial Cartesian state vector. Provide a list of lists for bulk conversions.

    Returns
    -------
    SMA, eccentricity, inclination, RAAN, perigee argument, mean anomaly, true anomaly, eccentric anomaly
    """

    if (isinstance(time, float) or isinstance(time, str)):
        single = True
        time, pv = [time], [pv]
    else:
        single = False

    if (isinstance(time[0], float)):
        resp = _conversion_stub.convertPvToElem(TransformFrameInput(
            src_frame=frame, time=time, pva=[DoubleArray(array=x) for x in pv]))
    else:
        resp = _conversion_stub.convertPvToElem(TransformFrameInput(
            src_frame=frame, UTC_time=time, pva=[DoubleArray(array=x) for x in pv]))
    return(resp.array[0].array if (single) else resp.array)

def get_UTC_string(j2000_offset: float, truncate: bool=True)->str:
    """Get ISO-8601 formatted UTC string corresponding to TT offset.

    Parameters
    ----------
    j2000_offset : Offset in TT from J2000 epoch [s] or list of offsets.
    truncate : Truncate to milliseconds level accuracy if True (default).

    Returns
    -------
    ISO-8601 formatted UTC string or list of strings.
    """

    if (isinstance(j2000_offset, float)):
        return(_conversion_stub.getUTCString(DoubleArray(array=[float(truncate), j2000_offset])).value)
    return(_conversion_stub.getUTCString(DoubleArray(array=[float(truncate), *j2000_offset])).value.split(" "))

def get_J2000_epoch_offset(utc: str)->float:
    """Get TT offset corresponding to ISO-8601 formatted UTC string.

    Parameters
    ----------
    utc : ISO-8601 formatted UTC string or list of strings.

    Returns
    -------
    Offset in TT from J2000 epoch [s] or list of offsets.
    """

    resp = _conversion_stub.getJ2000EpochOffset(StringValue(value=utc if isinstance(utc, str) else " ".join(utc)))
    return(resp.array[0] if (len(resp.array) == 1) else resp.array)

def get_epoch_difference(from_epoch: int, to_epoch: int)->str:
    """Get the constant time difference between two epochs.

    Parameters
    ----------
    from_epoch : From epoch; a value from the Epoch enumeration.
    to_epoch : To epoch; a value from the Epoch enumeration.

    Returns
    -------
    to_epoch - from_epoch in seconds.
    """

    return(_conversion_stub.getEpochDifference(IntegerArray(array=[from_epoch, to_epoch])).value)

if (__name__ != '__main__'):
    __pdoc__ = {m: False for m in ("AnglesInput", "DoubleArray", "IntegerArray", "TransformFrameInput")}
    _conversion_stub = ConversionStub(RemoteServer.channel())
