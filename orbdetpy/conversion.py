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

from typing import List, Tuple
from google.protobuf.wrappers_pb2 import DoubleValue, StringValue
from orbdetpy.rpc.conversion_pb2_grpc import ConversionStub
from orbdetpy.rpc.messages_pb2 import AnglesInput, TransformFrameInput
from orbdetpy.rpc.server import RemoteServer

def transform_frame(src_frame: int, time: float, pva: List[float], dest_frame: int)->List[float]:
    """Transform a state vector from one frame to another.

    Parameters
    ----------
    src_frame : Source reference frame; a constant from Frame.
    time : Offset in TT from J2000 epoch [s].
    pva : State vector to transform, can be pos or pos+vel or pos+vel+acc.
    dest_frame : Destination reference frame; a constant from Frame.

    Returns
    -------
    State vector transformed to the destination frame.
    """

    with RemoteServer.channel() as chan:
        resp = ConversionStub(chan).transformFrame(TransformFrameInput(
            src_frame=src_frame, time=time, pva=pva, dest_frame=dest_frame))
    return(resp.array)

def azel_to_radec(time: float, az: float, el: float, lat: float,
                  lon: float, alt: float, frame: int)->Tuple[float, float]:
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

    with RemoteServer.channel() as chan:
        resp = ConversionStub(chan).convertAzElToRaDec(AnglesInput(
            time=[time], angle1=[az], angle2=[el], latitude=lat,
            longitude=lon, altitude=alt, frame=frame))
    return(resp.array)

def radec_to_azel(frame: int, time: float, ra: float, dec: float,
                  lat: float, lon: float, alt: float)->Tuple[float, float]:
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

    with RemoteServer.channel() as chan:
        resp = ConversionStub(chan).convertRaDecToAzEl(AnglesInput(
            time=[time], angle1=[ra], angle2=[dec], latitude=lat,
            longitude=lon, altitude=alt, frame=frame))
    return(resp.array)

def pos_to_lla(frame: int, time: float, pos: List[float])->List[float]:
    """Convert an inertial position to WGS-84 lat/lon/alt.

    Parameters
    ----------
    frame : Inertial reference frame; a constant from Frame.
    time : Offset in TT from J2000 epoch [s].
    pos : Inertial position vector.

    Returns
    -------
    WGS-84 latitude [rad], longitude [rad], altitude [m].
    """

    with RemoteServer.channel() as chan:
        resp = ConversionStub(chan).convertPosToLLA(TransformFrameInput(
            src_frame=frame, time=time, pva=pos))
    return(resp.array)

def elem_to_pv(frame: int, time: float, sma: float, ecc: float, inc: float,
               raan: float, argp: float, mean_anom: float)->List[float]:
    """Convert osculating orbital elements to Cartesian state vector.

    Parameters
    ----------
    frame : Inertial reference frame; a constant from Frame.
    time : Offset in TT from J2000 epoch [s].
    sma : Semi-major axis [m].
    ecc : Eccentricity.
    inc : Inclination [rad].
    raan : RA of ascending node [rad].
    argp : Argument of perigee [rad].
    mean_anom : Mean anomaly [rad].

    Returns
    -------
    Cartesian state vector.
    """

    with RemoteServer.channel() as chan:
        resp = ConversionStub(chan).convertElemToPv(TransformFrameInput(
            src_frame=frame, time=time, pva=[sma,ecc,inc,raan,argp,mean_anom]))
    return(resp.array)

def pv_to_elem(frame: int, time: float, pv: List[float])->List[float]:
    """Convert Cartesian state vector to osculating orbital elements.

    Parameters
    ----------
    frame : Inertial reference frame; a constant from Frame.
    time : Offset in TT from J2000 epoch [s].
    pv : Inertial Cartesian state vector.

    Returns
    -------
    SMA,eccentricity,inclination,RAAN,perigee argument,mean anomaly,true anomaly
    """

    with RemoteServer.channel() as chan:
        resp = ConversionStub(chan).convertPvToElem(TransformFrameInput(
            src_frame=frame, time=time, pva=pv))
    return(resp.array)

def get_UTC_string(j2000_offset: float)->str:
    """Get ISO-8601 formatted UTC string corresponding to TT offset.

    Parameters
    ----------
    j2000_offset : Offset in TT from J2000 epoch [s].

    Returns
    -------
    ISO-8601 formatted UTC string.
    """

    with RemoteServer.channel() as chan:
        resp = ConversionStub(chan).getUTCString(DoubleValue(value=j2000_offset))
    return(resp.value)

def get_J2000_epoch_offset(utc_time: str)->float:
    """Get TT offset corresponding to ISO-8601 formatted UTC string.

    Parameters
    ----------
    utc_time : ISO-8601 formatted UTC string.

    Returns
    -------
    Offset in TT from J2000 epoch [s].
    """

    with RemoteServer.channel() as chan:
        resp = ConversionStub(chan).getJ2000EpochOffset(StringValue(value=utc_time))
    return(resp.value)
