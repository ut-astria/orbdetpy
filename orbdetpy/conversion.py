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

from orbdetpy.rpc import messages_pb2, conversion_pb2_grpc
from orbdetpy.rpc.server import RemoteServer

def transform_frame(srcframe, time, pva, destframe):
    """ Transforms a state vector from one frame to another.

    Args:
        srcframe: Source reference frame; a constant from 
                  the enum Frame
        time: State vector epoch (ISO-8601 formatted UTC string)
        pva: State vector to transform, can be pos or pos + vel or
             pos + vel + acc
        destframe: Destination reference frame; a constant from 
                   the enum Frame

    Returns:
        State vector transformed to the destination frame.
    """

    with RemoteServer.channel() as channel:
        stub = conversion_pb2_grpc.ConversionStub(channel)
        resp = stub.transformFrame(messages_pb2.TransformFrameInput(
            src_frame=srcframe.name, time=time,
            pva=pva, dest_frame=destframe.name))

    return(list(resp.array))

def azel_to_radec(time, az, el, lat, lon, alt, frame):
    """ Convert Azimuth/Elevation to Right Ascension/Declination.

    Args:
        time: ISO-8601 formatted UTC string
        az: Azimuth
        el: Elevation
        lat: Observer geodetic latitude
        lon: Observer geodetic longitude
        alt: Observer altitude
        frame: Destination reference frame; Frame.EME2000 or Frame.GCRF

    Returns:
        Right Ascension and Declination.
    """

    with RemoteServer.channel() as channel:
        stub = conversion_pb2_grpc.ConversionStub(channel)
        resp = stub.convertAzElToRaDec(messages_pb2.AnglesInput(
            time=[time], angle1=[az], angle2=[el], latitude=lat,
            longitude=lon, altitude=alt, frame=frame.name))

    return(list(resp.array))

def radec_to_azel(frame, time, ra, dec, lat, lon, alt):
    """ Convert Right Ascension/Declination to Azimuth/Elevation.

    Args:
        frame: Source reference frame; Frame.EME2000 or Frame.GCRF
        time: ISO-8601 formatted UTC string
        ra: Right Ascension
        dec: Declination
        lat: Observer geodetic latitude
        lon: Observer geodetic longitude
        alt: Observer altitude

    Returns:
        Azimuth and Elevation.
    """

    with RemoteServer.channel() as channel:
        stub = conversion_pb2_grpc.ConversionStub(channel)
        resp = stub.convertRaDecToAzEl(messages_pb2.AnglesInput(
            time=[time], angle1=[ra], angle2=[dec], latitude=lat,
            longitude=lon, altitude=alt, frame=frame.name))

    return(list(resp.array))
