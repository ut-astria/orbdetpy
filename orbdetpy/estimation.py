# estimation.py - Orbit estimation functions.
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

from traceback import format_exc
from orbdetpy import write_output_file
from orbdetpy.rpc import messages_pb2, estimation_pb2_grpc
from orbdetpy.rpc.server import RemoteServer
from orbdetpy.rpc.tools import (build_settings, build_measurements,
                                convert_estimation)

def determine_orbit(config, meas, output_file = None):
    """ Performs orbit determination given config and measurements.

    Args:
        config: OD configuration (Dictionary, file name, text
                file-like object, or JSON encoded string)
        meas: List of measurements (List, file name, text
                file-like object, or JSON encoded string)
        output_file: If specified, the orbit fit will be written to
                     the file name or text file-like object given

    Returns:
        Orbit determination results.
    """

    if (isinstance(config, list)):
        od_output = []
    else:
        od_output = None
        config = [config]
    if (od_output is None):
        meas = [meas]
    if (output_file and not isinstance(output_file, list)):
        output_file = [output_file]

    with RemoteServer.channel() as channel:
        stub = estimation_pb2_grpc.EstimationStub(channel)
        requests = [stub.determineOrbit.future(messages_pb2.DetermineOrbitInput(
            config = build_settings(c),
            measurements = build_measurements(m))) for c, m in zip(config, meas)]

        for idx, req in enumerate(requests):
            try:
                fit_data = convert_estimation(req.result().array)
            except Exception as exc:
                fit_data = format_exc()

            if (od_output is None):
                od_output = fit_data
            else:
                od_output.append(fit_data)
            if (output_file):
                write_output_file(output_file[idx], fit_data)

    return(od_output)

def iod_laplace(frame, lat, lon, alt, time, ra, dec):
    """ Estimate orbit from 3 RA/Dec angles using the Laplace method.

    Args:
        frame: Measurement and estimation reference frame ("EME2000", "GCRF")
        lat: Observer geodetic latitude
        lon: Observer geodetic longitude
        alt: Observer altitude
        time: List of 3 ISO-8601 formatted UTC strings [t1, t2, t3]
        ra: List of 3 Right Ascensions
        dec: List of 3 Declinations

    Returns:
        Position and velocity estimate at time t2.
    """

    with RemoteServer.channel() as channel:
        stub = estimation_pb2_grpc.EstimationStub(channel)
        resp = stub.iodLaplace(messages_pb2.AnglesInput(
            time=time, angle1=ra, angle2=dec, latitude=lat,
            longitude=lon, altitude=alt, frame=frame))

    return(list(resp.array))
