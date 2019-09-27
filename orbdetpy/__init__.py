# __init__.py - orbdetpy package initialization.
# Copyright (C) 2018-2019 University of Texas
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

import io
import os
import glob
import json
import grpc
import atexit
import psutil
from subprocess import DEVNULL
from orbdetpy.protobuf import (conversion_pb2, conversion_pb2_grpc,
                               measurements_pb2, measurements_pb2_grpc,
                               estimation_pb2, estimation_pb2_grpc)

class ServerProcess:
    rpc_host = "localhost"
    rpc_port = 50051
    rpc_uri = "{}:{}".format(rpc_host, rpc_port)

    @classmethod
    def start(cls, datadir, jarfile):
        # Check for running server instances
        jar = os.path.split(jarfile)[-1]
        for p in psutil.process_iter(attrs = ["name", "cmdline"]):
            if (p.info["name"] == "java" and
                any(x.endswith(jar) for x in p.info["cmdline"])):
                return

        # Start server
        cmdline = ["java", "-Xmx2G", "-jar", jarfile,
                   str(cls.rpc_port), datadir]
        cls.rpc_server_proc = psutil.Popen(cmdline, stdout = DEVNULL,
                                           stderr = DEVNULL)

        atexit.register(ServerProcess.stop)
        with grpc.insecure_channel(cls.rpc_uri) as chan:
            grpc.channel_ready_future(chan).result(timeout = 10.0)

    @classmethod
    def stop(cls):
        if (cls.rpc_server_proc):
            cls.rpc_server_proc.terminate()

def read_param(p):
    if (isinstance(p, str)):
        if (os.path.isfile(p)):
            with open(p, "r") as fp:
                data = fp.read()
        else:
            data = p
    elif (isinstance(p, io.TextIOBase)):
        data = p.read()
    else:
        data = json.dumps(p)

    return(data)

def write_output_file(outfile, data):
    if (isinstance(outfile, str)):
        with open(outfile, "w") as fp:
            fp.write(data)
    elif (isinstance(outfile, io.TextIOBase)):
        outfile.write(data)

def simulateMeasurements(config, output_file = None):
    """ Simulates measurement data given a configuration.

    Args:
        config: Simulation configuration (Dictionary, file name, text
                file-like object, or JSON encoded string).
        output_file: If specified, the measurements will be written to
                     the file name or text file-like object given. 

    Returns:
        Simulated measurements in the same format as config (Dictionary 
        or JSON encoded string).

    """

    with grpc.insecure_channel(ServerProcess.rpc_uri) as chan:
        stub = measurements_pb2_grpc.MeasurementsStub(chan)
        resp = stub.generateMeasurements(
            measurements_pb2.GenerateMeasurementsInput(
                config_json = read_param(config)))

    if (output_file):
        write_output_file(output_file, resp.measurements_json)
    if (isinstance(config, dict)):
        return(json.loads(resp.measurements_json))
    return(resp.measurements_json)

def determineOrbit(config, meas, output_file = None):
    """ Performs orbit determination given config and measurements.

    Args:
        config: OD configuration (Dictionary, file name, text
                file-like object, or JSON encoded string).
        meas: List of measurements (List, file name, text
                file-like object, or JSON encoded string).
        output_file: If specified, the orbit fit will be written to
                     the file name or text file-like object given. 

    Returns:
        Orbit determination results in the same format as config 
        (Dictionary or JSON encoded string).

    """

    with grpc.insecure_channel(ServerProcess.rpc_uri) as chan:
        stub = estimation_pb2_grpc.EstimationStub(chan)
        resp = stub.determineOrbit(estimation_pb2.DetermineOrbitInput(
            config_json = read_param(config),
            measurements_json = read_param(meas)))

    if (output_file):
        write_output_file(output_file, resp.estimation_json)
    if (isinstance(config, dict)):
        return(json.loads(resp.estimation_json))
    return(resp.estimation_json)

def transformFrame(srcframe, time, pva, destframe):
    """ Transforms a state vector from one frame to another.

    Args:
        srcframe: Source reference frame ("EME2000", "GCRF",
                  "ITRF", "MOD", "TOD", or "TEME").
        time: State vector epoch (ISO-8601 formatted UTC string).
        pva: State vector to transform [m, m, m, m/s, m/s, m/s] or 
             [m, m, m, m/s, m/s, m/s, m/s^2, m/s^2, m/s^2].
        destframe: Destination reference frame ("EME2000", "GCRF",
                   "ITRF", "MOD", "TOD", or "TEME")..

    Returns:
        State vector transformed to the destination frame.

    """

    with grpc.insecure_channel(ServerProcess.rpc_uri) as chan:
        stub = conversion_pb2_grpc.ConversionStub(chan)
        resp = stub.transformFrame(conversion_pb2.TransformFrameInput(
            src_frame=srcframe, time=time, pva=pva, dest_frame=destframe))
    return(list(resp.pva))

if (__name__ != '__main__'):
    _rootdir = os.path.dirname(os.path.abspath(__file__))
    _datadir = os.path.join(_rootdir, "data")
    _libsdir = os.path.join(_rootdir, "target")
    _jarfile = glob.glob(os.path.join(_libsdir, "astria*.jar"))[0]
    ServerProcess.start(_datadir, _jarfile)
