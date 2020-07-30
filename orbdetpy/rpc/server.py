# server.py - RPC server connection management.
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

import os
import atexit
import psutil
from grpc import channel_ready_future, insecure_channel

class RemoteServer:
    rpc_host = "localhost"
    rpc_port = "50051"
    rpc_uri = f"{rpc_host}:{rpc_port}"
    rpc_server = None
    rpc_channel = None

    @classmethod
    def connect(cls, data_dir, jar_file):
        running = False
        jar = os.path.split(jar_file)[-1]
        for p in psutil.process_iter(attrs=["name", "cmdline"]):
            running = p.info["name"] == "java" and any(x.endswith(jar) for x in p.info["cmdline"])
            if (running):
                break

        if (not running):
            cls.rpc_server = psutil.Popen(["java", "-Xmx2G", "-XX:+UseG1GC", "-jar", jar_file, cls.rpc_port, data_dir])
        atexit.register(RemoteServer.disconnect)
        cls.rpc_channel = insecure_channel(cls.rpc_uri, options=[("grpc.max_send_message_length", 2147483647),
                                                                 ("grpc.max_receive_message_length", 2147483647)])
        channel_ready_future(cls.rpc_channel).result(timeout=30.0)

    @classmethod
    def channel(cls):
        return(cls.rpc_channel)

    @classmethod
    def disconnect(cls):
        try:
            if (cls.rpc_channel is not None):
                cls.rpc_channel.close()
            if (cls.rpc_server is not None):
                cls.rpc_server.terminate()
        except Exception as exc:
            pass
