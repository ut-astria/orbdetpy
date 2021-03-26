#!/usr/bin/env bash
# start_server.sh - Start the orbdetpy RPC server.
# Copyright (C) 2019-2021 University of Texas
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

orbd_path=$(dirname "$0")
orbd_jars=( $(ls $orbd_path/target/orbdetpy-server*.jar) )
orbd_port=$(python -c "import socket;s=socket.socket();s.bind(('127.0.0.1', 0));print(s.getsockname()[1]);s.close();")
if [ -d ${JAVA_HOME:-/dev/null}/bin ]
then
   java_exec=$JAVA_HOME/bin/java
else
   java_exec=java
fi

$java_exec -Xmx2G -XX:+UseG1GC -jar ${orbd_jars[0]} $orbd_port $orbd_path/orekit-data/
