#!/bin/bash
# start_RPC_server.sh - Start the orbdetpy RPC server.
# Copyright (C) 2019 University of Texas
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

BASE_DIR=$(dirname "$0")

for jar_file in $BASE_DIR/target/orbdetpy-server*.jar
do
    java -Xmx2G -XX:+UseG1GC -jar $jar_file 50051 data/
done
