#!/usr/bin/env bash
# generate_stubs.sh - Generate Python gRPC stubs.
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

base_dir=$(dirname "$0")
pushd $base_dir >/dev/null
python -m grpc_tools.protoc -I src/main/proto/ --python_out=./rpc --grpc_python_out=./rpc src/main/proto/*.proto
# Fix for module import path issue in gRPC generated Python code
sed -i "s/^import /import orbdetpy.rpc./g" ./rpc/*_grpc.py
sed -i "s/^import orbdetpy.rpc.grpc/import grpc/g" ./rpc/*_grpc.py
sed -i "s/^import messages_pb2/import orbdetpy.rpc.messages_pb2/g" ./rpc/*_pb2.py
popd >/dev/null
