#!/bin/bash
# build.sh - Script to build Java source code.
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

pushd $(dirname "$0") > /dev/null

export CLASSPATH=`find ../lib -name "*.jar" | tr "\n" ":"`

find . -name "*.class" | xargs rm
find . -name "*.java" | xargs javac -Xlint:unchecked -Xlint:deprecation

CLASS=`find . -name "*.class" | tr "\n" " "`
jar cf ../lib/astria.jar $CLASS

popd > /dev/null
