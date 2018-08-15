# __init__.py - Package initialization routines.
# Copyright (C) 2018 Shiva Iyer <shiva.iyer AT utexas DOT edu>
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

from os import environ, path, walk, sep

def init():
    global _rootdir, _libsdir, _datadir

    _rootdir = path.dirname(path.dirname(path.abspath(__file__)))
    _libsdir = path.join(_rootdir, "lib")
    _datadir = path.join(_rootdir, "data")

    cpath = ""
    csc = ":" if sep == "/" else ";"
    for r, d, files in walk(_libsdir):
        for file in files:
            if (file.endswith(".jar")):
                cpath += path.join(r, file) + csc
    environ["CLASSPATH"] = cpath
