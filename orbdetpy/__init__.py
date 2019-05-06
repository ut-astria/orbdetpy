# __init__.py - Package initialization routines.
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

from os import environ, path, walk, sep

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

from jnius import autoclass

_DataManager = autoclass("org.astria.DataManager")
_Estimation = autoclass("org.astria.Estimation")
_Simulation = autoclass("org.astria.Simulation")
_Utilities = autoclass("org.astria.Utilities")

_DataManager.initialize(_datadir)

def determineOrbit(config, meas):
    return(_Estimation(config, meas).determineOrbit())

def simulateMeasurements(config):
    return(_Simulation(config).simulateMeasurements())

def transformFrame(srcframe, time, pv, destframe):
    return(_Utilities.transformFrame(srcframe, time, pv, destframe))

def iodGooding(gslat, gslon, gsalt, tmstr, azi, ele, rho1init, rho3init):
    return(_Utilities.iodGooding(gslat, gslon, gsalt, tmstr, azi, ele,
                                 rho1init, rho3init))
