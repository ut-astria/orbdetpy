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

from glob import iglob
from os import environ, path, sep

_rootdir = path.dirname(path.abspath(__file__))
_libsdir = path.join(_rootdir, "lib")
_datadir = path.join(_rootdir, "data")

cpath = ""
csc = ":" if sep == "/" else ";"
jars = iglob(path.join(_libsdir, "**"), recursive = True)
for jar in jars:
    if (jar.endswith(".jar")):
        cpath += jar + csc

environ["CLASSPATH"] = cpath
from jnius import autoclass

_DataManager = autoclass("org.astria.DataManager")
_Estimation = autoclass("org.astria.Estimation")
_Simulation = autoclass("org.astria.Simulation")
_Utilities = autoclass("org.astria.Utilities")
_DataManager.initialize(_datadir)

def determineOrbit(config, meas):
    """ Performs orbit determination given config and measurements.

    Args:
        config: OD configuration as JSON encoded string.
        meas: List of measurements as JSON encoded string.

    Returns:
        Orbit determination results as JSON encoded string.

    """
    return(_Estimation(config, meas).determineOrbit())

def simulateMeasurements(config):
    """ Simulates measurement data given a configuration.

    Args:
        config: Simulation configuration as JSON encoded string.

    Returns:
        Simulated measurements as JSON encoded string.

    """
    return(_Simulation(config).simulateMeasurements())

def transformFrame(srcframe, time, pv, destframe):
    """ Transforms a state vector from one frame to another.

    Args:
        srcframe: Source reference frame.
        time: Epoch for the state vector.
        pv: State vector to transform.
        destframe: Destination reference frame.

    Returns:
        State vector transformed to the destination frame.

    """
    return(_Utilities.transformFrame(srcframe, time, pv, destframe))

def iodGooding(gslat, gslon, gsalt, tmstr, azi, ele, rho1init, rho3init):
    """ Performs Gooding initial orbit determination. """
    return(_Utilities.iodGooding(gslat, gslon, gsalt, tmstr, azi, ele,
                                 rho1init, rho3init))
