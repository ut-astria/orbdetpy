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

import io
import os
import glob
import json

_rootdir = os.path.dirname(os.path.abspath(__file__))
_libsdir = os.path.join(_rootdir, "lib")
_datadir = os.path.join(_rootdir, "data")

cpath = ""
csc = ":" if os.sep == "/" else ";"
for jar in glob.iglob(os.path.join(_libsdir, "**"), recursive = True):
    if (jar.endswith(".jar")):
        cpath += jar + csc

os.environ["CLASSPATH"] = cpath
import jnius

_DataManager = jnius.autoclass("org.astria.DataManager")
_Estimation = jnius.autoclass("org.astria.Estimation")
_Simulation = jnius.autoclass("org.astria.Simulation")
_Utilities = jnius.autoclass("org.astria.Utilities")
_DataManager.initialize(_datadir)

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

    obs = _Simulation(read_param(config)).simulateMeasurements()
    if (output_file):
        write_output_file(output_file, obs)
    if (isinstance(config, dict)):
        obs = json.loads(obs)
    return(obs)

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

    fit = _Estimation(read_param(config),
                      read_param(meas)).determineOrbit()
    if (output_file):
        write_output_file(output_file, fit)
    if (isinstance(config, dict)):
        fit = json.loads(fit)
    return(fit)

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
    return(_Utilities.transformFrame(srcframe, time, pv, destframe))

def iodGooding(gslat, gslon, gsalt, tmstr, azi, ele, rho1init, rho3init):
    """ Performs Gooding initial orbit determination. """
    return(_Utilities.iodGooding(gslat, gslon, gsalt, tmstr, azi, ele,
                                 rho1init, rho3init))
