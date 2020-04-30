# __init__.py - orbdetpy package initialization.
# Copyright (C) 2018-2020 University of Texas
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
import json
from os import path
from enum import Enum
from .version import __version__
from orbdetpy.rpc.server import RemoteServer

class Frame(Enum):
    """ Constants to identify Orekit reference frames.
    """

    CIRF_CONVENTIONS_1996_ACCURATE_EOP = "CIRF/1996 accurate EOP"
    CIRF_CONVENTIONS_1996_SIMPLE_EOP = "CIRF/1996 simple EOP"
    CIRF_CONVENTIONS_2003_ACCURATE_EOP = "CIRF/2003 accurate EOP"
    CIRF_CONVENTIONS_2003_SIMPLE_EOP = "CIRF/2003 simple EOP"
    CIRF_CONVENTIONS_2010_ACCURATE_EOP = "CIRF/2010 accurate EOP"
    CIRF_CONVENTIONS_2010_SIMPLE_EOP = "CIRF/2010 simple EOP"
    ECLIPTIC_CONVENTIONS_1996 = "Ecliptic/1996"
    ECLIPTIC_CONVENTIONS_2003 = "Ecliptic/2003"
    ECLIPTIC_CONVENTIONS_2010 = "Ecliptic/2010"
    EME2000 = "EME2000"
    GCRF = "GCRF"
    GTOD_CONVENTIONS_1996_ACCURATE_EOP = "GTOD/1996 accurate EOP"
    GTOD_CONVENTIONS_1996_SIMPLE_EOP = "GTOD/1996 simple EOP"
    GTOD_CONVENTIONS_2003_ACCURATE_EOP = "GTOD/2003 accurate EOP"
    GTOD_CONVENTIONS_2003_SIMPLE_EOP = "GTOD/2003 simple EOP"
    GTOD_CONVENTIONS_2010_ACCURATE_EOP = "GTOD/2010 accurate EOP"
    GTOD_CONVENTIONS_2010_SIMPLE_EOP = "GTOD/2010 simple EOP"
    GTOD_WITHOUT_EOP_CORRECTIONS = "GTOD/1996 without EOP"
    ICRF = "ICRF"
    ITRF_CIO_CONV_1996_ACCURATE_EOP = "CIO/1996-based ITRF accurate EOP"
    ITRF_CIO_CONV_1996_SIMPLE_EOP = "CIO/1996-based ITRF simple EOP"
    ITRF_CIO_CONV_2003_ACCURATE_EOP = "CIO/2003-based ITRF accurate EOP"
    ITRF_CIO_CONV_2003_SIMPLE_EOP = "CIO/2003-based ITRF simple EOP"
    ITRF_CIO_CONV_2010_ACCURATE_EOP = "CIO/2010-based ITRF accurate EOP"
    ITRF_CIO_CONV_2010_SIMPLE_EOP = "CIO/2010-based ITRF simple EOP"
    ITRF_EQUINOX_CONV_1996_ACCURATE_EOP = "Equinox/1996-based ITRF accurate EOP"
    ITRF_EQUINOX_CONV_1996_SIMPLE_EOP = "Equinox/1996-based ITRF simple EOP"
    ITRF_EQUINOX_CONV_2003_ACCURATE_EOP = "Equinox/2003-based ITRF accurate EOP"
    ITRF_EQUINOX_CONV_2003_SIMPLE_EOP = "Equinox/2003-based ITRF simple EOP"
    ITRF_EQUINOX_CONV_2010_ACCURATE_EOP = "Equinox/2010-based ITRF accurate EOP"
    ITRF_EQUINOX_CONV_2010_SIMPLE_EOP = "Equinox/2010-based ITRF simple EOP"
    MOD_CONVENTIONS_1996 = "MOD/1996"
    MOD_CONVENTIONS_2003 = "MOD/2003"
    MOD_CONVENTIONS_2010 = "MOD/2010"
    MOD_WITHOUT_EOP_CORRECTIONS = "MOD/1996 without EOP"
    PZ90_11 = "PZ90.11"
    TEME = "TEME"
    TIRF_CONVENTIONS_1996_ACCURATE_EOP = "TIRF/1996 accurate EOP"
    TIRF_CONVENTIONS_1996_SIMPLE_EOP = "TIRF/1996 simple EOP"
    TIRF_CONVENTIONS_2003_ACCURATE_EOP = "TIRF/2003 accurate EOP"
    TIRF_CONVENTIONS_2003_SIMPLE_EOP = "TIRF/2003 simple EOP"
    TIRF_CONVENTIONS_2010_ACCURATE_EOP = "TIRF/2010 accurate EOP"
    TIRF_CONVENTIONS_2010_SIMPLE_EOP = "TIRF/2010 simple EOP"
    TOD_CONVENTIONS_1996_ACCURATE_EOP = "TOD/1996 accurate EOP"
    TOD_CONVENTIONS_1996_SIMPLE_EOP = "TOD/1996 simple EOP"
    TOD_CONVENTIONS_2003_ACCURATE_EOP = "TOD/2003 accurate EOP"
    TOD_CONVENTIONS_2003_SIMPLE_EOP = "TOD/2003 simple EOP"
    TOD_CONVENTIONS_2010_ACCURATE_EOP = "TOD/2010 accurate EOP"
    TOD_CONVENTIONS_2010_SIMPLE_EOP = "TOD/2010 simple EOP"
    TOD_WITHOUT_EOP_CORRECTIONS = "TOD/1996 without EOP"
    VEIS_1950 = "VEIS1950"

def read_param(param):
    if (isinstance(param, str)):
        if (path.isfile(param)):
            with open(param, "r") as fp:
                data = json.load(fp)
        else:
            data = json.loads(param)
    elif (isinstance(param, io.TextIOBase)):
        data = json.load(param)
    else:
        data = param

    return(data)

def write_output_file(outfile, data):
    if (isinstance(outfile, str)):
        with open(outfile, "w") as fp:
            if (isinstance(data, str)):
                fp.write(data)
            else:
                json.dump(data, fp)
    elif (isinstance(outfile, io.TextIOBase)):
        if (isinstance(data, str)):
            outfile.write(data)
        else:
            json.dump(data, outfile)

if (__name__ != '__main__'):
    _rootdir = path.dirname(path.abspath(__file__))
    _datadir = path.join(_rootdir, "data")
    _libsdir = path.join(_rootdir, "target")
    _jarfile = path.join(_libsdir, "orbdetpy-server-{}.jar".format(__version__))
    RemoteServer.connect(_datadir, _jarfile)
