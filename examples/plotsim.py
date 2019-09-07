# plotsim.py - Plot simulated data.
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

import sys
import math
import json
import numpy
import argparse
from numpy.linalg import norm
import dateutil.parser
import matplotlib.pyplot as plt

def main(args):
    with open(args.config, "r") as f:
        cfg = json.load(f)
    with open(args.data, "r") as f:
        out = json.load(f)

    mu = 398600.4418
    tstamp, hvec, hmag, ener, sma, ecc, inc, raan, argp = [], [], [], [], [], [], [], [], []
    for o in out:
        rv = [x/1000.0 for x in o["TrueState"]["Cartesian"][:6]]
        r = norm(rv[:3])
        v = norm(rv[3:])
        h = numpy.cross(rv[:3], rv[3:])
        hn = norm(h)
        a = 1.0/(2.0/r-v*v/mu)
        e = numpy.cross(rv[3:], h)/mu - rv[:3]/r
        en = norm(e)
        i = math.acos(h[2]/hn)*180.0/math.pi
        n = numpy.cross([0, 0, 1], h)
        nn = norm(n)
        O = math.acos(numpy.dot(n, [1, 0, 0])/nn)*180.0/math.pi
        if (numpy.dot([0, 1, 0], n) < 0.0):
            O = 360.0 - O
        w = math.acos(numpy.dot(n, e)/(nn*en))*180.0/math.pi
        if (numpy.dot([0, 0, 1], e) < 0.0):
            w = 360.0 - w

        tstamp.append(dateutil.parser.isoparse(o["Time"]))
        hvec.append(h)
        hmag.append(hn)
        ener.append(v**2/2.0 - mu/r)
        sma.append(a)
        ecc.append(en)
        inc.append(i)
        raan.append(O)
        argp.append(w)

    hvec = numpy.array(hvec)
    tim = [(t - tstamp[0]).total_seconds()/3600 for t in tstamp]

    plt.figure(0)
    plt.subplot(511)
    plt.plot(tim, sma, "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("SMA [km]")
    plt.subplot(512)
    plt.plot(tim, ecc, "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("Eccentricity")
    plt.subplot(513)
    plt.plot(tim, inc, "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("Inclination [deg]")
    plt.subplot(514)
    plt.plot(tim, raan, "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("RAAN [deg]")
    plt.subplot(515)
    plt.plot(tim, argp, "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("Perigee arg. [deg]")

    plt.figure(1)
    plt.subplot(211)
    plt.plot(tim, hmag, "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("Angular momentum")
    plt.subplot(212)
    plt.plot(tim, ener, "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("Specific energy")

    plt.figure(2)
    plt.subplot(311)
    plt.plot(tim, hvec[:,0], "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("h(x)")
    plt.subplot(312)
    plt.plot(tim, hvec[:,1], "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("h(y)")
    plt.subplot(313)
    plt.plot(tim, hvec[:,2], "-b")
    plt.xlabel("Time [hr]")
    plt.ylabel("h(z)")
    plt.suptitle("Components of angular momentum")
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description='Plotting utility for sim data.')
    parser.add_argument('config', help='Path to config file.')
    parser.add_argument('data', help='Path to simulated data file.')
    args = parser.parse_args()
    main(args)
