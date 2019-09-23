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
import json
import argparse
import dateutil.parser
import matplotlib.pyplot as plt
from math import acos, pi
from numpy import array, cross, dot
from numpy.linalg import norm

def main(args):
    with open(args.config, "r") as f:
        cfg = json.load(f)
    with open(args.data, "r") as f:
        out = json.load(f)

    mu = 398600.4418
    tstamp,hvec,hmag,ener,sma,ecc,inc,raan,argp,tran = [],[],[],[],[],[],[],[],[],[]
    for o in out:
        rv = [x/1000.0 for x in o["TrueState"]["Cartesian"][:6]]
        rn = norm(rv[:3])
        vn = norm(rv[3:])
        h = cross(rv[:3], rv[3:])
        hn = norm(h)
        a = 1.0/(2.0/rn-vn**2/mu)
        e = cross(rv[3:], h)/mu - rv[:3]/rn
        en = norm(e)
        i = acos(h[2]/hn)*180.0/pi
        n = cross([0, 0, 1], h)
        nn = norm(n)
        O = acos(dot(n, [1, 0, 0])/nn)*180.0/pi
        if (n[1] < 0.0):
            O = 360.0 - O
        w = acos(dot(n, e)/(nn*en))*180.0/pi
        if (e[2] < 0.0):
            w = 360.0 - w
        theta = acos(dot(rv[:3], e)/(rn*en))*180.0/pi
        if (dot(rv[:3], rv[3:]) < 0):
            theta = 360.0 - theta

        tstamp.append(dateutil.parser.isoparse(o["Time"]))
        hvec.append(h)
        hmag.append(hn)
        ener.append(0.5*vn**2 - mu/rn)
        sma.append(a)
        ecc.append(en)
        inc.append(i)
        raan.append(O)
        argp.append(w)
        tran.append(theta)

    hvec = array(hvec)
    tim = [(t - tstamp[0]).total_seconds()/3600 for t in tstamp]

    plt.figure(0)
    plt.subplot(611)
    plt.plot(tim, sma, "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("SMA [km]")
    plt.subplot(612)
    plt.plot(tim, ecc, "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("Eccentricity")
    plt.subplot(613)
    plt.plot(tim, inc, "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("Inclination [deg]")
    plt.subplot(614)
    plt.plot(tim, raan, "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("RAAN [deg]")
    plt.subplot(615)
    plt.plot(tim, argp, "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("Perigee arg. [deg]")
    plt.subplot(616)
    plt.plot(tim, tran, "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("True anom. [deg]")

    plt.figure(1)
    plt.subplot(211)
    plt.plot(tim, hmag, "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("Angular momentum")
    plt.subplot(212)
    plt.plot(tim, ener, "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("Specific energy")

    plt.figure(2)
    plt.subplot(311)
    plt.plot(tim, hvec[:,0], "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("h(x)")
    plt.subplot(312)
    plt.plot(tim, hvec[:,1], "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("h(y)")
    plt.subplot(313)
    plt.plot(tim, hvec[:,2], "ob")
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
