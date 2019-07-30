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
from numpy.linalg import norm
import dateutil.parser
import matplotlib.pyplot as plt

if (len(sys.argv) < 3):
    print("Usage: python %s config_file simulated_data_file" % sys.argv[0])
    exit()

with open(sys.argv[1], "r") as f:
    cfg = json.load(f)
with open(sys.argv[2], "r") as f:
    out = json.load(f)

mu = 398600.4418
tstamp, hvec, hmag, ener, sma, ecc, inc, added = [], [], [], [], [], [], [], []
for o in out:
    if (o["Time"] in added):
        continue
    added.append(o["Time"])

    rv = [x/1000.0 for x in o["TrueState"]["Cartesian"][:6]]
    r, v = norm(rv[:3]), norm(rv[3:])
    h = numpy.cross(rv[:3], rv[3:])
    hn = norm(h)
    a = 1.0/(2.0/r-v*v/mu)
    e = math.sqrt(1 - hn*hn/(mu*a))
    i = math.acos(h[2]/hn)*180.0/math.pi

    tstamp.append(dateutil.parser.isoparse(o["Time"]))
    hvec.append(h)
    hmag.append(hn)
    ener.append(v**2/2.0 - mu/r)
    sma.append(a)
    ecc.append(e)
    inc.append(i)

hvec = numpy.array(hvec)
tim = [(t - tstamp[0]).total_seconds()/3600 for t in tstamp]

fig = plt.figure(0)
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

fig = plt.figure(1)
plt.subplot(211)
plt.plot(tim, hmag, "-b")
plt.xlabel("Time [hr]")
plt.ylabel("Angular momentum")
plt.subplot(212)
plt.plot(tim, ener, "-b")
plt.xlabel("Time [hr]")
plt.ylabel("Specific energy")

fig = plt.figure(2)
plt.subplot(311)
plt.plot(tim, sma, "-b")
plt.xlabel("Time [hr]")
plt.ylabel("SMA [km]")
plt.subplot(312)
plt.plot(tim, ecc, "-b")
plt.xlabel("Time [hr]")
plt.ylabel("Eccentricity")
plt.subplot(313)
plt.plot(tim, inc, "-b")
plt.xlabel("Time [hr]")
plt.ylabel("Inclination [deg]")
plt.show()
