# plotodet.py - Module to plot OD output.
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

import sys
import math
import json
import numpy
from datetime import datetime
import matplotlib.pyplot as plt

if (len(sys.argv) < 4):
    print("Usage: python %s config_file measurement_file output_file"
          % sys.argv[0])
    exit()

with open(sys.argv[1], "r") as f:
    cfg = json.load(f)
with open(sys.argv[2], "r") as f:
    inp = json.load(f)
with open(sys.argv[3], "r") as f:
    out = json.load(f)["Estimation"]

key = list(cfg["Measurements"].keys())
dmcrun = (cfg["Estimation"].get("DMCCorrTime", 0.0) > 0.0 and
          cfg["Estimation"].get("DMCSigmaPert", 0.0) > 0.0)

tstamp, prefit, posfit, inocov, params, estmacc = [], [], [], [], [], []
for i, o in zip(inp, out):
    tstamp.append(datetime.strptime(i["Time"], "%Y-%m-%dT%H:%M:%S.%fZ"))

    prefit.append([i[key[0]] - o["PreFit"][key[0]],
                   i[key[1]] - o["PreFit"][key[1]]])
    posfit.append([i[key[0]] - o["PostFit"][key[0]],
                   i[key[1]] - o["PostFit"][key[1]]])

    inocov.append([3.0*numpy.sqrt(o["InnovationCovariance"][0][0]),
                   3.0*numpy.sqrt(o["InnovationCovariance"][1][1])])

    if (len(o["EstimatedState"]) > 6):
        if (dmcrun):
            params.append(o["EstimatedState"][6:-3])
        else:
            params.append(o["EstimatedState"][6:])

    if (dmcrun):
        estmacc.append(o["EstimatedState"][-3:])

pre = numpy.array(prefit)
pos = numpy.array(posfit)
cov = numpy.array(inocov)
par = numpy.array(params)
estmacc = numpy.array(estmacc)
tim = [(t - tstamp[0]).total_seconds()/3600 for t in tstamp]

angles = ["Azimuth", "Elevation", "RightAscension", "Declination"]
if (key[0] in angles and key[1] in angles):
    pre *= 648000.0/math.pi
    pos *= 648000.0/math.pi
    cov *= 648000.0/math.pi
    units = ["arcsec", "arcsec"]
else:
    units = ["m", "m/s"]

fig = plt.figure(0)
plt.subplot(211)
plt.semilogx(tim, pre[:,0], "ob")
plt.xlabel("Time [hr]")
plt.ylabel("%s residual [%s]" % (key[0], units[0]))
plt.subplot(212)
plt.semilogx(tim, pre[:,1], "ob")
plt.xlabel("Time [hr]")
plt.ylabel("%s residual [%s]" % (key[1], units[1]))
plt.suptitle("Pre-fit residuals")

fig = plt.figure(1)
plt.subplot(211)
plt.semilogx(tim, pos[:,0], "ob")
plt.semilogx(tim, -cov[:,0], "-r")
plt.semilogx(tim,  cov[:,0], "-r", label = r"Innov. 3$\sigma$")
plt.xlabel("Time [hr]")
plt.ylabel("%s residual [%s]" % (key[0], units[0]))
plt.legend(loc = "upper left")
plt.subplot(212)
plt.semilogx(tim, pos[:,1], "ob")
plt.semilogx(tim, -cov[:,1], "-r")
plt.semilogx(tim,  cov[:,1], "-r", label = r"Innov. 3$\sigma$")
plt.xlabel("Time [hr]")
plt.ylabel("%s residual [%s]" % (key[1], units[1]))
plt.legend(loc = "upper left")
plt.suptitle("Post-fit residuals")

for i in range(par.shape[-1]):
    if (i == 0):
        fig = plt.figure(2)
    plt.subplot(par.shape[1], 1, i + 1)
    plt.semilogx(tim, par[:,i], "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("Parameter %d" % (i + 1))

if (dmcrun):
    fig = plt.figure(3)
    lab = [r"Radial [$\frac{m}{s^2}$]", r"In track [$\frac{m}{s^2}$]",
           r"Cross track [$\frac{m}{s^2}$]"]
    for i in range(3):
        plt.subplot(3, 1, i+1)
        plt.semilogx(tim, estmacc[:,i], "-b")
        plt.xlabel("Time [hr]")
        plt.ylabel(lab[i])

    plt.suptitle("DMC estimated acceleration")

plt.show()
