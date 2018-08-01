# plotres.py - Module to plot OD output.
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

tstamp, prefit, posfit, inocov = [], [], [], []
for i, o in zip(inp, out):
    tstamp.append(datetime.strptime(i["Time"], "%Y-%m-%dT%H:%M:%S.%fZ"))

    p = []
    for k, v in o["PreFit"].items():
        p.append(i[k] - v)
    prefit.append(p)

    p = []
    for k, v in o["PostFit"].items():
        p.append(i[k] - v)
    posfit.append(p)

    if ("InnovationCovariance" in o):
        p = []
        for j in range(len(o["InnovationCovariance"])):
            p.append(3.0*numpy.sqrt(o["InnovationCovariance"][j][j]))
        inocov.append(p)

pre = numpy.array(prefit)
pos = numpy.array(posfit)
cov = numpy.array(inocov)
key = list(cfg["Measurements"].keys())
tim = [(t - tstamp[0]).total_seconds()/3600 for t in tstamp]

fig = plt.figure(0)
plt.subplot(211)
plt.semilogx(tim, pre[:,0], "ob")
plt.xlabel("Time [hr]")
plt.ylabel("%s residual" % key[0])
plt.subplot(212)
plt.semilogx(tim, pre[:,1], "ob")
plt.xlabel("Time [hr]")
plt.ylabel("%s residual" % key[1])
plt.suptitle("Pre-fit residuals")

fig = plt.figure(1)
plt.subplot(211)
plt.semilogx(tim, pos[:,0], "ob")
if (len(cov) > 0):
    plt.semilogx(tim, -cov[:,0], "-r")
    plt.semilogx(tim,  cov[:,0], "-r")
plt.xlabel("Time [hr]")
plt.ylabel("%s residual" % key[0])
plt.subplot(212)
plt.semilogx(tim, pos[:,1], "ob")
if (len(cov) > 0):
    plt.semilogx(tim, -cov[:,1], "-r")
    plt.semilogx(tim,  cov[:,1], "-r")
plt.xlabel("Time [hr]")
plt.ylabel("%s residual" % key[1])
plt.suptitle("Post-fit residuals")

plt.show()
