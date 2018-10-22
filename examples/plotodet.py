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

def plot(cfgfile, inpfile, outfile, interactive = False, filepath = None):
    with open(cfgfile, "r") as fp:
        cfg = json.load(fp)
    with open(inpfile, "r") as fp:
        inp = json.load(fp)
    with open(outfile, "r") as fp:
        out = json.load(fp)["Estimation"]

    key = list(cfg["Measurements"].keys())
    dmcrun = (cfg["Estimation"].get("DMCCorrTime", 0.0) > 0.0 and
              cfg["Estimation"].get("DMCSigmaPert", 0.0) > 0.0)

    tstamp, prefit, posfit, inocov, params, estmacc, estmcov = [], [], [], [], [], [], []
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

    #        p = []
    #        for m in range(6, len(o["EstimatedState"])):
    #            p.append(3.0*numpy.sqrt(o["EstimatedCovariance"][m][m]))
    #        estmcov.append(p)

        if (dmcrun):
            estmacc.append(o["EstimatedState"][-3:])

    pre = numpy.array(prefit)
    pos = numpy.array(posfit)
    cov = numpy.array(inocov)
    par = numpy.array(params)
    estmacc = numpy.array(estmacc)
    #estmcov = numpy.array(estmcov)
    tim = [(t - tstamp[0]).total_seconds()/3600 for t in tstamp]

    angles = ["Azimuth", "Elevation", "RightAscension", "Declination"]
    if (key[0] in angles and key[1] in angles):
        pre *= 648000.0/math.pi
        pos *= 648000.0/math.pi
        cov *= 648000.0/math.pi
        units = ["arcsec", "arcsec"]
    else:
        units = ["m", "m/s"]

    outfiles = []

    plt.figure(0)
    plt.suptitle("Pre-fit residuals")
    plt.subplot(211)
    plt.semilogx(tim, pre[:,0], "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("%s [%s]" % (key[0], units[0]))
    plt.subplot(212)
    plt.semilogx(tim, pre[:,1], "ob")
    plt.xlabel("Time [hr]")
    plt.ylabel("%s [%s]" % (key[1], units[1]))
    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (filepath is not None):
        outfiles.append(filepath + "_prefit.png")
        plt.savefig(outfiles[-1], format = "png")

    plt.figure(1) 
    plt.suptitle("Post-fit residuals")
    plt.subplot(211)
    plt.semilogx(tim, pos[:,0], "ob")
    plt.semilogx(tim, -cov[:,0], "-r")
    plt.semilogx(tim,  cov[:,0], "-r", label = r"Innov. 3$\sigma$")
    plt.xlabel("Time [hr]")
    plt.ylabel("%s [%s]" % (key[0], units[0]))
    plt.legend(loc = "best")
    plt.ylim(-cov[1,0], cov[1,0])
    plt.subplot(212)
    plt.semilogx(tim, pos[:,1], "ob")
    plt.semilogx(tim, -cov[:,1], "-r")
    plt.semilogx(tim,  cov[:,1], "-r", label = r"Innov. 3$\sigma$")
    plt.xlabel("Time [hr]")
    plt.ylabel("%s [%s]" % (key[1], units[1]))
    plt.legend(loc = "best")
    plt.ylim(-cov[1,1], cov[1,1])
    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (filepath is not None):
        outfiles.append(filepath + "_postfit.png")
        plt.savefig(outfiles[-1], format = "png")

    parnames, parmvals = [], []
    if (cfg["Drag"]["Coefficient"]["Estimation"] == "Estimate"):
        parnames.append(r"$C_D$")
#        parmvals.append(cfg["Drag"]["Coefficient"]["Value"])

    if (cfg["RadiationPressure"]["Creflection"]["Estimation"] == "Estimate"):
        parnames.append(r"$C_R$")
#        parmvals.append(cfg["RadiationPressure"]["Creflection"]["Value"])

    for i in range(par.shape[-1]):
        if (i == 0):
            plt.figure(2)
            plt.suptitle("Estimated parameters")

        plt.subplot(par.shape[1], 1, i + 1)
        plt.semilogx(tim, par[:,i], "ob")
    #    plt.semilogx(tim, -parmvals[i] - estmcov[:,i], "-r")
    #    plt.semilogx(tim,  parmvals[i] + estmcov[:,i], "-r", label = r"Covariance 3$\sigma$")
        plt.xlabel("Time [hr]")
        plt.ylabel(parnames[i])

    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (filepath is not None):
        outfiles.append(filepath + "_estpar.png")
        plt.savefig(outfiles[-1], format = "png")

    if (dmcrun):
        plt.figure(3)
        plt.suptitle("DMC estimated accelerations")

        lab = [r"Radial [$\frac{m}{s^2}$]", r"In track [$\frac{m}{s^2}$]",
               r"Cross track [$\frac{m}{s^2}$]"]
        for i in range(3):
            plt.subplot(3, 1, i+1)
            plt.semilogx(tim, estmacc[:,i], "ob")
    #        plt.semilogx(tim, -estmcov[:,i-3], "-r")
    #        plt.semilogx(tim,  estmcov[:,i-3], "-r", label = r"Covariance 3$\sigma$")
            plt.xlabel("Time [hr]")
            plt.ylabel(lab[i])

        plt.tight_layout(rect = [0, 0.03, 1, 0.95])
        if (filepath is not None):
            outfiles.append(filepath + "_estacc.png")
            plt.savefig(outfiles[-1], format = "png")

    if (interactive):
        plt.show()

    plt.close("all")
    return(outfiles)

if (__name__ == "__main__"):
    if (len(sys.argv) < 4):
        print("Usage: python %s config_file measurement_file output_file" % sys.argv[0])
        exit()

    plot(sys.argv[1], sys.argv[2], sys.argv[3], interactive = True)
