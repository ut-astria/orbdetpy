# plotodet.py - Module to plot OD output.
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


def plot(cfgfile, inpfile, outfile, interactive = False, filepath = None):
    with open(cfgfile, "r") as fp:
        cfg = json.load(fp)
    with open(inpfile, "r") as fp:
        inp = json.load(fp)
    with open(outfile, "r") as fp:
        out = json.load(fp)["Estimation"]

    key = tuple(cfg["Measurements"].keys())
    dmcrun = (cfg["Estimation"].get("DMCCorrTime", 0.0) > 0.0 and
              cfg["Estimation"].get("DMCSigmaPert", 0.0) > 0.0)

    tstamp, prefit, posfit, inocov, params, estmacc, estmcov = [], [], [], [], [], [], []
    for i, o in zip(inp, out):
        tstamp.append(dateutil.parser.isoparse(i["Time"]))

        if ("PositionVelocity" in key):
            prefit.append([ix - ox for ix, ox in zip(i["PositionVelocity"],
                                                     o["PreFit"]["PositionVelocity"])])
            posfit.append([ix - ox for ix, ox in zip(i["PositionVelocity"],
                                                     o["PostFit"]["PositionVelocity"])])
        else:
            prefit.append([i[key[0]] - o["PreFit"][key[0]][-1],
                           i[key[1]] - o["PreFit"][key[1]][-1]])
            posfit.append([i[key[0]] - o["PostFit"][key[0]][-1],
                           i[key[1]] - o["PostFit"][key[1]][-1]])

        p = []
        for m in range(len(o["InnovationCovariance"])):
            p.append(3.0*numpy.sqrt(o["InnovationCovariance"][m][m]))
        inocov.append(p)

        if (len(o["EstimatedState"]) > 6):
            if (dmcrun):
                params.append(o["EstimatedState"][6:-3])
            else:
                params.append(o["EstimatedState"][6:])

        if (dmcrun):
            r = numpy.array(o["EstimatedState"][:3])
            r /= norm(r)
            v = numpy.array(o["EstimatedState"][3:6])
            v /= norm(v)
            h = numpy.cross(r, v)
            rot = numpy.vstack((r, numpy.cross(h, r), h))
            estmacc.append(rot.dot(o["EstimatedState"][-3:]))

    pre = numpy.array(prefit)
    pos = numpy.array(posfit)
    cov = numpy.array(inocov)
    par = numpy.array(params)
    estmacc = numpy.array(estmacc)
    tim = [(t - tstamp[0]).total_seconds()/3600 for t in tstamp]

    angles = ("Azimuth", "Elevation", "RightAscension", "Declination")
    if (key[0] in angles and key[1] in angles):
        pre *= 648000.0/math.pi
        pos *= 648000.0/math.pi
        cov *= 648000.0/math.pi
        units = ("arcsec", "arcsec")
    else:
        if ("PositionVelocity" in key):
            units = ("m", "m", "m", "m/s", "m/s", "m/s")
        else:
            units = ("m", "m/s")

    if ("PositionVelocity" in key):
        ylabs = (r"$\Delta x$", r"$\Delta y$", r"$\Delta z$",
                 r"$\Delta v_x$", r"$\Delta v_y$", r"$\Delta v_z$")
        order = (1, 3, 5, 2, 4, 6)
    else:
        ylabs = key

    outfiles = []

    plt.figure(0)
    plt.suptitle("Pre-fit residuals")
    for i in range(pre.shape[-1]):
        if ("PositionVelocity" in key):
            plt.subplot(3, 2, order[i])
        else:
            plt.subplot(2, 1, i + 1)
        plt.semilogx(tim, pre[:,i], "ob")
        plt.xlabel("Time [hr]")
        plt.ylabel("%s [%s]" % (ylabs[i], units[i]))

    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (filepath is not None):
        outfiles.append(filepath + "_prefit.png")
        plt.savefig(outfiles[-1], format = "png")

    plt.figure(1) 
    plt.suptitle("Post-fit residuals")
    for i in range(pre.shape[-1]):
        if ("PositionVelocity" in key):
            plt.subplot(3, 2, order[i])
        else:
            plt.subplot(2, 1, i + 1)
        plt.semilogx(tim, pos[:,i], "ob")
        plt.semilogx(tim, -cov[:,i], "-r")
        plt.semilogx(tim,  cov[:,i], "-r", label = r"Innov. 3$\sigma$")
        plt.xlabel("Time [hr]")
        plt.ylabel("%s [%s]" % (ylabs[i], units[i]))
        if ("PositionVelocity" not in key):
            plt.ylim(-cov[i,0], cov[i,0])

    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (filepath is not None):
        outfiles.append(filepath + "_postfit.png")
        plt.savefig(outfiles[-1], format = "png")

    parnames, parmvals = [], []
    if (cfg["Drag"]["Coefficient"]["Estimation"] == "Estimate"):
        parnames.append(r"$C_D$")

    if (cfg["RadiationPressure"]["Creflection"]["Estimation"] == "Estimate"):
        parnames.append(r"$C_R$")

    for i in range(par.shape[-1]):
        if (i == 0):
            plt.figure(2)
            plt.suptitle("Estimated parameters")

        plt.subplot(par.shape[1], 1, i + 1)
        plt.semilogx(tim, par[:,i], "ob")
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
            plt.semilogx(tim, estmacc[:,i], "-b")
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


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description='Plots orbit determination results.')
    parser.add_argument('config', help='Path to config file.')
    parser.add_argument('input', help='Path to input file.')
    parser.add_argument('output', help='Path to output file.')
    parser.add_argument('--interactive', default=True,
                        help='Whether to show plots interactively.')
    args = parser.parse_args()
    plot(args.config, args.input, args.output, interactive=args.interactive)
