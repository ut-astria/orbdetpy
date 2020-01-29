# estimation.py - Plot orbit determination results.
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

import math
import json
import numpy
import dateutil.parser
from numpy.linalg import norm
import matplotlib.pyplot as plt
import matplotlib.patches as patch
from orbdetpy import read_param
from orbdetpy.plotting import maximize_plot

def plot(config, measurements, orbit_fit, interactive = False, output_file_path = None):
    cfg = read_param(config)
    inp = [x for x in read_param(measurements) if (
        "Station" in x or "Position" in x or "PositionVelocity" in x)]
    out = [x for x in read_param(orbit_fit)
           if ("PreFit" in x and "PostFit" in x)]

    key = list(cfg["Measurements"].keys())
    dmcrun = (cfg["Estimation"].get("DMCCorrTime", 0.0) > 0.0 and
              cfg["Estimation"].get("DMCSigmaPert", 0.0) > 0.0)
    dmcidx = 6

    cd = cfg.get("Drag", {}).get("Coefficient", {}).get("Estimation", "")
    cr = cfg.get("RadiationPressure", {}).get("Creflection", {}).get("Estimation", "")
    if (cd == "Estimate"):
        dmcidx += 1
    if (cr == "Estimate"):
        dmcidx += 1

    parnames = []
    if (cd == "Estimate" or (cd == "Consider" and cr == "Consider")):
        parnames.extend([r"$C_D$", r"$C_R$"])
    else:
        parnames.extend([r"$C_R$", r"$C_D$"])

    idx = 0
    patches = []
    cycle = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    cmap = {None : cycle[0]}
    for sk, sv in cfg.get("Stations", {}).items():
        if (sk not in cmap):
            cmap[sk] = cycle[idx]
            patches.append(patch.Patch(color = cycle[idx], label = sk))
            idx = (idx + 1) % len(cycle)

        if (sv.get("BiasEstimation", "") in ["Estimate", "Consider"]):
            for m in cfg.get("Measurements", {}).keys():
                parnames.append(sk + m)

    tstamp,prefit,posfit,inocov,params,estmacc,estmcov,colors = [],[],[],[],[],[],[],[]
    for i, o in zip(inp, out):
        colors.append(cmap[o.get("Station")])
        tstamp.append(dateutil.parser.isoparse(i["Time"]))

        if (key[0] == "Position" or key[0] == "PositionVelocity"):
            prefit.append([ix-ox for ix, ox
                           in zip(i[key[0]], o["PreFit"][key[0]])])
            posfit.append([ix-ox for ix, ox
                           in zip(i[key[0]], o["PostFit"][key[0]])])
        else:
            if (len(key) == 2):
                prefit.append([i[key[0]]-o["PreFit"][key[0]][-1],
                               i[key[1]]-o["PreFit"][key[1]][-1]])
                posfit.append([i[key[0]]-o["PostFit"][key[0]][-1],
                               i[key[1]]-o["PostFit"][key[1]][-1]])
            else:
                prefit.append([i[key[0]]-o["PreFit"][key[0]][-1]])
                posfit.append([i[key[0]]-o["PostFit"][key[0]][-1]])

        p = []
        for m in range(len(o.get("InnovationCovariance", []))):
            p.append(3.0*numpy.sqrt(o["InnovationCovariance"][m][m]))
        if (len(p) > 0):
            inocov.append(p)
        else:
            inocov.append([0.0]*len(prefit[0]))

        if (len(o["EstimatedState"]) > 6):
            if (dmcrun):
                params.append(o["EstimatedState"][6:dmcidx] +
                              o["EstimatedState"][dmcidx+3:])
            else:
                params.append(o["EstimatedState"][6:])

        if (dmcrun):
            r = numpy.array(o["EstimatedState"][:3])
            r /= norm(r)
            v = numpy.array(o["EstimatedState"][3:6])
            v /= norm(v)
            h = numpy.cross(r, v)
            rot = numpy.vstack((r, numpy.cross(h, r), h))
            estmacc.append(rot.dot(o["EstimatedState"][dmcidx:dmcidx+3]))

    pre = numpy.array(prefit)
    pos = numpy.array(posfit)
    cov = numpy.array(inocov)
    par = numpy.array(params)
    estmacc = numpy.array(estmacc)
    start = tstamp[0] if (tstamp[0] < tstamp[-1]) else tstamp[-1]
    tim = [(t - start).total_seconds()/3600 for t in tstamp]

    angles = ["Azimuth", "Elevation", "RightAscension", "Declination"]
    if (key[0] in angles and key[1] in angles):
        pre *= 648000.0/math.pi
        pos *= 648000.0/math.pi
        cov *= 648000.0/math.pi
        units = ["arcsec", "arcsec"]
    else:
        if ("Position" in key):
            units = ["m", "m", "m"]
        elif ("PositionVelocity" in key):
            units = ["m", "m", "m", "m/s", "m/s", "m/s"]
        else:
            units = ["m", "m/s"]

    if ("Position" in key):
        ylabs = [r"$\Delta x$", r"$\Delta y$", r"$\Delta z$"]
        order = [1, 2, 3]
    elif ("PositionVelocity" in key):
        ylabs = [r"$\Delta x$", r"$\Delta y$", r"$\Delta z$",
                 r"$\Delta v_x$", r"$\Delta v_y$", r"$\Delta v_z$"]
        order = [1, 3, 5, 2, 4, 6]
    else:
        ylabs = key

    outfiles = []
    plt.figure(0)
    maximize_plot()
    plt.suptitle("Pre-fit residuals")
    for i in range(pre.shape[-1]):
        if ("Position" in key):
            plt.subplot(3, 1, order[i])
        elif ("PositionVelocity" in key):
            plt.subplot(3, 2, order[i])
        else:
            plt.subplot(pre.shape[-1], 1, i + 1)
        plt.scatter(tim, pre[:,i], color = colors, marker = "o", s = 7)
        plt.legend(handles = patches, loc = "best")
        plt.xlabel("Time [hr]")
        plt.ylabel("%s [%s]" % (ylabs[i], units[i]))

    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (output_file_path):
        outfiles.append(output_file_path + "_prefit.png")
        plt.savefig(outfiles[-1], format = "png")

    plt.figure(1) 
    maximize_plot()
    plt.suptitle("Post-fit residuals")
    for i in range(pre.shape[-1]):
        if ("Position" in key):
            plt.subplot(3, 1, order[i])
        elif ("PositionVelocity" in key):
            plt.subplot(3, 2, order[i])
        else:
            plt.subplot(pre.shape[-1], 1, i + 1)
        plt.scatter(tim, pos[:,i], color = colors, marker = "o", s = 7)
        plt.legend(handles = patches, loc = "best")
        plt.plot(tim, -cov[:,i], "-r")
        plt.plot(tim,  cov[:,i], "-r")
        plt.xlabel("Time [hr]")
        plt.ylabel("%s [%s]" % (ylabs[i], units[i]))

    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (output_file_path):
        outfiles.append(output_file_path + "_postfit.png")
        plt.savefig(outfiles[-1], format = "png")

    for i in range(par.shape[-1]):
        if (i == 0):
            plt.figure(2)
            maximize_plot()
            plt.suptitle("Estimated parameters")
        plt.subplot(par.shape[1], 1, i + 1)
        plt.scatter(tim, par[:,i], marker = "o", s = 7)
        plt.xlabel("Time [hr]")
        plt.ylabel(parnames[i])

    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (output_file_path):
        outfiles.append(output_file_path + "_estpar.png")
        plt.savefig(outfiles[-1], format = "png")

    if (dmcrun):
        lab = [r"Radial [$\frac{m}{s^2}$]", r"In track [$\frac{m}{s^2}$]",
               r"Cross track [$\frac{m}{s^2}$]"]
        plt.figure(3)
        maximize_plot()
        plt.suptitle("DMC estimated accelerations")
        for i in range(3):
            plt.subplot(3, 1, i+1)
            plt.plot(tim, estmacc[:,i])
            plt.xlabel("Time [hr]")
            plt.ylabel(lab[i])

        plt.tight_layout(rect = [0, 0.03, 1, 0.95])
        if (output_file_path):
            outfiles.append(output_file_path + "_estacc.png")
            plt.savefig(outfiles[-1], format = "png")

    if (interactive):
        plt.show()
    plt.close("all")
    return(outfiles)
