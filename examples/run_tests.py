# run_tests.py - Program to run simulation and OD tests.
# Copyright (C) 2019-2020 University of Texas
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

import os
import sys
import json
import math
import time
from orbdetpy.estimation import determine_orbit
from orbdetpy.simulation import simulate_measurements

print("run_tests start : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
odpdir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
datdir = os.path.join(odpdir, "examples", "data")

sim_cfg, sim_out = [], []
for root, dirs, files in os.walk(datdir):
    outdir = os.path.join(root, "output")
    if ("output" in dirs):
        dirs.remove("output")
    else:
        os.makedirs(outdir)

    for fname in files:
        idx = fname.find("sim_cfg.json")
        if (idx != -1):
            print("Simulating {}".format(fname))
            sim_cfg.append(os.path.join(root, fname))
            sim_out.append(os.path.join(outdir, fname[:idx] + "obs.json"))

simulate_measurements(sim_cfg, output_file = sim_out)

od_cfg, od_obs, od_out = [], [], []
for root, dirs, files in os.walk(datdir):
    outdir = os.path.join(root, "output")
    if ("output" in dirs):
        dirs.remove("output")

    for fname in files:
        idx = fname.find("od_cfg.json")
        if (idx == -1):
            continue

        print("Fitting {}".format(fname))
        for algo in ["EKF", "UKF"]:
            with open(os.path.join(root, fname), "r") as fp:
                config = json.load(fp)
                config["estmFilter"] = algo
            od_cfg.append(config)
            od_obs.append(os.path.join(outdir, fname[:idx] + "obs.json"))
            od_out.append(os.path.join(outdir, fname[:idx] + algo + "_fit.json"))

determine_orbit(od_cfg, od_obs, output_file = od_out)

for root, dirs, files in os.walk(datdir):
    outdir = os.path.join(root, "output")
    if ("output" in dirs):
        dirs.remove("output")

    for fname in files:
        idx = fname.find("od_cfg.json")
        if (idx == -1):
            continue

        with open(os.path.join(outdir, fname[:idx] + "obs.json"), "r") as fp:
            obs = json.load(fp)
        exact = [x for x in obs if ("station" in x or "positionVelocity" in x)]

        for algo in ["EKF", "UKF"]:
            output = []
            with open(os.path.join(outdir, fname[:idx] + algo + "_fit.json"), "r") as fp:
                fit = json.load(fp)
            estim = [x for x in fit if ("preFit" in x and "postFit" in x)]

            for exa, est in zip(exact, estim):
                if ("trueStateCartesian" in exa):
                    X = exa["trueStateCartesian"][:6]
                else:
                    X = exa["positionVelocity"][:6]
                diff = [X[i] - est["estimatedState"][i] for i in range(6)]
                sigm = [3.0*math.sqrt(est["estimatedCovariance"][i][i]) for i in range(6)]
                output.append({"Time" : exa["time"], "Station" : exa.get("station"),
                               "StateResidual" : diff, "Covariance3Sigma" : sigm,
                               "WithinBounds" : [diff[i] > -sigm[i] and diff[i] < sigm[i] for i in range(6)]})

            with open(os.path.join(outdir, fname[:idx] + algo + "_diff.json"), "w") as fp:
                json.dump(output, fp, indent = 1)

print("run_tests end   : %s" % time.strftime("%Y-%m-%d %H:%M:%S"))
