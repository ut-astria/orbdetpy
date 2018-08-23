# ekf.py - Wrapper for Orekit's Extended Kalman Filter.
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

if __name__ == "__main__":
    exit()

import jnius
from .orekit import *
from .utils import *

class observer(jnius.PythonJavaClass):
    __javacontext__ = "app"
    __javainterfaces__ = [
        "org/orekit/estimation/sequential/KalmanObserver"]

    def __init__(self, conf, m, par):
        super().__init__()
        self.config = conf
        self.meas = m
        self.pest = par
        self.results = []

    @jnius.java_method(
        "(Lorg/orekit/estimation/sequential/KalmanEstimation;)V")
    def evaluationPerformed(self, est):
        n = (est.getCurrentMeasurementNumber() - 1)//2
        plst = est.getEstimatedPropagationParameters().getDrivers().toArray() + \
               est.getEstimatedMeasurementsParameters().getDrivers().toArray()

        if (len(self.results) <= n):
            k = list(self.config["Measurements"].keys())[0]
            self.results.append({"Time" : self.meas[n]["Time"],
                                 "PreFit" : {}, "PostFit" : {}})
        else:
            k = list(self.config["Measurements"].keys())[1]

        res = self.results[n]
        res["PreFit"][k] = est.getPredictedMeasurement().getEstimatedValue()[0]
        res["PostFit"][k] = est.getCorrectedMeasurement().getEstimatedValue()[0]
        res["EstimatedState"] = pvtolist(est.getPredictedSpacecraftStates()[0].getPVCoordinates())
        for p in self.pest:
            for l in plst:
                if l.getName() == p:
                    res["EstimatedState"].append(l.getValue())
        res["EstimatedCovariance"] = est.getPhysicalEstimatedCovarianceMatrix().getData()

def estimate(config, meas):
    frame = FramesFactory.getEME2000()
    gsta = stations(config)

    X0 = config["Propagation"]["InitialState"]
    X0 = CartesianOrbit(PVCoordinates(Vector3D(X0[:3]), Vector3D(X0[3:6])),
                        frame, strtodate(config["Propagation"]["Start"]),
                        Constants.EGM96_EARTH_MU)

    prop = NumericalPropagatorBuilder(X0, DormandPrince853IntegratorBuilder(
        1E-3, 300.0, 1.0), PositionAngle.MEAN, 10.0)
    prop.setMass(config["SpaceObject"]["Mass"])
    for f in forces(config, False):
        prop.addForceModel(f)

    sdim, parm, pest = estparms(config)
    plst = prop.getPropagationParametersDrivers()
    for n, l in zip(pest, parm.tolist()):
        pdrv = ParameterDriver(String(n), l[2], 1.0, l[0], l[1])
        pdrv.setSelected(True)
        plst.add(pdrv)

    build = KalmanEstimatorBuilder()
    build.addPropagationConfiguration(prop, ConstantProcessNoise(
        DiagonalMatrix(config["Estimation"]["Covariance"]),
        DiagonalMatrix(config["Estimation"]["ProcessNoise"])))

    cbak = observer(config, meas, pest)
    filt = build.build()
    filt.setObserver(cbak)

    allobs = ArrayList()
    for m in meas:
        tm = strtodate(m["Time"])
        for k,v in config["Measurements"].items():
            if (k == "Range"):
                obs = Range(gsta[m["Station"]], tm, m[k],
                            v["Error"], 1.0, v["TwoWay"])
            elif (k == "RangeRate"):
                obs = RangeRate(gsta[m["Station"]], tm, m[k],
                                v["Error"], 1.0, v["TwoWay"])

            allobs.add(obs)

    est = filt.processMeasurements(allobs)[0]
    pvend = pvtolist(est.getPVCoordinates(
        strtodate(config["Propagation"]["End"]), frame))
    if (sdim > 6):
        pvend.extend(cbak.results[-1]["EstimatedState"][6:])

    return({"Estimation" : cbak.results,
            "Propagation" : {"Time" : config["Propagation"]["End"],
                             "State" : pvend}})
