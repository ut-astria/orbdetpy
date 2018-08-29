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
        self.angles = ["Azimuth", "Elevation",
                       "RightAscension", "Declination"]

    @jnius.java_method(
        "(Lorg/orekit/estimation/sequential/KalmanEstimation;)V")
    def evaluationPerformed(self, est):
        keys = list(self.config["Measurements"].keys())
        if (keys[0] in self.angles):
            n = est.getCurrentMeasurementNumber() - 1
        else:
            n = (est.getCurrentMeasurementNumber() - 1)//len(keys)

        if (len(self.results) <= n):
            k = keys[0]
            self.results.append({"Time" : self.meas[n]["Time"],
                                 "PreFit" : {}, "PostFit" : {}})
        else:
            k = keys[1]

        res = self.results[n]
        fitv = est.getPredictedMeasurement().getEstimatedValue()
        if (keys[0] in self.angles):
            res["PreFit"][keys[0]] = fitv[0]
            res["PreFit"][keys[1]] = fitv[1]
        else:
            res["PreFit"][k] = fitv[0]

        fitv = est.getCorrectedMeasurement().getEstimatedValue()
        if (keys[0] in self.angles):
            res["PostFit"][keys[0]] = fitv[0]
            res["PostFit"][keys[1]] = fitv[1]
        else:
            res["PostFit"][k] = fitv[0]

        res["EstimatedState"] = pvtolist(est.getPredictedSpacecraftStates()[0].getPVCoordinates())
        plst = est.getEstimatedPropagationParameters().getDrivers().toArray() + \
               est.getEstimatedMeasurementsParameters().getDrivers().toArray()
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
    for mea in meas:
        tm = strtodate(mea["Time"])
        for key, val in config["Measurements"].items():
            if (key == "Range"):
                allobs.add(Range(gsta[mea["Station"]], tm, mea[key],
                                 val["Error"], 1.0, val["TwoWay"]))
            elif (key == "RangeRate"):
                allobs.add(RangeRate(gsta[mea["Station"]], tm, mea[key],
                                     val["Error"], 1.0, val["TwoWay"]))
            elif (key in ["Azimuth", "Elevation"]):
                allobs.add(AngularAzEl(gsta[mea["Station"]], tm,
                                       [mea["Azimuth"], mea["Elevation"]],
                                       [config["Measurements"]["Azimuth"]["Error"],
                                        config["Measurements"]["Elevation"]["Error"]],
                                       [1.0, 1.0]))
                break
            elif (key in ["RightAscension", "Declination"]):
                allobs.add(AngularRaDec(gsta[mea["Station"]], frame, tm,
                                        [mea["RightAscension"], mea["Declination"]],
                                        [config["Measurements"]["RightAscension"]["Error"],
                                         config["Measurements"]["Declination"]["Error"]],
                                        [1.0, 1.0]))
                break

    est = filt.processMeasurements(allobs)[0]
    pvend = pvtolist(est.getPVCoordinates(
        strtodate(config["Propagation"]["End"]), frame))
    if (sdim > 6):
        pvend.extend(cbak.results[-1]["EstimatedState"][6:])

    return({"Estimation" : cbak.results,
            "Propagation" : {"Time" : config["Propagation"]["End"],
                             "State" : pvend}})
