# ukf.py - Unscented Kalman filter implementation.
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

from numpy import *
from numpy.linalg import *
from orbdetpy import config
from .orekit import *
from .utils import *

def estimate(meas):
    frame = FramesFactory.getEME2000()
    gsta = stations()

    sdim = 6
    epoch = strtodate(config["Propagation"]["Start"])
    mass = config["SpaceObject"]["Mass"]
    prop = PropUtil(epoch, mass, frame, forces(True), sdim)

    X0 = array([config["Propagation"]["InitialState"]]).T
    P = diag(config["Estimation"]["Covariance"])
    Q = diag(config["Estimation"]["ProcessNoise"])
    R = zeros([len(config["Measurements"]), len(config["Measurements"])])
    for i, m in enumerate(config["Measurements"]):
        R[i,i] = config["Measurements"][m]["Error"]**2

    xhat = X0.copy()
    tm = epoch
    tofm = 0.0
    wfil = 0.5/sdim
    dt = 30
    Qst = block([[Q[:3,:3]*dt**2/4, Q[:3,:3]*dt/2],
                 [Q[:3,:3]*dt/2, Q[:3,:3]]])*dt**2
    sigm = zeros([sdim, 2*sdim])
    supd = zeros([len(config["Measurements"]), 2*sdim])
    results = []

    for midx in range(len(meas) + 1):
        if (midx < len(meas)):
            mea = meas[midx]
            tm, t0 = strtodate(mea["Time"]), tm
            res = {"Time" : mea["Time"], "PreFit" : {}, "PostFit" : {}}

            tmlt = AbsoluteDate(tm, -tofm)
            pvs = TimeStampedPVCoordinates(tmlt,
                                           Vector3D(xhat[0,0], xhat[1,0], xhat[2,0]),
                                           Vector3D(xhat[3,0], xhat[4,0], xhat[5,0]))
            pvi = gsta[mea["Station"]].getBaseFrame().getPVCoordinates(tm, frame)
            tofm, tof0 = Range.signalTimeOfFlight(pvs, pvi.getPosition(), tm), tofm
        else:
            tm, t0 = strtodate(config["Propagation"]["End"]), tm
            tofm, tof0 = 0.0, tofm

        sqrP = cholesky(P*sdim)
        for i in range(sdim):
            sigm[:,[i]] = xhat + sqrP[:,[i]]
            sigm[:,[sdim+i]] = xhat - sqrP[:,[i]]

        sppr = array([prop.propagate(t0.durationFrom(epoch) - tof0,
                                     sigm.ravel(order = "F").tolist(),
                                     tm.durationFrom(epoch) - tofm)]).reshape(
                                         (sdim, -1), order = "F")
        xhatpre = sppr.sum(axis = 1, keepdims = True)*wfil
        if (midx == len(meas)):
            break

        raw, obs = [], []
        for key, val in config["Measurements"].items():
            raw.append(mea[key])
            if (key == "Range"):
                obs.append(Range(gsta[mea["Station"]], tm, mea[key],
                                 val["Error"], 1.0, val["TwoWay"]))
            elif (key == "RangeRate"):
                obs.append(RangeRate(gsta[mea["Station"]], tm, mea[key],
                                     val["Error"], 1.0, val["TwoWay"]))

        Ppre = Qst.copy()
        tmlt = AbsoluteDate(tm, -tofm)
        for i in range(2*sdim):
            y = sppr[:,[i]] - xhatpre
            Ppre += outer(y, y)*wfil

            ssta = [SpacecraftState(CartesianOrbit(PVCoordinates(
                Vector3D(sppr[0,i], sppr[1,i], sppr[2,i]),
                Vector3D(sppr[3,i], sppr[4,i], sppr[5,i])),
                frame, tmlt, Constants.EGM96_EARTH_MU), mass)]
            for j, o in enumerate(obs):
                supd[j,i] = o.estimate(1, 1, ssta).getEstimatedValue()[0]

        Pyy = R.copy()
        Pxy = zeros([sdim, len(config["Measurements"])])
        yhatpre = supd.sum(axis = 1, keepdims = True)*wfil
        for i in range(2*sdim):
            y = supd[:,[i]] - yhatpre
            Pyy += outer(y, y)*wfil
            Pxy += outer(sppr[:,[i]] - xhatpre, y)*wfil

        K = Pxy.dot(inv(Pyy))
        xhat = xhatpre + K.dot(array([raw]).T - yhatpre)
        P = Ppre - K.dot(Pyy.dot(K.T))

        ssta = [SpacecraftState(CartesianOrbit(PVCoordinates(
            Vector3D(xhat[0,0], xhat[1,0], xhat[2,0]),
            Vector3D(xhat[3,0], xhat[4,0], xhat[5,0])),
            frame, tmlt, Constants.EGM96_EARTH_MU), mass)]
        res["EstimatedState"] = xhat[:,0].tolist()
        res["EstimatedCovariance"] = P.tolist()
        res["InnovationCovariance"] = Pyy.tolist()
        for i, m in enumerate(config["Measurements"]):
            res["PreFit"][m] = yhatpre[i,0]
            res["PostFit"][m] = obs[i].estimate(1, 1, ssta).getEstimatedValue()[0]
        results.append(res)

    return({"Estimation" : results,
            "Propagation" : {"Time" : config["Propagation"]["End"],
                             "State" : xhatpre[:,0].tolist()}})
