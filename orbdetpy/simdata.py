# simdata.py - Simulate data for a spacecraft's orbit.
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

from numpy.random import randn
import jnius
from orbdetpy import config
from .orekit import *
from .utils import *

def simulate():
    gsta = stations()
    tm = strtodate(config["Propagation"]["Start"])
    prend = strtodate(config["Propagation"]["End"])
    tstep = config["Propagation"]["Step"]
    X0 = config["Propagation"]["InitialState"]

    prop = NumericalPropagator(DormandPrince853Integrator(
        1E-3, 300.0, 1E-14, 1E-12))
    prop.setInitialState(SpacecraftState(CartesianOrbit(PVCoordinates(
        Vector3D(X0[:3]), Vector3D(X0[3:6])), FramesFactory.getEME2000(),
        tm, Constants.EGM96_EARTH_MU), config["SpaceObject"]["Mass"]))
    for f in forces(False):
        prop.addForceModel(f)

    res = []
    while True:
        ssta = [prop.propagate(tm)]
        meas = {"Time" : tm.toString() + "Z",
                "State" : pvtolist(ssta[0].getPVCoordinates())}
        for gstr, gs in gsta.items():
            azel = AngularAzEl(gs, tm, [0.0, 0.0], [0.1, 0.1], [1.0,1.0])
            if (azel.estimate(1, 1, ssta).getEstimatedValue()[1] <= 0.0):
                continue

            meas["Station"] = gstr
            for key, val in config["Measurements"].items():
                if (key == "Range"):
                    obs = Range(gs, tm, 0.0, val["Error"],
                                1.0, val["TwoWay"])
                elif (key == "RangeRate"):
                    obs = RangeRate(gs, tm, 0.0, val["Error"],
                                    1.0, val["TwoWay"])

                meas[key] = obs.estimate(1, 1, ssta).getEstimatedValue()[0] + \
                            randn()*val["Error"]
            break

        if (len(meas) > 2):
            res.append(meas)

        dt = prend.durationFrom(tm)
        tm = AbsoluteDate(tm, min(dt, tstep))
        if (dt <= 0):
            break

    return(res)
    
