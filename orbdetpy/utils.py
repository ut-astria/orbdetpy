# utils.py - Miscellaneous utility routines.
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

from .orekit import *

def forces(config, pointmass):
    ut1scale = TimeScalesFactory.getUT1(IERSConventions.IERS_2010, False)

    fmod = []
    if (pointmass):
        fmod.append(NewtonianAttraction(Constants.EGM96_EARTH_MU))

    grav = GravityFieldFactory.getNormalizedProvider(
        config["Gravity"]["Degree"], config["Gravity"]["Order"])
    fmod.append(HolmesFeatherstoneAttractionModel(_itrf, grav))

    if (config["OceanTides"]["Degree"] >= 0 and
        config["OceanTides"]["Order"] >= 0):
        fmod.append(OceanTides(_itrf,
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            Constants.EGM96_EARTH_MU, config["OceanTides"]["Degree"],
            config["OceanTides"]["Order"],
            IERSConventions.IERS_2010, ut1scale))

    if (config["Drag"]["Model"].lower() == "exponential"):
        fmod.append(DragForce(SimpleExponentialAtmosphere(
            OneAxisEllipsoid(
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            Constants.WGS84_EARTH_FLATTENING, _itrf),
            config["Drag"]["ExpRho0"], config["Drag"]["ExpH0"],
            config["Drag"]["ExpHScale"]),
            IsotropicDrag(config["SpaceObject"]["Area"],
            config["Drag"]["Coeff"]["Value"])))

    if (config["SolidTides"]["Sun"] or config["SolidTides"]["Moon"]):
        fmod.append(PropUtil.solidtides(_itrf,
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            Constants.EGM96_EARTH_MU, grav.getTideSystem(),
            IERSConventions.IERS_2010, ut1scale,
            config["SolidTides"]["Sun"], config["SolidTides"]["Moon"]))

    if (config["ThirdBodies"]["Sun"]):
        fmod.append(ThirdBodyAttraction(CelestialBodyFactory.getSun()))

    if (config["ThirdBodies"]["Moon"]):
        fmod.append(ThirdBodyAttraction(CelestialBodyFactory.getMoon()))

    if (config["RadiationPressure"]["Sun"]):
        fmod.append(SolarRadiationPressure(
            149597870000.0, 4.56E-6, CelestialBodyFactory.getSun(),
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            IsotropicRadiationClassicalConvention(
            config["SpaceObject"]["Area"], config["SpaceObject"]["Cabs"],
            config["SpaceObject"]["Cspec"])))

    mans = config.get("Maneuvers", [])
    for man in mans:
        fmod.append(ConstantThrustManeuver(strtodate(man["Time"]),
            man["Duration"], man["Thrust"], man["Isp"],
            Vector3D(man["Direction"])))

    return(fmod)

def stations(config):
    gsta = {}
    dt = strtodate("2000-01-01T12:00:00Z")
    for k, v in config["Stations"].items():
        gsta[k] = GroundStation(TopocentricFrame(OneAxisEllipsoid(
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            Constants.WGS84_EARTH_FLATTENING, _itrf),
            GeodeticPoint(v["Latitude"], v["Longitude"], v["Altitude"]), k))
        gsta[k].getPrimeMeridianOffsetDriver().setReferenceDate(dt)
        gsta[k].getPolarOffsetXDriver().setReferenceDate(dt)
        gsta[k].getPolarOffsetYDriver().setReferenceDate(dt)

    return(gsta)

def strtodate(s):
    return(AbsoluteDate(DateTimeComponents.parseDateTime(String(s)),
                        TimeScalesFactory.getUTC()))

def pvtolist(okpv):
    return(okpv.getPosition().toArray() + okpv.getVelocity().toArray())

_itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)
