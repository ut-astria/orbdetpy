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

from numpy import array
from orbdetpy import _datadir
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

    dragmod = config["Drag"]["Model"]
    if (dragmod == "Exponential"):
        atmos = SimpleExponentialAtmosphere(
            OneAxisEllipsoid(
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
            Constants.WGS84_EARTH_FLATTENING, _itrf),
            config["Drag"]["ExpRho0"], config["Drag"]["ExpH0"],
            config["Drag"]["ExpHScale"])
    elif (dragmod == "MSISE"):
        mindt, maxdt, swdata = spaceweather()
        msise = MSISEInputs(mindt, maxdt, swdata)

        atmos = NRLMSISE00(msise, CelestialBodyFactory.getSun(),
            OneAxisEllipsoid(
            Constants.WGS84_EARTH_EQUATORIAL_RADIUS,
                Constants.WGS84_EARTH_FLATTENING, _itrf)).withSwitch(9, -1)
        for s in config["Drag"].get("MSISEDisable", []):
            atmos = atmos.withSwitch(s, 0)
    else:
        atmos = None

    if (atmos is not None):
        fmod.append(DragForce(atmos,
            IsotropicDrag(config["SpaceObject"]["Area"],
            config["Drag"]["Coefficient"]["Value"])))

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
            IsotropicRadiationSingleCoefficient(
            config["SpaceObject"]["Area"],
            config["RadiationPressure"]["Creflection"]["Value"],
            config["RadiationPressure"]["Creflection"]["Min"],
            config["RadiationPressure"]["Creflection"]["Max"])))

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

def estparms(config):
    sdim = 6
    parm, pest = [], []
    grps = [["Drag", "Coefficient", DragSensitive.DRAG_COEFFICIENT],
            ["RadiationPressure", "Creflection",
             RadiationSensitive.REFLECTION_COEFFICIENT]]
    for g in grps:
        c = config.get(g[0], {}).get(g[1], {})
        if (c.get("Estimation", "") != "Estimate"):
            continue
        sdim += 1
        parm.append([c["Min"], c["Max"], c["Value"]])
        pest.append(g[2])

    return(sdim, array(parm), pest) 

def spaceweather():
    with open(os.path.join(_datadir, "SpaceWeather.dat"), "r") as fp:
        toks = []
        for l in fp:
            toks.append([t.strip() for t in l.split(",")])

    mindt = strtodate("%s-%s-%sT00:00:00Z" % (
        toks[0][0], toks[0][1], toks[0][2]))
    maxdt = strtodate("%s-%s-%sT00:00:00Z" % (
        toks[-1][0], toks[-1][1], toks[-1][2]))

    swdata = Hashtable()
    for tok in toks:
        val = [float(t) if len(t) > 0 else 0.0 for t in tok]
        swdata.put(String("%s%s%s" % (tok[0], tok[1], tok[2])), val)

    return(mindt, maxdt, swdata)

def strtodate(s):
    return(AbsoluteDate(DateTimeComponents.parseDateTime(String(s)),
                        TimeScalesFactory.getUTC()))

def pvtolist(okpv):
    return(okpv.getPosition().toArray() + okpv.getVelocity().toArray())

_itrf = FramesFactory.getITRF(IERSConventions.IERS_2010, False)
