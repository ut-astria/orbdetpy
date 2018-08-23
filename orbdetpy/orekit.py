# orekit.py - Module to instantiate Orekit classes using Pyjnius.
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

import os
from jnius import autoclass
from orbdetpy import _datadir

# Java classes
File = autoclass("java.io.File")
String = autoclass("java.lang.String")
ArrayList = autoclass("java.util.ArrayList")

# Hipparchus classes
Vector3D = autoclass("org.hipparchus.geometry.euclidean.threed.Vector3D")
Array2DRowRealMatrix = autoclass("org.hipparchus.linear.Array2DRowRealMatrix")
DiagonalMatrix = autoclass("org.hipparchus.linear.DiagonalMatrix")
DormandPrince853Integrator = autoclass("org.hipparchus.ode.nonstiff.DormandPrince853Integrator")

# Orekit classes
CelestialBodyFactory = autoclass("org.orekit.bodies.CelestialBodyFactory")
GeodeticPoint = autoclass("org.orekit.bodies.GeodeticPoint")
OneAxisEllipsoid = autoclass("org.orekit.bodies.OneAxisEllipsoid")

DirectoryCrawler = autoclass("org.orekit.data.DirectoryCrawler")
DataProvidersManager = autoclass("org.orekit.data.DataProvidersManager")

AngularAzEl = autoclass("org.orekit.estimation.measurements.AngularAzEl")
AngularRaDec = autoclass("org.orekit.estimation.measurements.AngularRaDec")
GroundStation = autoclass("org.orekit.estimation.measurements.GroundStation")
Range = autoclass("org.orekit.estimation.measurements.Range")
RangeRate = autoclass("org.orekit.estimation.measurements.RangeRate")
ConstantProcessNoise = autoclass("org.orekit.estimation.sequential.ConstantProcessNoise")
KalmanEstimatorBuilder = autoclass("org.orekit.estimation.sequential.KalmanEstimatorBuilder")

DragForce = autoclass("org.orekit.forces.drag.DragForce")
DragSensitive = autoclass("org.orekit.forces.drag.DragSensitive")
IsotropicDrag = autoclass("org.orekit.forces.drag.IsotropicDrag")
SimpleExponentialAtmosphere = autoclass("org.orekit.forces.drag.atmosphere.SimpleExponentialAtmosphere")
GravityFieldFactory = autoclass("org.orekit.forces.gravity.potential.GravityFieldFactory")
HolmesFeatherstoneAttractionModel = autoclass("org.orekit.forces.gravity.HolmesFeatherstoneAttractionModel")
OceanTides = autoclass("org.orekit.forces.gravity.OceanTides")
SolidTides = autoclass("org.orekit.forces.gravity.SolidTides")
NewtonianAttraction = autoclass("org.orekit.forces.gravity.NewtonianAttraction")
ThirdBodyAttraction = autoclass("org.orekit.forces.gravity.ThirdBodyAttraction")
ConstantThrustManeuver = autoclass("org.orekit.forces.maneuvers.ConstantThrustManeuver")
IsotropicRadiationClassicalConvention = autoclass("org.orekit.forces.radiation.IsotropicRadiationClassicalConvention")
RadiationSensitive = autoclass("org.orekit.forces.radiation.RadiationSensitive")
SolarRadiationPressure = autoclass("org.orekit.forces.radiation.SolarRadiationPressure")

Frame = autoclass("org.orekit.frames.Frame")
FramesFactory = autoclass("org.orekit.frames.FramesFactory")
TopocentricFrame = autoclass("org.orekit.frames.TopocentricFrame")
Transform = autoclass("org.orekit.frames.Transform")

CartesianOrbit = autoclass("org.orekit.orbits.CartesianOrbit")
PositionAngle = autoclass("org.orekit.orbits.PositionAngle")

DormandPrince853IntegratorBuilder = autoclass("org.orekit.propagation.conversion.DormandPrince853IntegratorBuilder")
NumericalPropagatorBuilder = autoclass("org.orekit.propagation.conversion.NumericalPropagatorBuilder")
NumericalPropagator = autoclass("org.orekit.propagation.numerical.NumericalPropagator")
SpacecraftState = autoclass("org.orekit.propagation.SpacecraftState")

AbsoluteDate = autoclass("org.orekit.time.AbsoluteDate")
DateTimeComponents = autoclass("org.orekit.time.DateTimeComponents")
TimeScalesFactory = autoclass("org.orekit.time.TimeScalesFactory")

Constants = autoclass("org.orekit.utils.Constants")
IERSConventions = autoclass("org.orekit.utils.IERSConventions")
ParameterDriver = autoclass("org.orekit.utils.ParameterDriver")
ParameterDriversList = autoclass("org.orekit.utils.ParameterDriversList")
PVCoordinates = autoclass("org.orekit.utils.PVCoordinates")
TimeStampedPVCoordinates = autoclass("org.orekit.utils.TimeStampedPVCoordinates")

# Orbdetpy Java utilities class
PropUtil = autoclass("org.astria.PropUtil")

DataProvidersManager.getInstance().addProvider(DirectoryCrawler(File(String(_datadir))))
