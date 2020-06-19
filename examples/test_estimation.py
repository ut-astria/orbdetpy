# test_estimation.py - Run measurement simulation and OD tests.
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

from numpy import diag
from numpy.random import multivariate_normal
from orbdetpy import (configure, add_facet, add_maneuver, add_station,
                      AttitudeType, Filter, ManeuverTrigger, ManeuverType,
                      MeasurementType, Constant)
from orbdetpy.conversion import get_J2000_epoch_offset
from orbdetpy.estimation import determine_orbit
from orbdetpy.propagation import propagate_orbits
from orbdetpy.plotting.estimation import plot

# Set up configuration for measurement generation
cfg = configure(prop_start=get_J2000_epoch_offset("2019-05-01T00:00:00Z"),
                prop_initial_state=[-23183898.259, 35170229.755, 43425.075,
                                    -2566.938, -1692.19, 138.948],
                prop_end=get_J2000_epoch_offset("2019-05-01T01:00:00Z"),
                prop_step=300.0,
                sim_measurements=True)

# Uncomment code block below to define box-wing model
#add_facet(cfg, Constant.PLUS_I, 3.0)
#add_facet(cfg, Constant.PLUS_J, 3.0)
#add_facet(cfg, Constant.PLUS_K, 3.0)
#add_facet(cfg, Constant.MINUS_I, 3.0)
#add_facet(cfg, Constant.MINUS_J, 3.0)
#add_facet(cfg, Constant.MINUS_K, 3.0)
#cfg.rso_solar_array_area = 20.0
#cfg.rso_solar_array_axis[:] = Constant.PLUS_J

# Uncomment code block below to define initial attitude
#cfg.rso_attitude_provider = AttitudeType.FIXED_RATE
#cfg.rso_spin_velocity[:] = [30*Constant.DEGREE, 0, 0]
#cfg.rso_spin_acceleration[:] = Constant.ZERO

# Uncomment code block below to add a maneuver
#add_maneuver(cfg, get_J2000_epoch_offset("2019-05-01T00:10:00Z"), ManeuverTrigger.DATE_TIME,
#             None, ManeuverType.CONSTANT_THRUST, [*Constant.MINUS_J, 30, 100, 250])

# Define ground stations
add_station(cfg, "Maui", 0.3614, -2.7272, 3059.0)
add_station(cfg, "Millstone", 0.7438, -1.2652, 100.0)

# Uncomment sections below to choose measurement types
#cfg.measurements[MeasurementType.AZIMUTH].error[:] = [Constant.ARC_SECOND]
#cfg.measurements[MeasurementType.ELEVATION].error[:] = [Constant.ARC_SECOND]

#cfg.measurements[MeasurementType.RANGE].error[:] = [10.0]
#cfg.measurements[MeasurementType.RANGE_RATE].error[:] = [0.1]

cfg.measurements[MeasurementType.RIGHT_ASCENSION].error[:] = [Constant.ARC_SECOND]
cfg.measurements[MeasurementType.DECLINATION].error[:] = [Constant.ARC_SECOND]

# Propagate orbits and generate measurements
meas = propagate_orbits([cfg])[0].array

# Perturb truth initial state before running OD
cfg.prop_initial_state[:] = multivariate_normal(
    cfg.prop_initial_state, diag([1E6, 1E6, 1E6, 1E2, 1E2, 1E2]))
cfg.prop_step = 0.0
cfg.estm_filter = Filter.EXTENDED_KALMAN

# Run OD on simulated measurements and plot residuals
fit = determine_orbit([cfg], [meas])[0]
plot(cfg, meas, fit, interactive=True, estim_param=False)
