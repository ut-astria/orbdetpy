# test_estimation.py - Run measurement simulation and OD tests.
# Copyright (C) 2020-2022 University of Texas
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
from orbdetpy import (configure, add_facet, add_maneuver, add_station, build_measurement, AttitudeType,
                      Filter, ManeuverTrigger, ManeuverType, MeasurementType, Constant)
from orbdetpy.ccsds import export_OEM, export_TDM
from orbdetpy.conversion import get_J2000_epoch_offset, get_UTC_string
from orbdetpy.estimation import determine_orbit
from orbdetpy.propagation import propagate_orbits
from orbdetpy.plotting.estimation import plot as estimation_plot
from orbdetpy.plotting.simulation import plot as simulation_plot

# Set up configuration for measurement generation
cfg = configure(rso_mass=2500.0, prop_start=get_J2000_epoch_offset("2019-05-01T00:00:00"),
                prop_initial_state=[-23183898.259, 35170229.755, 43425.075, -2566.938, -1692.19, 138.948],
                prop_end=get_J2000_epoch_offset("2019-05-01T01:00:00"), prop_step=300.0, sim_measurements=True)

# Uncomment to define box-wing model
#add_facet(cfg, Constant.PLUS_I, 3.0)
#add_facet(cfg, Constant.PLUS_J, 3.0)
#add_facet(cfg, Constant.PLUS_K, 3.0)
#add_facet(cfg, Constant.MINUS_I, 3.0)
#add_facet(cfg, Constant.MINUS_J, 3.0)
#add_facet(cfg, Constant.MINUS_K, 3.0)
#cfg.rso_solar_array_area = 20.0
#cfg.rso_solar_array_axis[:] = Constant.PLUS_J

# Uncomment to define initial attitude
#cfg.rso_attitude_provider = AttitudeType.FIXED_RATE
#cfg.rso_spin_velocity[:] = [30*Constant.DEGREE_TO_RAD, 0, 0]
#cfg.rso_spin_acceleration[:] = Constant.ZERO_VECTOR

# Uncomment to add a maneuver
#add_maneuver(cfg, get_J2000_epoch_offset("2019-05-01T00:10:00"), ManeuverTrigger.DATE_TIME,
#             None, ManeuverType.CONSTANT_THRUST, [*Constant.MINUS_J, 600, 0.01, 2500])

# ManeuverType.*_CHANGE maneuvers take a single "delta" argument representing change in the
# corresponding element; values are in meters for distances and radians for angles
#add_maneuver(cfg, get_J2000_epoch_offset("2019-05-01T00:10:00"), ManeuverTrigger.DATE_TIME, None, ManeuverType.PERIGEE_CHANGE, [5E4])

# Define ground stations
add_station(cfg, "Maui", 0.3614, -2.7272, 3059.0)
add_station(cfg, "Millstone", 0.7438, -1.2652, 100.0, fov_azimuth=-2.75, fov_elevation=0.62, fov_aperture=15*Constant.DEGREE_TO_RAD)

# Uncomment mutually exclusive sections below to choose different measurements 
# AZIMUTH and ELEVATION must be given
#cfg.measurements[MeasurementType.AZIMUTH].error[:] = [Constant.ARC_SECOND_TO_RAD]
#cfg.measurements[MeasurementType.ELEVATION].error[:] = [Constant.ARC_SECOND_TO_RAD]

# You can provide RANGE+RANGE_RATE or RANGE alone or RANGE_RATE alone
#cfg.measurements[MeasurementType.RANGE].error[:] = [10.0]
#cfg.measurements[MeasurementType.RANGE_RATE].error[:] = [0.1]

# RIGHT_ASCENSION and DECLINATION must be given
cfg.measurements[MeasurementType.RIGHT_ASCENSION].error[:] = [Constant.ARC_SECOND_TO_RAD]
cfg.measurements[MeasurementType.DECLINATION].error[:] = [Constant.ARC_SECOND_TO_RAD]

# Inertial position as measurement
#cfg.measurements[MeasurementType.POSITION].error[:] = [100.0, 200.0, 300.0]

# Inertial position/velocity as measurement
#cfg.measurements[MeasurementType.POSITION_VELOCITY].error[:] = [100.0, 200.0, 300.0, 3.0, 2.0, 1.0]

# Propagate orbits and generate measurements
meas = list(m for m in propagate_orbits([cfg])[0].array)

# Empty measurements are ignored in orbit estimation but they force the filters
# to output propagated states/covariances at arbitrary intermediate times
t0 = get_J2000_epoch_offset("2019-05-01T00:26:02.853")
for i in range(5):
    meas.append(build_measurement(t0 + i*60.0, "", []))

# Measurements must be sorted in ascending chronological order
meas.sort(key=lambda m: m.time)

# Plot orbital elements from the simulation
#simulation_plot(meas, interactive=True)

# Export simulated measurements to a CCSDS TDM file
#with open("orbdetpy_sim.tdm", "w") as fp:
#    fp.write(export_TDM(cfg, meas, "DEBRIS"))

# Perturb truth initial state before running OD
cfg.prop_initial_state[:] = multivariate_normal(cfg.prop_initial_state, diag([1E6, 1E6, 1E6, 1E2, 1E2, 1E2]))
# Use Filter.EXTENDED_KALMAN for the EKF
cfg.estm_filter = Filter.UNSCENTED_KALMAN
# Set the initial covariance; diagonal entries for the estimated state vector and parameters 
cfg.estm_covariance[:] = [1E6, 1E6, 1E6, 1E2, 1E2, 1E2, 1.00, 0.25, 1E-6, 1E-6, 1E-6]

# Run OD on simulated measurements
fit = determine_orbit([cfg], [meas])[0]

for f in fit:
    # print(f) to dump pre-fits/post-fits/covariances
    print(get_UTC_string(f.time), f.station, f.estimated_state[:6])

# Plot OD results
estimation_plot(cfg, meas, fit, interactive=True, estim_param=False)

# Export ephemeris to a CCSDS OEM file
#with open("orbdetpy_fit.oem", "w") as fp:
#    fp.write(export_OEM(cfg, fit, "2021-001A", "DEBRIS"))
