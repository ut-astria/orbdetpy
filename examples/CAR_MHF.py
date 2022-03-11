# CAR_MHF.py - Test CAR-MHF with real telescope measurements.
# Copyright (C) 2021-2022 University of Texas
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

import json
from orbdetpy import configure, add_station, build_measurement, DragModel, Frame, MeasurementType, Constant
from orbdetpy.conversion import get_J2000_epoch_offset
from orbdetpy.estimation import multi_target_OD
import os

config = [configure(prop_start=get_J2000_epoch_offset("2020-09-16T08:27:22.099Z"), rso_mass=2000.0, rso_area=20.0, drag_model=DragModel.UNDEFINED,
                    ocean_tides_degree=-1, ocean_tides_order=-1, solid_tides_sun=False, solid_tides_moon=False, estm_DMC_corr_time=0.0,
                    estm_DMC_sigma_pert=0.0, estm_process_noise=(1E-8,)*6, estm_covariance=[25E6, 25E6, 25E6, 1E4, 1E4, 1E4, 1.00, 0.25],
                    prop_initial_state=(1.7307E7, -2745771.4791, 2.0497E7, -896.8248, 3482.0158, 1274.9837))]
config[0].measurements[MeasurementType.RIGHT_ASCENSION].error[:] = [Constant.ARC_SECOND_TO_RAD]
config[0].measurements[MeasurementType.DECLINATION].error[:] = [Constant.ARC_SECOND_TO_RAD]

add_station(config[0], "NMSkies", 0.57426665348, -1.84183820339, 2225.04)

meas_obj = [[]]
with open(os.path.join(os.path.dirname(__file__), "CAR_MHF_obs.json"), "r") as fp:
    for m in json.load(fp):
        meas_obj[0].append(build_measurement(get_J2000_epoch_offset(m[0]), "NMSkies", m[1:3]))
        if (len(m) > 3):
           meas_obj[0][-1].angle_rates[:] = m[3:]

fit_data = multi_target_OD(config, meas_obj)
