# multi_target.py - Test multiple target OD.
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
from orbdetpy import configure, add_station, Constant, DragModel, EstimationType, MeasurementType
from orbdetpy.conversion import get_J2000_epoch_offset
from orbdetpy.estimation import multi_target_OD
from orbdetpy.propagation import propagate_orbits
from orbdetpy.plotting.estimation import plot
from orbdetpy.rpc.messages_pb2 import Parameter

noise = (1E-10,)*6
cov = (25E6, 25E6, 25E6, 1E2, 1E2, 1E2)
t0, t1 = get_J2000_epoch_offset(("2019-07-10T23:30:00", "2019-07-11T00:00:00"))

config = [configure(prop_start=t0, prop_end=t1, prop_step=60.0, prop_initial_state=(-20000151.1484, 21900.3443, 0.2039, 0, -4464.9556, 0),
                    drag_model=DragModel.UNDEFINED, ocean_tides_degree=-1, ocean_tides_order=-1, solid_tides_sun=False, solid_tides_moon=False,
                    drag_coefficient=Parameter(value=2.0, min=1.0, max=3.0, estimation=EstimationType.UNDEFINED),
                    rp_coeff_reflection=Parameter(value=1.5, min=1.0, max=2.0, estimation=EstimationType.UNDEFINED),
                    sim_measurements=True, estm_process_noise=noise, estm_covariance=cov),
          configure(prop_start=t0, prop_end=t1, prop_step=60.0, prop_initial_state=(-20000000, 0, 0, -100, -4464.3029, 0),
                    drag_model=DragModel.UNDEFINED, ocean_tides_degree=-1, ocean_tides_order=-1, solid_tides_sun=False, solid_tides_moon=False,
                    drag_coefficient=Parameter(value=2.0, min=1.0, max=3.0, estimation=EstimationType.UNDEFINED),
                    rp_coeff_reflection=Parameter(value=1.5, min=1.0, max=2.0, estimation=EstimationType.UNDEFINED),
                    sim_measurements=True, estm_process_noise=noise, estm_covariance=cov),
          configure(prop_start=t0, prop_end=t1, prop_step=60.0, prop_initial_state=(-20000000, -20000, -20000, -100, -4464.3029, 0),
                    drag_model=DragModel.UNDEFINED, ocean_tides_degree=-1, ocean_tides_order=-1, solid_tides_sun=False, solid_tides_moon=False,
                    drag_coefficient=Parameter(value=2.0, min=1.0, max=3.0, estimation=EstimationType.UNDEFINED),
                    rp_coeff_reflection=Parameter(value=1.5, min=1.0, max=2.0, estimation=EstimationType.UNDEFINED),
                    sim_measurements=True, estm_process_noise=noise, estm_covariance=cov)]
for cfg in config:
    add_station(cfg, "Arecibo", 0.3202, -1.1651, 497.0)
    add_station(cfg, "Kaena", 0.08683, -2.5872, 460.0)
    add_station(cfg, "Boston", 0.7495, -1.2504, 460.0)
    cfg.measurements[MeasurementType.RIGHT_ASCENSION].error[:] = [Constant.ARC_SECOND_TO_RAD]
    cfg.measurements[MeasurementType.DECLINATION].error[:] = [Constant.ARC_SECOND_TO_RAD]

obs_list = propagate_orbits(config)
merge_obs = [sorted((y for x in obs_list for y in list(x.array)), key=lambda x: f"{x.time}{x.station}")]

fit_data = multi_target_OD(config[:2], merge_obs)

print(f"Unassociated observations: {fit_data.unassociated_obs}")
for obs, fit in zip(fit_data.associated_obs, fit_data.est_output):
    assoc = [merge_obs[0][idx] for idx in obs.array]
    if (len(assoc) > 0):
        plot(config[0], assoc, fit.array, interactive=True, estim_param=False)
