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
from orbdetpy import configure, add_station, build_measurement, Constant, DragModel, MeasurementType
from orbdetpy.conversion import get_J2000_epoch_offset
from orbdetpy.estimation import multi_target_OD
from orbdetpy.propagation import propagate_orbits
from orbdetpy.plotting.estimation import plot as estimation_plot

t0, t1 = get_J2000_epoch_offset(("2019-07-10T23:30:00", "2019-07-11T00:00:00"))
config = [configure(prop_start=t0, prop_end=t1, prop_step=60.0, sim_measurements=True, drag_model=DragModel.UNDEFINED,
                    ocean_tides_degree=-1, ocean_tides_order=-1, solid_tides_sun=False, solid_tides_moon=False, 
                    estm_DMC_corr_time=0.0, estm_DMC_sigma_pert=0.0, estm_process_noise=[1E-10]*6,
                    estm_covariance=[25E6, 25E6, 25E6, 1E4, 1E4, 1E4, 1.00, 0.25],
                    prop_initial_state=(-20000151.1484369, 21900.344312363, 0.2039385729, 0, -4464.95561397887, 0)),
          configure(prop_start=t0, prop_end=t1, prop_step=60.0, sim_measurements=True, drag_model=DragModel.UNDEFINED,
                    ocean_tides_degree=-1, ocean_tides_order=-1, solid_tides_sun=False, solid_tides_moon=False,
                    estm_DMC_corr_time=0.0, estm_DMC_sigma_pert=0.0, estm_process_noise=[1E-10]*6,
                    estm_covariance=[25E6, 25E6, 25E6, 1E4, 1E4, 1E4, 1.00, 0.25],
                    prop_initial_state=(-20000000, 0, 0, -100, -4464.302857, 0)),
          configure(prop_start=t0, prop_end=t1, prop_step=60.0, sim_measurements=True, drag_model=DragModel.UNDEFINED,
                    ocean_tides_degree=-1, ocean_tides_order=-1, solid_tides_sun=False, solid_tides_moon=False,
                    estm_DMC_corr_time=0.0, estm_DMC_sigma_pert=0.0, estm_process_noise=[1E-10]*6,
                    estm_covariance=[25E6, 25E6, 25E6, 1E4, 1E4, 1E4, 1.00, 0.25],
                    prop_initial_state=(-20000000, -20000, -20000, -100, -4464.302857, 0))]

for cfg in config:
    add_station(cfg, "Arecibo", 0.32016610686, -1.16505575707, 497.0)
    add_station(cfg, "Kaena", 0.0868343901, -2.58718978041, 460.0)
    add_station(cfg, "Boston", 0.74949349786, -1.25044863878, 460.0)
    cfg.measurements[MeasurementType.RIGHT_ASCENSION].error[:] = [Constant.ARC_SECOND_TO_RAD]
    cfg.measurements[MeasurementType.DECLINATION].error[:] = [Constant.ARC_SECOND_TO_RAD]

obs_list = propagate_orbits(config)
merge_obs = [sorted((y for x in obs_list for y in list(x.array)), key=lambda x: f"{x.time}{x.station}")]

fit_data = multi_target_OD(config[:2], merge_obs)

print(fit_data.unassociated_obs)
for cfg, obs, fit in zip(config, fit_data.associated_obs, fit_data.est_output):
    assoc = [merge_obs[0][idx] for idx in obs.array]
    print()
    print(obs.array)
    print(assoc[-1].true_state)
    print(fit.array[-1].estimated_state[:6])
#    estimation_plot(cfg, assoc, fit.array, interactive=True)
