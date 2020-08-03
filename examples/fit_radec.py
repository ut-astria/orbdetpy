# fit_radec.py - Run OD on real angles data.
# Copyright (C) 2020 University of Texas
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

from orbdetpy import configure, add_station, build_measurement, Filter, Frame, MeasurementType, Constant
from orbdetpy.conversion import get_J2000_epoch_offset, get_UTC_string
from orbdetpy.estimation import determine_orbit, iod_laplace
from orbdetpy.plotting.estimation import plot

# "UTA-ASTRIA-NMSkies" ground station WGS-84 coordinates
lat, lon, alt = 0.5743, -1.8418, 2225.0

# Real measurements: UTC, RA (rad), dec (rad)
real_obs = [
    ["2020-07-24T03:19:36.035Z", 1.073579, 1.381180],
    ["2020-07-24T03:19:40.551Z", 0.950022, 1.386485],
    ["2020-07-24T03:19:45.226Z", 0.813143, 1.388843],
    ["2020-07-24T03:20:08.080Z", 0.177733, 1.348062],
    ["2020-07-24T03:20:12.593Z", 0.083064, 1.330833],
    ["2020-07-24T03:20:17.135Z", 6.282924, 1.311199],
    ["2020-07-24T03:20:35.638Z", 6.043317, 1.213738],
    ["2020-07-24T03:20:40.159Z", 6.002756, 1.186840],
    ["2020-07-24T03:20:44.824Z", 5.966159, 1.158158],
    ["2020-07-24T03:20:49.538Z", 5.933827, 1.128431],
    ["2020-07-24T03:20:54.236Z", 5.905503, 1.098115],
    ["2020-07-24T03:20:58.770Z", 5.881280, 1.068338],
    ["2020-07-24T03:21:07.908Z", 5.840156, 1.007179],
    ["2020-07-24T03:21:12.645Z", 5.822224, 0.975141],
    ["2020-07-24T03:21:17.240Z", 5.806584, 0.943872],
    ["2020-07-24T03:21:26.480Z", 5.779668, 0.880952],
    ["2020-07-24T03:21:31.131Z", 5.768043, 0.849405],
    ["2020-07-24T03:21:35.713Z", 5.757657, 0.818491],
    ["2020-07-24T03:21:40.249Z", 5.748296, 0.788056],
    ["2020-07-24T03:21:44.927Z", 5.739520, 0.757038],
    ["2020-07-24T03:21:54.265Z", 5.724259, 0.696125],
    ["2020-07-24T03:21:59.018Z", 5.717508, 0.665763],
    ["2020-07-24T03:22:03.679Z", 5.711479, 0.636510],
    ["2020-07-24T03:22:08.283Z", 5.706017, 0.608080],
    ["2020-07-24T03:22:12.904Z", 5.701007, 0.580095],
    ["2020-07-24T03:22:17.432Z", 5.696508, 0.553190],
    ["2020-07-24T03:22:22.073Z", 5.692283, 0.526236],
    ["2020-07-24T03:22:26.745Z", 5.688394, 0.499645],
    ["2020-07-24T03:22:31.417Z", 5.684841, 0.473698],
    ["2020-07-24T03:22:35.997Z", 5.681670, 0.448829],
    ["2020-07-24T03:22:40.696Z", 5.678703, 0.423978],
    ["2020-07-24T03:22:45.390Z", 5.675990, 0.399725],
    ["2020-07-24T03:22:49.930Z", 5.673622, 0.376921]
]

# Use Laplace IOD to estimate an initial state from 3 RA/dec
times, ra, dec = [], [], []
for i in range(3):
    times.append(get_J2000_epoch_offset(real_obs[i][0]))
    ra.append(real_obs[i][1])
    dec.append(real_obs[i][2])
initial_epoch = times[1]
initial_state = iod_laplace(Frame.EME2000, lat, lon, alt, times, ra, dec)

# Configure orbit determination
config = configure(prop_start=initial_epoch, prop_initial_state=initial_state,
                   prop_end=get_J2000_epoch_offset(real_obs[-1][0]))

# Define ground station(s)
add_station(config, "UTA-ASTRIA-NMSkies", lat, lon, alt)

# Define measurement types; RA/dec used here
config.measurements[MeasurementType.RIGHT_ASCENSION].error[:] = [2.0*Constant.ARC_SECOND_TO_RAD]
config.measurements[MeasurementType.DECLINATION].error[:] = [2.0*Constant.ARC_SECOND_TO_RAD]

# Use Filter.EXTENDED_KALMAN for the EKF
config.estm_filter = Filter.UNSCENTED_KALMAN

# Build a list of Measurement objects, one for each RA/dec pair
meas_obj = []
for o in real_obs:
    meas_obj.append(build_measurement(get_J2000_epoch_offset(o[0]), "UTA-ASTRIA-NMSkies", o[1:]))

# Run OD. Fitting single object and hence the [0]
fit = determine_orbit([config], [meas_obj])[0]
# Check for estimation errors
if (isinstance(fit, str)):
    print(fit)
    exit(1)

for f in fit:
    # print(f) to dump pre-fits/post-fits/covariances
    print(get_UTC_string(f.time), f.station, f.estimated_state[:6])
# Plot OD results
plot(config, meas_obj, fit, interactive=True, estim_param=False)
