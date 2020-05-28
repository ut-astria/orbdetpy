# tools.py - RPC utility functions.
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

from orbdetpy import read_param
from orbdetpy.rpc import messages_pb2

_settings_fields = {
    "rso_mass": ["rsoMass", 5.0],
    "rso_area": ["rsoArea", 0.1],
    "rso_solar_array_axis": ["rsoSolarArrayAxis", None],
    "rso_solar_array_area": ["rsoSolarArrayArea", None],
    "rso_attitude_provider": ["rsoAttitudeProvider", None],
    "rso_spin_velocity": ["rsoSpinVelocity", None],
    "rso_spin_acceleration": ["rsoSpinAcceleration", None],
    "gravity_degree": ["gravityDegree", 20],
    "gravity_order": ["gravityOrder", 20],
    "ocean_tides_degree": ["oceanTidesDegree", 20],
    "ocean_tides_order": ["oceanTidesOrder", 20],
    "third_body_sun": ["thirdBodySun", True],
    "third_body_moon": ["thirdBodyMoon", True],
    "solid_tides_sun": ["solidTidesSun", True],
    "solid_tides_moon": ["solidTidesMoon", True],
    "drag_model": ["dragModel", "MSISE"],
    "drag_exp_rho0": ["dragExpRho0", 3.614E-13],
    "drag_exp_H0": ["dragExpH0", 700000.0],
    "drag_exp_Hscale": ["dragExpHscale", 88667.0],
    "rp_sun": ["rpSun", True],
    "rp_coeff_absorption": ["rpCoeffAbsorption", 0.0],
    "prop_start": ["propStart", None],
    "prop_end": ["propEnd", None],
    "prop_step": ["propStep", 0.0],
    "prop_initial_state": ["propInitialState", None],
    "prop_initial_TLE": ["propInitialTLE", None],
    "prop_inertial_frame": ["propInertialFrame", "EME2000"],
    "prop_step_handler_start_time": ["propStepHandlerStartTime", None],
    "prop_step_handler_end_time": ["propStepHandlerEndTime", None],
    "integ_min_time_step": ["integMinTimeStep", 1.0E-3],
    "integ_max_time_step": ["integMaxTimeStep", 300.0],
    "integ_abs_tolerance": ["integAbsTolerance", 1.0E-14],
    "integ_rel_tolerance": ["integRelTolerance", 1.0E-12],
    "sim_measurements": ["simMeasurements", True],
    "sim_skip_unobservable": ["simSkipUnobservable", True],
    "sim_include_extras": ["simIncludeExtras", False],
    "sim_include_station_state": ["simIncludeStationState", False],
    "sim_include_angle_rates": ["simIncludeAngleRates", False],
    "estm_filter": ["estmFilter", None],
    "estm_covariance": ["estmCovariance", None],
    "estm_process_noise": ["estmProcessNoise", None],
    "estm_DMC_corr_time": ["estmDMCCorrTime", 0.0],
    "estm_DMC_sigma_pert": ["estmDMCSigmaPert", 0.0],
    "estm_outlier_sigma": ["estmOutlierSigma", 0.0],
    "estm_outlier_warmup": ["estmOutlierWarmup", 0],
    "estm_smoother_iterations": ["estmSmootherIterations", 10],
    "estm_enable_PDAF": ["estmEnablePDAF", False],
    "estm_enable_CAR_MHF": ["estmEnableCARMHF", False],
    "estm_detection_probability": ["estmDetectionProbability", 0.99],
    "estm_gating_probability": ["estmGatingProbability", 0.99],
    "estm_gating_threshold": ["estmGatingThreshold", 5.0]
}

_measurement_fields = {
    "time": "time",
    "station": "station",
    "azimuth": "azimuth",
    "elevation": "elevation",
    "range": "range",
    "range_rate": "rangeRate",
    "right_ascension": "rightAscension",
    "declination": "declination",
    "position": "position",
    "position_velocity": "positionVelocity",
    "angle_rates": "angleRates",
    "atmospheric_density": "atmDensity",
    "acceleration_gravity": "accGravity",
    "acceleration_drag": "accDrag",
    "acceleration_ocean_tides": "accOceanTides",
    "acceleration_solid_tides": "accSolidTides",
    "acceleration_third_bodies": "accThirdBodies",
    "acceleration_radiation_pressure": "accRadiationPressure",
    "acceleration_thrust": "accThrust",
    "station_state": "stationState",
    "true_state_cartesian": "trueStateCartesian",
    "true_state_sma": "trueStateSma",
    "true_state_ecc": "trueStateEcc",
    "true_state_inc": "trueStateInc",
    "true_state_raan": "trueStateRaan",
    "true_state_argp": "trueStateArgp",
    "true_state_mean_anom": "trueStateMeanAnom",
    "true_state_ex": "trueStateEx",
    "true_state_ey": "trueStateEy",
    "true_state_hx": "trueStateHx",
    "true_state_hy": "trueStateHy",
    "true_state_lm": "trueStateLm"
}

def build_settings(param):
    inp = read_param(param)
    cfg = messages_pb2.Settings()
    for dest, src in _settings_fields.items():
        fld = inp.get(src[0], src[1])
        if (fld is not None):
            if (isinstance(fld, list)):
                getattr(cfg, dest).extend(fld)
            else:
                setattr(cfg, dest, fld)

    for f in inp.get("rsoFacets", []):
        fac = messages_pb2.Facet(area = f["area"])
        fac.normal.extend(f["normal"])
        cfg.rso_facets.append(fac)

    coef = inp.get("dragCoefficient", {})
    cfg.drag_coefficient.CopyFrom(messages_pb2.Parameter(name="Cd", value=coef.get("value", 2.0),
                                                         min=coef.get("min", 1.0), max=coef.get("max", 3.0),
                                                         estimation=coef.get("estimation", "Estimate")))
    for f in inp.get("dragMSISEFlags", []):
        obj = messages_pb2.IntegerArray()
        obj.array.extend(f)
        cfg.drag_MSISE_flags.append(obj)

    coef = inp.get("rpCoeffReflection", {})
    cfg.rp_coeff_reflection.CopyFrom(messages_pb2.Parameter(name="Cr", value=coef.get("value", 1.5),
                                                            min=coef.get("min", 1.0), max=coef.get("max", 2.0),
                                                            estimation=coef.get("estimation", "Estimate")))

    for m in inp.get("cfgManeuvers", []):
        man = messages_pb2.Maneuver(time = m["time"], trigger_event = m["triggerEvent"],
                                    maneuver_type = m["maneuverType"])
        if ("triggerParams" in m):
            man.trigger_params.extend(m["triggerParams"])
        if ("maneuverParams" in m):
            man.maneuver_params.extend(m["maneuverParams"])
        cfg.maneuvers.append(man)

    for k, v in inp.get("cfgStations", {}).items():
        sta =  messages_pb2.Station(latitude=v["latitude"], longitude=v["longitude"],
                                    altitude=v["altitude"], azimuth_bias=v.get("azimuthBias", 0.0),
                                    elevation_bias=v.get("elevationBias", 0.0),
                                    range_bias=v.get("rangeBias", 0.0),
                                    range_rate_bias=v.get("rangeRateBias", 0.0),
                                    right_ascension_bias=v.get("rightAscensionBias", 0.0),
                                    declination_bias=v.get("declinationBias", 0.0),
                                    bias_estimation=v.get("biasEstimation", ""))
        sta.position_bias.extend(v.get("positionBias", [0.0]*3))
        sta.position_velocity_bias.extend(v.get("positionVelocityBias", [0.0]*6))
        cfg.stations[k].CopyFrom(sta)

    for k, v in inp.get("cfgMeasurements", {}).items():
        mea = messages_pb2.MeasurementSetting(two_way = v.get("twoWay", True))
        mea.error.extend(v["error"])
        cfg.measurements[k].CopyFrom(mea)

    coef = inp.get("estmDMCAcceleration", {})
    cfg.estm_DMC_acceleration.CopyFrom(messages_pb2.Parameter(name="", value=coef.get("value", 0.0),
                                                              min=coef.get("min", -1.0E-3), max=coef.get("max", 1.0E-3),
                                                              estimation=coef.get("estimation", "Estimate")))
    return(cfg)

def build_measurements(param):
    output = []
    inputs = read_param(param)
    for inp in inputs:
        obj = messages_pb2.Measurement()
        output.append(obj)
        for dest, src in _measurement_fields.items():
            fld = inp.get(src)
            if (fld):
                if (isinstance(fld, list)):
                    getattr(obj, dest).extend(fld)
                else:
                    setattr(obj, dest, fld)

    return(output)

def convert_propagation(inputs):
    output = [[inp.time, [list(i.array) for i in inp.states]] for inp in inputs]
    return(output)

def convert_measurements(inputs):
    output = []
    for inp in inputs:
        out = {}
        output.append(out)
        for src, dest in _measurement_fields.items():
            fld = getattr(inp, src)
            if (fld):
                out[dest] = list(fld) if (dest in ["stationState", "trueStateCartesian"]) else fld

    return(output)

def convert_estimation(inputs):
    output = []
    for inp in inputs:
        out = {}
        output.append(out)
        out["time"] = inp.time
        if (inp.station):
            out["station"] = inp.station
        out["estimatedState"] = list(inp.estimated_state)
        if (inp.propagated_covariance):
            out["propagatedCovariance"] = [list(i.array) for i in inp.propagated_covariance]
        if (inp.innovation_covariance):
            out["innovationCovariance"] = [list(i.array) for i in inp.innovation_covariance]
        if (inp.estimated_covariance):
            out["estimatedCovariance"] = [list(i.array) for i in inp.estimated_covariance]
        if (inp.pre_fit):
            out["preFit"] = {}
            for key, val in inp.pre_fit.items():
                out["preFit"][key] = list(val.array)
        if (inp.post_fit):
            out["postFit"] = {}
            for key, val in inp.post_fit.items():
                out["postFit"][key] = list(val.array)
        if (inp.clutter_probability):
            out["clutterProbability"] = inp.clutter_probability

    return(output)
