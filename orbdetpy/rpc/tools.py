# tools.py - RPC utility functions.
# Copyright (C) 2019 University of Texas
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
    "rso_mass": [["SpaceObject", "Mass"], 1.0],
    "rso_area": [["SpaceObject", "Area"], 1.0],
    "rso_solar_array_axis": [["SpaceObject", "SolarArray", "Axis"], None],
    "rso_solar_array_area": [["SpaceObject", "SolarArray", "Area"], None],
    "rso_attitude_provider": [["SpaceObject", "Attitude", "Provider"], None],
    "rso_spin_velocity": [["SpaceObject", "Attitude", "SpinVelocity"], None],
    "rso_spin_acceleration": [["SpaceObject", "Attitude", "SpinAcceleration"], None],
    "gravity_degree": [["Gravity", "Degree"], 20],
    "gravity_order": [["Gravity", "Order"], 20],
    "ocean_tides_degree": [["OceanTides", "Degree"], 20],
    "ocean_tides_order": [["OceanTides", "Order"], 20],
    "third_body_sun": [["ThirdBodies", "Sun"], True],
    "third_body_moon": [["ThirdBodies", "Moon"], True],
    "solid_tides_sun": [["SolidTides", "Sun"], True],
    "solid_tides_moon": [["SolidTides", "Moon"], True],
    "drag_model": [["Drag", "Model"], "MSISE"],
    "drag_exp_rho0": [["Drag", "ExpRho0"], 0.0],
    "drag_exp_H0": [["Drag", "ExpH0"], 0.0],
    "drag_exp_Hscale": [["Drag", "ExpHScale"], 0.0],
    "rp_sun": [["RadiationPressure", "Sun"], True],
    "rp_coeff_absorption": [["RadiationPressure", "Cabsorption"], 0.1],
    "prop_start": [["Propagation", "Start"], None],
    "prop_end": [["Propagation", "End"], None],
    "prop_step": [["Propagation", "Step"], 0.0],
    "prop_initial_state": [["Propagation", "InitialState"], None],
    "prop_initial_TLE": [["Propagation", "InitialTLE"], None],
    "prop_inertial_frame": [["Propagation", "InertialFrame"], None],
    "prop_step_handler_start_time": [["Propagation", "StepHandlerStartTime"], None],
    "prop_step_handler_end_time": [["Propagation", "StepHandlerEndTime"], None],
    "integ_min_time_step": [["Integration", "MinTimeStep"], 1.0E-3],
    "integ_max_time_step": [["Integration", "MaxTimeStep"], 300.0],
    "integ_abs_tolerance": [["Integration", "AbsTolerance"], 1.0E-14],
    "integ_rel_tolerance": [["Integration", "RelTolerance"], 1.0E-12],
    "sim_measurements": [["Simulation", "SimulateMeasurements"], True],
    "sim_skip_unobservable": [["Simulation", "SkipUnobservable"], True],
    "sim_include_extras": [["Simulation", "IncludeExtras"], False],
    "sim_include_station_state": [["Simulation", "IncludeStationState"], False],
    "estm_filter": [["Estimation", "Filter"], None],
    "estm_covariance": [["Estimation", "Covariance"], None],
    "estm_process_noise": [["Estimation", "ProcessNoise"], None],
    "estm_DMC_corr_time": [["Estimation", "DMCCorrTime"], 0.0],
    "estm_DMC_sigma_pert": [["Estimation", "DMCSigmaPert"], 0.0],
    "estm_outlier_sigma": [["Estimation", "OutlierSigma"], 0.0],
    "estm_outlier_warmup": [["Estimation", "OutlierWarmup"], 0]
}

_measurement_fields = {
    "time": "Time",
    "station": "Station",
    "azimuth": "Azimuth",
    "elevation": "Elevation",
    "range": "Range",
    "range_rate": "RangeRate",
    "right_ascension": "RightAscension",
    "declination": "Declination",
    "position": "Position",
    "position_velocity": "PositionVelocity",
    "atmospheric_density": "AtmDensity",
    "acceleration_gravity": "AccGravity",
    "acceleration_drag": "AccDrag",
    "acceleration_ocean_tides": "AccOceanTides",
    "acceleration_solid_tides": "AccSolidTides",
    "acceleration_third_bodies": "AccThirdBodies",
    "acceleration_radiation_pressure": "AccRadiationPressure",
    "acceleration_thrust": "AccThrust",
    "station_state": "StationState",
    "true_state_cartesian": "TrueStateCartesian",
    "true_state_sma": "TrueStateSMA",
    "true_state_ecc": "TrueStateEcc",
    "true_state_inc": "TrueStateInc",
    "true_state_raan": "TrueStateRAAN",
    "true_state_argp": "TrueStateArgp",
    "true_state_mean_anom": "TrueStateMeanAnom",
    "true_state_ex": "TrueStateEx",
    "true_state_ey": "TrueStateEy",
    "true_state_hx": "TrueStateHx",
    "true_state_hy": "TrueStateHy",
    "true_state_lm": "TrueStateLm"
}

def build_settings(param):
    inp = read_param(param)
    cfg = messages_pb2.Settings()
    for dest, src in _settings_fields.items():
        fld = inp
        for s in src[0]:
            fld = fld.get(s, src[1]) if s == src[0][-1] else fld.get(s)
            if (fld is None):
                fld = src[1]
                break
        if (fld is None):
            continue
        if (isinstance(fld, list)):
            getattr(cfg, dest).extend(fld)
        else:
            setattr(cfg, dest, fld)

    for f in inp.get("SpaceObject", {}).get("Facets", []):
        fac = messages_pb2.Facet(area = f["Area"])
        fac.normal.extend(f["Normal"])
        cfg.rso_facets.append(fac)

    coef = inp.get("Drag", {}).get("Coefficient", {})
    cfg.drag_coefficient.CopyFrom(messages_pb2.Parameter(name="Cd", value=coef.get("Value", 2.0),
                                                         min=coef.get("Min", 1.0), max=coef.get("Max", 3.0),
                                                         estimation=coef.get("Estimation", "Estimate")))
    for f in inp.get("Drag", {}).get("MSISEFlags", []):
        obj = messages_pb2.IntegerArray()
        obj.array.extend(f)
        cfg.drag_MSISE_flags.append(obj)

    coef = inp.get("RadiationPressure", {}).get("Creflection", {})
    cfg.rp_coeff_reflection.CopyFrom(messages_pb2.Parameter(name="Cr", value=coef.get("Value", 1.5),
                                                            min=coef.get("Min", 1.0), max=coef.get("Max", 2.0),
                                                            estimation=coef.get("Estimation", "Estimate")))

    for m in inp.get("Maneuvers", []):
        man = messages_pb2.Maneuver(time = m["Time"], trigger_event = m["TriggerEvent"],
                                    maneuver_type = m["ManeuverType"])
        if ("TriggerParams" in m):
            man.trigger_params.extend(m["TriggerParams"])
        if ("ManeuverParams" in m):
            man.maneuver_params.extend(m["ManeuverParams"])
        cfg.maneuvers.append(man)

    for k, v in inp.get("Stations", {}).items():
        sta =  messages_pb2.Station(latitude=v["Latitude"], longitude=v["Longitude"],
                                    altitude=v["Altitude"], azimuth_bias=v.get("AzimuthBias", 0.0),
                                    elevation_bias=v.get("ElevationBias", 0.0),
                                    range_bias=v.get("RangeBias", 0.0),
                                    range_rate_bias=v.get("RangeRateBias", 0.0),
                                    right_ascension_bias=v.get("RightAscensionBias", 0.0),
                                    declination_bias=v.get("DeclinationBias", 0.0),
                                    bias_estimation=v.get("BiasEstimation", ""))
        sta.position_bias.extend(v.get("PositionBias", [0.0]*3))
        sta.position_velocity_bias.extend(v.get("PositionVelocityBias", [0.0]*6))
        cfg.stations[k].CopyFrom(sta)

    for k, v in inp.get("Measurements", {}).items():
        mea = messages_pb2.MeasurementSetting(two_way = v.get("TwoWay", True))
        mea.error.extend(v["Error"])
        cfg.measurements[k].CopyFrom(mea)

    coef = inp.get("Estimation", {}).get("DMCAcceleration", {})
    cfg.estm_DMC_acceleration.CopyFrom(messages_pb2.Parameter(name="", value=coef.get("Value", 0.0),
                                                              min=coef.get("Min", -1.0E-3), max=coef.get("Max", 1.0E-3),
                                                              estimation=coef.get("Estimation", "Estimate")))
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
    output = [[inp.time, [list(i.array) for i in inp.states]]
              for inp in inputs]
    return(output)

def convert_measurements(inputs):
    output = []
    for inp in inputs:
        out = {}
        output.append(out)
        for src, dest in _measurement_fields.items():
            fld = getattr(inp, src)
            if (fld):
                out[dest] = list(fld) if dest == "TrueStateCartesian" else fld

    return(output)

def convert_estimation(inputs):
    output = []
    for inp in inputs:
        out = {}
        output.append(out)
        out["Time"] = inp.time
        if (inp.station):
            out["Station"] = inp.station
        out["EstimatedState"] = list(inp.estimated_state)
        if (inp.propagated_covariance):
            out["PropagatedCovariance"] = [list(i.array) for i in inp.propagated_covariance]
        if (inp.innovation_covariance):
            out["InnovationCovariance"] = [list(i.array) for i in inp.innovation_covariance]
        if (inp.estimated_covariance):
            out["EstimatedCovariance"] = [list(i.array) for i in inp.estimated_covariance]
        if (inp.pre_fit):
            out["PreFit"] = {}
            for key, val in inp.pre_fit.items():
                out["PreFit"][key] = list(val.array)
        if (inp.post_fit):
            out["PostFit"] = {}
            for key, val in inp.post_fit.items():
                out["PostFit"][key] = list(val.array)

    return(output)
