# __init__.py - orbdetpy package initialization.
# Copyright (C) 2018-2020 University of Texas
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

"""
.. include:: ../README.md
"""

from os import path
from math import pi, sqrt
from typing import List, Optional, Tuple
from .version import __version__
from orbdetpy.rpc.messages_pb2 import Facet, Maneuver, Measurement, Parameter, Settings, Station
from orbdetpy.rpc.server import RemoteServer

class AttitudeType():
    """Orekit attitude providers.
    """

    UNDEFINED = 0
    NADIR_POINTING = 1
    BODY_CENTER_POINTING = 2
    FIXED_RATE = 3

class DragModel():
    """Atmospheric drag density models.
    """

    UNDEFINED = 0
    EXPONENTIAL = 1
    MSISE2000 = 2
    WAM = 3

class Epoch():
    """Orekit AbsoluteDate epochs.
    """

    BEIDOU_EPOCH = 0
    CCSDS_EPOCH = 1
    FIFTIES_EPOCH = 2
    GALILEO_EPOCH = 3
    GLONASS_EPOCH = 4
    GPS_EPOCH = 5
    IRNSS_EPOCH = 6
    J2000_EPOCH = 7
    UNIX_EPOCH = 8
    JULIAN_EPOCH = 9
    MODIFIED_JULIAN_EPOCH = 10
    QZSS_EPOCH = 11

class EstimationType():
    """Parameter estimation types.
    """

    UNDEFINED = 0
    CONSIDER = 1
    ESTIMATE = 2

class Filter():
    """Estimation filters.
    """

    EXTENDED_KALMAN = 0
    UNSCENTED_KALMAN = 1

class Frame():
    """Orekit reference frames.
    """

    CIRF_CONVENTIONS_1996_ACCURATE_EOP = "CIRF_CONVENTIONS_1996_ACCURATE_EOP"
    CIRF_CONVENTIONS_1996_SIMPLE_EOP = "CIRF_CONVENTIONS_1996_SIMPLE_EOP"
    CIRF_CONVENTIONS_2003_ACCURATE_EOP = "CIRF_CONVENTIONS_2003_ACCURATE_EOP"
    CIRF_CONVENTIONS_2003_SIMPLE_EOP = "CIRF_CONVENTIONS_2003_SIMPLE_EOP"
    CIRF_CONVENTIONS_2010_ACCURATE_EOP = "CIRF_CONVENTIONS_2010_ACCURATE_EOP"
    CIRF_CONVENTIONS_2010_SIMPLE_EOP = "CIRF_CONVENTIONS_2010_SIMPLE_EOP"
    ECLIPTIC_CONVENTIONS_1996 = "ECLIPTIC_CONVENTIONS_1996"
    ECLIPTIC_CONVENTIONS_2003 = "ECLIPTIC_CONVENTIONS_2003"
    ECLIPTIC_CONVENTIONS_2010 = "ECLIPTIC_CONVENTIONS_2010"
    EME2000 = "EME2000"
    GCRF = "GCRF"
    GTOD_CONVENTIONS_1996_ACCURATE_EOP = "GTOD_CONVENTIONS_1996_ACCURATE_EOP"
    GTOD_CONVENTIONS_1996_SIMPLE_EOP = "GTOD_CONVENTIONS_1996_SIMPLE_EOP"
    GTOD_CONVENTIONS_2003_ACCURATE_EOP = "GTOD_CONVENTIONS_2003_ACCURATE_EOP"
    GTOD_CONVENTIONS_2003_SIMPLE_EOP = "GTOD_CONVENTIONS_2003_SIMPLE_EOP"
    GTOD_CONVENTIONS_2010_ACCURATE_EOP = "GTOD_CONVENTIONS_2010_ACCURATE_EOP"
    GTOD_CONVENTIONS_2010_SIMPLE_EOP = "GTOD_CONVENTIONS_2010_SIMPLE_EOP"
    GTOD_WITHOUT_EOP_CORRECTIONS = "GTOD_WITHOUT_EOP_CORRECTIONS"
    ICRF = "ICRF"
    ITRF_CIO_CONV_1996_ACCURATE_EOP = "ITRF_CIO_CONV_1996_ACCURATE_EOP"
    ITRF_CIO_CONV_1996_SIMPLE_EOP = "ITRF_CIO_CONV_1996_SIMPLE_EOP"
    ITRF_CIO_CONV_2003_ACCURATE_EOP = "ITRF_CIO_CONV_2003_ACCURATE_EOP"
    ITRF_CIO_CONV_2003_SIMPLE_EOP = "ITRF_CIO_CONV_2003_SIMPLE_EOP"
    ITRF_CIO_CONV_2010_ACCURATE_EOP = "ITRF_CIO_CONV_2010_ACCURATE_EOP"
    ITRF_CIO_CONV_2010_SIMPLE_EOP = "ITRF_CIO_CONV_2010_SIMPLE_EOP"
    ITRF_EQUINOX_CONV_1996_ACCURATE_EOP = "ITRF_EQUINOX_CONV_1996_ACCURATE_EOP"
    ITRF_EQUINOX_CONV_1996_SIMPLE_EOP = "ITRF_EQUINOX_CONV_1996_SIMPLE_EOP"
    ITRF_EQUINOX_CONV_2003_ACCURATE_EOP = "ITRF_EQUINOX_CONV_2003_ACCURATE_EOP"
    ITRF_EQUINOX_CONV_2003_SIMPLE_EOP = "ITRF_EQUINOX_CONV_2003_SIMPLE_EOP"
    ITRF_EQUINOX_CONV_2010_ACCURATE_EOP = "ITRF_EQUINOX_CONV_2010_ACCURATE_EOP"
    ITRF_EQUINOX_CONV_2010_SIMPLE_EOP = "ITRF_EQUINOX_CONV_2010_SIMPLE_EOP"
    MOD_CONVENTIONS_1996 = "MOD_CONVENTIONS_1996"
    MOD_CONVENTIONS_2003 = "MOD_CONVENTIONS_2003"
    MOD_CONVENTIONS_2010 = "MOD_CONVENTIONS_2010"
    MOD_WITHOUT_EOP_CORRECTIONS = "MOD_WITHOUT_EOP_CORRECTIONS"
    PZ90_11 = "PZ90_11"
    TEME = "TEME"
    TIRF_CONVENTIONS_1996_ACCURATE_EOP = "TIRF_CONVENTIONS_1996_ACCURATE_EOP"
    TIRF_CONVENTIONS_1996_SIMPLE_EOP = "TIRF_CONVENTIONS_1996_SIMPLE_EOP"
    TIRF_CONVENTIONS_2003_ACCURATE_EOP = "TIRF_CONVENTIONS_2003_ACCURATE_EOP"
    TIRF_CONVENTIONS_2003_SIMPLE_EOP = "TIRF_CONVENTIONS_2003_SIMPLE_EOP"
    TIRF_CONVENTIONS_2010_ACCURATE_EOP = "TIRF_CONVENTIONS_2010_ACCURATE_EOP"
    TIRF_CONVENTIONS_2010_SIMPLE_EOP = "TIRF_CONVENTIONS_2010_SIMPLE_EOP"
    TOD_CONVENTIONS_1996_ACCURATE_EOP = "TOD_CONVENTIONS_1996_ACCURATE_EOP"
    TOD_CONVENTIONS_1996_SIMPLE_EOP = "TOD_CONVENTIONS_1996_SIMPLE_EOP"
    TOD_CONVENTIONS_2003_ACCURATE_EOP = "TOD_CONVENTIONS_2003_ACCURATE_EOP"
    TOD_CONVENTIONS_2003_SIMPLE_EOP = "TOD_CONVENTIONS_2003_SIMPLE_EOP"
    TOD_CONVENTIONS_2010_ACCURATE_EOP = "TOD_CONVENTIONS_2010_ACCURATE_EOP"
    TOD_CONVENTIONS_2010_SIMPLE_EOP = "TOD_CONVENTIONS_2010_SIMPLE_EOP"
    TOD_WITHOUT_EOP_CORRECTIONS = "TOD_WITHOUT_EOP_CORRECTIONS"
    VEIS_1950 = "VEIS_1950"

class ManeuverTrigger():
    """Maneuver trigger events.
    """

    UNDEFINED = 0
    DATE_TIME = 1
    LONGITUDE_CROSSING = 2
    APSIDE_CROSSING = 3

class ManeuverType():
    """Maneuver types.
    """

    UNDEFINED = 0
    CONSTANT_THRUST = 1
    NORTH_SOUTH_STATIONING = 2
    EAST_WEST_STATIONING = 3
    SEMI_MAJOR_AXIS_CHANGE = 4
    PERIGEE_CHANGE = 5
    ECCENTRICITY_CHANGE = 6
    INCLINATION_CHANGE = 7
    RAAN_CHANGE = 8
    ARG_PERIGEE_CHANGE = 9
    STOP_PROPAGATION = 10

class MeasurementType():
    """Measurement types.
    """

    AZIMUTH = 0
    ELEVATION = 1
    RANGE = 2
    RANGE_RATE = 3
    RIGHT_ASCENSION = 4
    DECLINATION = 5
    POSITION = 6
    POSITION_VELOCITY = 7

class PositionAngle():
    """Orekit position angles.
    """

    MEAN = 0
    ECCENTRIC = 1
    TRUE = 2

class TDMFileFormat():
    """CCSDS TDM file formats.
    """

    KEYVALUE = 0
    XML = 1
    UNKNOWN = 2

class OutputFlag():
    """Bit flags that control propagation and estimation output.
    """

    OUTPUT_ESTM_COV = 1
    OUTPUT_PROP_COV = 2
    OUTPUT_INNO_COV = 4
    OUTPUT_RESIDUALS = 8
    OUTPUT_DENSITY = 16

class Constant():
    """Miscellaneous constants.
    """

    DEGREE_TO_RAD = pi/180
    ARC_SECOND_TO_RAD = pi/648000
    PLUS_I = (1, 0, 0)
    PLUS_J = (0, 1, 0)
    PLUS_K = (0, 0, 1)
    MINUS_I = (-1, 0, 0)
    MINUS_J = (0, -1, 0)
    MINUS_K = (0, 0, -1)
    ZERO_VECTOR = (0, 0, 0)

def configure(**kwargs)->Settings:
    """Create Settings object with default values.

    Parameters
    ----------
    **kwargs : Keyword arguments for non-default parameters.

    Returns
    -------
    Settings object initialized with defaults and given values.
    """

    cfg = Settings(rso_mass=5.0, rso_area=0.1, gravity_degree=20, gravity_order=20, ocean_tides_degree=20, ocean_tides_order=20,
                   third_body_sun=True, third_body_moon=True, solid_tides_sun=True, solid_tides_moon=True, drag_model=DragModel.MSISE2000,
                   drag_coefficient=Parameter(value=2.0, min=1.0, max=3.0, estimation=EstimationType.ESTIMATE),
                   drag_exp_rho0=3.614E-13, drag_exp_H0=700000.0, drag_exp_Hscale=88667.0, rp_sun=True,
                   rp_coeff_reflection=Parameter(value=1.5, min=1.0, max=2.0, estimation=EstimationType.ESTIMATE),
                   prop_inertial_frame=Frame.GCRF, integ_min_time_step=1E-3, integ_max_time_step=300.0, integ_abs_tolerance=1E-14,
                   integ_rel_tolerance=1E-12, output_flags=OutputFlag.OUTPUT_ESTM_COV|OutputFlag.OUTPUT_PROP_COV|OutputFlag.OUTPUT_INNO_COV|
                   OutputFlag.OUTPUT_RESIDUALS, estm_filter=Filter.UNSCENTED_KALMAN, estm_DMC_corr_time=40.0, estm_DMC_sigma_pert=5E-9,
                   estm_DMC_acceleration=Parameter(value=0.0, min=-1E-3, max=1E-3, estimation=EstimationType.ESTIMATE),
                   estm_smoother_iterations=10, estm_enable_PDAF=False, estm_detection_probability=0.99,
                   estm_gating_probability=0.99, estm_gating_threshold=5.0)
    for key, value in kwargs.items():
        attr = getattr(cfg, key)
        if (hasattr(attr, "CopyFrom")):
            attr.CopyFrom(value)
        elif (hasattr(attr, "MergeFrom")):
            attr[:] = value
        else:
            setattr(cfg, key, value)
    return(cfg)

def add_facet(cfg: Settings, normal: Tuple[float, float, float], area: float)->Facet:
    """Add a facet to a box-wing spacecraft model.

    Parameters
    ----------
    cfg : Settings object.
    normal : Facet unit outward normal vector in the spacecraft's co-moving LVLH frame.
    area : Facet area [m^2].

    Returns
    -------
    Added Facet object
    """

    cfg.rso_facets.append(Facet(normal=normal, area=area))
    return(cfg.rso_facets[-1])

def add_maneuver(cfg: Settings, time: float, trigger_event: int, trigger_params: Optional[List[float]],
                 maneuver_type: int, maneuver_params: Optional[List[float]])->Maneuver:
    """Add maneuver to propagation force models.

    Parameters
    ----------
    cfg : Settings object.
    time : TT offset from J2000 epoch [s].
    trigger_event : Constant from ManeuverTrigger.
    trigger_params : Additional data for trigger event.
    maneuver_type : Constant from ManeuverType.
    maneuver_params : Additional maneuver specific data.

    Returns
    -------
    Added Maneuver object.
    """

    cfg.maneuvers.append(Maneuver(time=time, trigger_event=trigger_event, trigger_params=trigger_params,
                                  maneuver_type=maneuver_type, maneuver_params=maneuver_params))
    return(cfg.maneuvers[-1])

def add_station(cfg: Settings, name: str, latitude: float, longitude: float, altitude: float,
                fov_azimuth: float=None, fov_elevation: float=None, fov_aperture: float=None,
                bias: Optional[List[float]]=None, bias_estimation: int=None)->Station:
    """Add a ground station.

    Parameters
    ----------
    cfg : Settings object.
    name : Station name.
    latitude : Station WGS-84 latitude [rad].
    longitude : Station WGS-84 longitude [rad].
    altitude : Station height above WGS-84 reference geoid [m].
    fov_azimuth : FOV center azimuth [rad].
    fov_elevation : FOV center elevation [rad].
    fov_aperture : FOV aperture angle [rad].
    bias : Measurement biases.
    bias_estimation : Constant from EstimationType.

    Returns
    -------
    Added Station object.
    """

    cfg.stations[name].latitude = latitude
    cfg.stations[name].longitude = longitude
    cfg.stations[name].altitude = altitude
    if (fov_azimuth is not None and fov_elevation is not None and fov_aperture is not None):
        cfg.stations[name].fov_azimuth = fov_azimuth
        cfg.stations[name].fov_elevation = fov_elevation
        cfg.stations[name].fov_aperture = fov_aperture
    if (bias is not None):
        cfg.stations[name].bias[:] = bias
    if (bias_estimation is not None):
        cfg.stations[name].bias_estimation = bias_estimation
    return(cfg.stations[name])

def build_measurement(time: float, station_name: str, values: List[float])->Measurement:
    """Build Measurement object from observed or simulated data.

    Parameters
    ----------
    time : Measurement time; TT offset from J2000 epoch [s].
    station_name : Station name.
    values : Measurement data corresponding to types configured in Settings object.

    Returns
    -------
    Measurement object.
    """

    return(Measurement(time=time, station=station_name, values=values))

def ltr_to_matrix(lower_triangle: List[float])->List[List[float]]:
    """Construct a symmetric matrix from its lower triangle.

    Parameters
    ----------
    lower_triangle : Lower triangle of a symmetric matrix in row-major order.

    Returns
    -------
    Symmetric matrix as a list of lists or None on error.
    """

    m = (sqrt(8*len(lower_triangle) + 1) - 1)/2
    n = int(m)
    if (m == n):
        k, matrix = 0, [[0]*n for _ in range(n)]
        for i in range(n):
            for j in range(i+1):
                matrix[i][j], matrix[j][i] = lower_triangle[k], lower_triangle[k]
                k += 1
        return(matrix)
    return(None)

if (__name__ != '__main__'):
    __pdoc__ = {m: False for m in ("Facet", "Maneuver", "Measurement", "Parameter", "Settings",
                                   "Station", "TDMFileFormat", "rpc", "version")}
    _root_dir = path.dirname(path.abspath(path.realpath(__file__)))
    _data_dir = path.join(_root_dir, "orekit-data")
    _libs_dir = path.join(_root_dir, "target")
    _jar_file = path.join(_libs_dir, "orbdetpy-server-{}.jar".format(__version__))
    RemoteServer.connect(_data_dir, _jar_file)
