# estimation.py - Orbit estimation functions.
# Copyright (C) 2019-2022 University of Texas
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

from typing import List, Tuple
from traceback import format_exc
from orbdetpy.rpc.estimation_pb2_grpc import EstimationStub
from orbdetpy.rpc.messages_pb2 import AnglesInput, DetermineOrbitInput, MeasurementArray, MultiTargetInput, Settings
from orbdetpy.rpc.server import RemoteServer

def determine_orbit(config: List[Settings], meas):
    """Run orbit determination for the given objects and measurements.

    Parameters
    ----------
    config : List of Settings objects.
    meas : List of list of Measurement objects.

    Returns
    -------
    Orbit determination results.
    """

    od_output, requests = [], []
    for c, m in zip(config, meas):
        inp = DetermineOrbitInput(config=c)
        inp.measurements.extend(m)
        requests.append(_estimation_stub.determineOrbit.future(inp))

    for req in requests:
        try:
            fit_data = req.result().array
        except Exception as exc:
            fit_data = format_exc()
        od_output.append(fit_data)
    return(od_output)

def multi_target_OD(config_list: List[Settings], meas_list):
    """ Run multiple target orbit determination using CAR/MHF and JPDA.

    Parameters
    ----------
    config_list : List of Settings objects.
    meas_list : List of list of Measurement objects.

    Returns
    -------
    Multiple target orbit determination results.
    """

    # Disable DMC because CAR-MHF currently only handles static process noise
    for cfg in config_list:
        cfg.estm_DMC_corr_time = 0.0
        cfg.estm_DMC_sigma_pert = 0.0

    resp = _estimation_stub.multiTargetOD(MultiTargetInput(
        config=config_list, measurements=(MeasurementArray(array=m) for m in meas_list)))
    return(resp)

def iod_laplace(frame: int, lat: float, lon: float, alt: float, time: Tuple[float, float, float],
                ra: Tuple[float, float, float], dec: Tuple[float, float, float])->List[float]:
    """Estimate orbit from 3 RA/Dec angles using the Laplace method.

    Parameters
    ----------
    frame : Estimation reference frame; a constant from Frame.
    lat : Observer WGS-84 latitude [rad].
    lon : Observer WGS-84 longitude [rad].
    alt : Observer height above WGS-84 reference ellipsoid [m].
    time : Times [t1, t2, t3]; each a TT offset from J2000 epoch [s].
    ra : List of 3 Right Ascensions.
    dec : List of 3 Declinations.

    Returns
    -------
    Position and velocity estimate at time t2.
    """

    resp = _estimation_stub.iodLaplace(AnglesInput(time=time, angle1=ra, angle2=dec, latitude=lat,
                                                   longitude=lon, altitude=alt, frame=frame))
    return(resp.array)

if (__name__ != "__main__"):
    _estimation_stub = EstimationStub(RemoteServer.channel())
