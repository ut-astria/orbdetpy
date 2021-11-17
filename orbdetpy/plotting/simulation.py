# simulation.py - Plot simulation results.
# Copyright (C) 2018-2021 University of Texas
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

from math import acos, pi, sqrt
from numpy import array, cross, dot
from numpy.linalg import norm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

def plot(sim_data, interactive: bool=False, output_file_path: str=None):
    """Plot simulated orbital elements, angular momenta, and specific energy.

    Parameters
    ----------
    sim_data : Return value from propagate_orbits().
    interactive : Show interactive plots if True.
    output_file_path : File path and name prefixes if plots are to be saved.

    Returns
    -------
    List of plot files if they were saved to disk.
    """

    mu, earth_equ_r = 398600.4418, 6378.1363
    tstamp,hvec,hmag,ener,sma,ecc,inc,raan,argp,tran,period,altitude = [],[],[],[],[],[],[],[],[],[],[],[]
    for o in sim_data:
        if (hasattr(o, "true_state")):
            rv = [x/1000.0 for x in o.true_state[:6]]
        else:
            rv = [x/1000.0 for x in o.estimated_state[:6]]
        if (not rv):
            continue
        rn = norm(rv[:3])
        vn = norm(rv[3:])
        h = cross(rv[:3], rv[3:])
        hn = norm(h)
        a = 1.0/(2.0/rn-vn**2/mu)
        e = cross(rv[3:], h)/mu - rv[:3]/rn
        en = norm(e)
        i = acos(h[2]/hn)*180.0/pi
        n = cross([0, 0, 1], h)
        nn = norm(n)
        O = acos(dot(n, [1, 0, 0])/nn)*180.0/pi
        if (n[1] < 0.0):
            O = 360.0 - O
        w = acos(dot(n, e)/(nn*en))*180.0/pi
        if (e[2] < 0.0):
            w = 360.0 - w
        theta = acos(dot(rv[:3], e)/(rn*en))*180.0/pi
        if (dot(rv[:3], rv[3:]) < 0):
            theta = 360.0 - theta
        tstamp.append(o.time)
        hvec.append(h)
        hmag.append(hn)
        ener.append(0.5*vn**2 - mu/rn)
        sma.append(a)
        ecc.append(en)
        inc.append(i)
        raan.append(O)
        argp.append(w)
        tran.append(theta)
        period.append(2.0*pi*sqrt(a**3/mu)/60.0)
        altitude.append(rn - earth_equ_r)

    outfiles = []
    hvec = array(hvec)
    tim = [(t - tstamp[0])/3600 for t in tstamp]

    plt.figure(0)
    plt.suptitle("Orbital Elements")
    for i in range(6):
        axis = plt.subplot(6, 1, i + 1)
        axis.ticklabel_format(useOffset=False)
        plt.scatter(tim, (sma, ecc, inc, raan, argp, tran)[i], c="b", marker="o", s=3)
        plt.xlabel("Time [hr]")
        plt.ylabel(("a [km]", "e", "i [deg]", r"$\Omega$ [deg]", r"$\omega$ [deg]", r"$\theta$ [deg]")[i])
    if (output_file_path):
        outfiles.append(output_file_path + "_elements.png")
        plt.savefig(outfiles[-1], format="png")

    plt.figure(1)
    plt.suptitle("Specific Angular Momentum & Energy")
    for i in range(2):
        axis = plt.subplot(2, 1, i + 1)
        axis.ticklabel_format(useOffset=False)
        plt.scatter(tim, (hmag, ener)[i], c="b", marker="o", s=3)
        plt.xlabel("Time [hr]")
        plt.ylabel((r"h [$km^2/s$]", r"E [$km^2/s^2$]")[i])
    plt.tight_layout(rect=(0, 0.03, 1, 0.95))
    if (output_file_path):
        outfiles.append(output_file_path + "_momentum_energy.png")
        plt.savefig(outfiles[-1], format="png")

    fig = plt.figure(2)
    plt.suptitle("Specific Angular Momentum")
    axis = fig.add_subplot(111, projection="3d")
    axis.ticklabel_format(useOffset=False)
    axis.scatter(hvec[:,0], hvec[:,1], hvec[:,2], c="b", marker="o", s=3)
    axis.grid(True)
    axis.set_xlabel(r"h(x) [$km^2/s$]")
    axis.set_ylabel(r"h(y) [$km^2/s$]")
    axis.set_zlabel(r"h(z) [$km^2/s$]")
    plt.tight_layout(rect=(0, 0.03, 1, 0.95))
    if (output_file_path):
        outfiles.append(output_file_path + "_momentum_3D.png")
        plt.savefig(outfiles[-1], format="png")

    fig = plt.figure(3)
    plt.suptitle("Gabbard Plot")
    plt.ticklabel_format(useOffset=False)
    plt.scatter(period, altitude, c="b", marker="o", s=3)
    plt.grid(True)
    plt.xlabel("Time [min]")
    plt.ylabel("Altitude [km]")
    plt.tight_layout(rect=(0, 0.03, 1, 0.95))
    if (output_file_path):
        outfiles.append(output_file_path + "_gabbard.png")
        plt.savefig(outfiles[-1], format="png")

    if (interactive):
        plt.show()
    plt.close("all")
    return(outfiles)
