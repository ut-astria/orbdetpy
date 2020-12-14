# simulation.py - Plot simulation results.
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

import matplotlib.pyplot as plt
from math import acos, pi
from numpy import array, cross, dot
from numpy.linalg import norm

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

    mu = 398600.4418
    tstamp,hvec,hmag,ener,sma,ecc,inc,raan,argp,tran = [],[],[],[],[],[],[],[],[],[]
    for o in sim_data:
        if (hasattr(o, "true_state")):
            rv = [x/1000.0 for x in o.true_state[:6]]
        else:
            rv = [x/1000.0 for x in o.estimated_state[:6]]
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

    outfiles = []
    hvec = array(hvec)
    tim = [(t - tstamp[0])/3600 for t in tstamp]

    plt.figure(0)
    plt.subplot(611)
    plt.scatter(tim, sma, marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel("a [km]")
    plt.subplot(612)
    plt.scatter(tim, ecc, marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel("e")
    plt.subplot(613)
    plt.scatter(tim, inc, marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel("i [deg]")
    plt.subplot(614)
    plt.scatter(tim, raan, marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel(r"$\Omega$ [deg]")
    plt.subplot(615)
    plt.scatter(tim, argp, marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel(r"$\omega$ arg. [deg]")
    plt.subplot(616)
    plt.scatter(tim, tran, marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel(r"$\theta$ [deg]")
    plt.suptitle("Orbital Elements")
    if (output_file_path is not None):
        outfiles.append(output_file_path + "_elements.png")
        plt.savefig(outfiles[-1], format = "png")

    plt.figure(1)
    plt.subplot(211)
    plt.scatter(tim, hmag, marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel("Angular Momentum")
    plt.subplot(212)
    plt.scatter(tim, ener, marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel("Specific Energy")
    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (output_file_path is not None):
        outfiles.append(output_file_path + "_momentum_energy.png")
        plt.savefig(outfiles[-1], format = "png")

    plt.figure(2)
    plt.subplot(311)
    plt.scatter(tim, hvec[:,0], marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel("h(x)")
    plt.subplot(312)
    plt.scatter(tim, hvec[:,1], marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel("h(y)")
    plt.subplot(313)
    plt.scatter(tim, hvec[:,2], marker = "o", s = 7)
    plt.xlabel("Time [hr]")
    plt.ylabel("h(z)")
    plt.suptitle("Angular Momentum Components")
    plt.tight_layout(rect = [0, 0.03, 1, 0.95])
    if (output_file_path is not None):
        outfiles.append(output_file_path + "_momentum_components.png")
        plt.savefig(outfiles[-1], format = "png")

    if (interactive):
        plt.show()
    plt.close("all")
    return(outfiles)
