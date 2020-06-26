# predict_passes.py - Predict satellite passes for ground stations.
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

import json
from os import path
import datetime as dt
import tkinter as tk
from orbdetpy import configure, add_station, MeasurementType, Constant
from orbdetpy.conversion import get_J2000_epoch_offset, get_UTC_string
from orbdetpy.propagation import propagate_orbits

class Application(tk.Frame):
    def __init__(self, master=None):
        super().__init__(master)
        self.cfg_file = path.join(path.dirname(path.realpath(__file__)), "pass_config.json")
        self.master = master
        self.master.minsize(800, 600)
        self.master.title("Satellite Pass Prediction")
        self.master.protocol("WM_DELETE_WINDOW", self.save_settings)
        self.pack()
        self.create_widgets()
        self.load_settings()

    def create_widgets(self):
        tk.Label(self, text="Latitude [deg]").grid(row=0, sticky="E")
        self.latitude = tk.Entry(self, width=19)
        self.latitude.grid(row=0, column=1, sticky="W")
        tk.Label(self, text="Longitude [deg]").grid(row=0, column=2, sticky="E")
        self.longitude = tk.Entry(self, width=19)
        self.longitude.grid(row=0, column=3, sticky="W")
        tk.Label(self, text="Altitude [m]").grid(row=0, column=4, sticky="E")
        self.altitude = tk.Entry(self, width=12)
        self.altitude.grid(row=0, column=5, sticky="W")

        tk.Label(self, text="FOV azimuth [deg]").grid(row=1, sticky="E")
        self.fov_azimuth = tk.Entry(self, width=19)
        self.fov_azimuth.grid(row=1, column=1, sticky="W")
        tk.Label(self, text="FOV elevation [deg]").grid(row=1, column=2, sticky="E")
        self.fov_elevation = tk.Entry(self, width=19)
        self.fov_elevation.grid(row=1, column=3, sticky="W")
        tk.Label(self, text="FOV aperture [deg]").grid(row=1, column=4, sticky="E")
        self.fov_aperture = tk.Entry(self, width=12)
        self.fov_aperture.grid(row=1, column=5, sticky="W")

        tk.Label(self, text="Start UTC").grid(row=2, sticky="E")
        self.start_time = tk.Entry(self, width=19)
        self.start_time.grid(row=2, column=1, sticky="W")
        tk.Label(self, text="End UTC").grid(row=2, column=2, sticky="E")
        self.end_time = tk.Entry(self, width=19)
        self.end_time.grid(row=2, column=3, sticky="W")
        tk.Label(self, text="Step size [s]").grid(row=2, column=4, sticky="E")
        self.step_size = tk.Entry(self, width=12)
        self.step_size.grid(row=2, column=5, sticky="W")

        tk.Label(self, text="Two line elements").grid(row=3, sticky="E")
        self.tle = tk.Text(self, height=9, width=72)
        self.tle.grid(row=3, column=1, columnspan=6, sticky="W")

        self.predict = tk.Button(self, text="Predict passes", command=self.predict)
        self.predict.grid(row=4, columnspan=6)

        tk.Label(self, text="UTC, Azimuth [deg], Elevation [deg]").grid(row=5, columnspan=6)
        self.output = tk.Text(self, height=18, width=90)
        self.output.grid(row=6, columnspan=6)

    def load_settings(self):
        cfg = {}
        if (path.isfile(self.cfg_file)):
            with open(self.cfg_file, "r") as fp:
                cfg = json.load(fp)

        d0 = dt.datetime.today()
        d1 = d0 + dt.timedelta(days=1)
        self.latitude.insert(0, cfg.get("latitude", 30.6714))
        self.longitude.insert(0, cfg.get("longitude", -104.0219))
        self.altitude.insert(0, cfg.get("altitude", 2070.0))
        self.fov_azimuth.insert(0, cfg.get("fov_azimuth", 0.0))
        self.fov_elevation.insert(0, cfg.get("fov_elevation", 0.0))
        self.fov_aperture.insert(0, cfg.get("fov_aperture", 0.0))
        self.start_time.insert(0, dt.datetime(d0.year, d0.month, d0.day).strftime("%Y-%m-%dT%H:%M:%SZ"))
        self.end_time.insert(0, dt.datetime(d1.year, d1.month, d1.day).strftime("%Y-%m-%dT%H:%M:%SZ"))
        self.step_size.insert(0, 60.0)
        self.tle.insert(tk.END, cfg.get("tle", ""))
        self.output.insert(tk.END, cfg.get("output", ""))

    def save_settings(self):
        cfg = {}
        for f in ["latitude", "longitude", "altitude", "fov_azimuth", "fov_elevation",
                  "fov_aperture", "tle", "output"]:
            cfg[f] = getattr(self, f).get("0.0", "end-1c") if (
                f in ["tle", "output"]) else getattr(self, f).get()

        self.master.destroy()
        with open(self.cfg_file, "w") as fp:
            json.dump(cfg, fp, indent=1)

    def predict(self):
        self.predict["state"] = "disabled"
        start = get_J2000_epoch_offset(self.start_time.get())
        end = get_J2000_epoch_offset(self.end_time.get())
        data, tle = {}, [l for l in self.tle.get("0.0", "end-1c").splitlines()
                         if l.startswith("1") or l.startswith("2")]
        for f in ["latitude", "longitude", "altitude", "fov_azimuth", "fov_elevation",
                  "fov_aperture", "step_size"]:
            data[f] = float(getattr(self, f).get())
            if (f not in ["altitude", "step_size"]):
                data[f] *= Constant.DEGREE

        cfg_list = []
        for i in range(0, len(tle), 2):
            cfg_list.append(configure(prop_start=start, prop_initial_TLE=tle[i:i+2],
                                      prop_end=end, prop_step=data["step_size"],
                                      sim_measurements=True))
            add_station(cfg_list[-1], "Sensor", data["latitude"], data["longitude"],
                        data["altitude"], data["fov_azimuth"], data["fov_elevation"],
                        data["fov_aperture"])
            cfg_list[-1].measurements[MeasurementType.AZIMUTH].error[:] = [0.0]
            cfg_list[-1].measurements[MeasurementType.ELEVATION].error[:] = [0.0]

        if (len(cfg_list)):
            i = 0
            self.output.delete("0.0", tk.END)
            for o in propagate_orbits(cfg_list):
                self.output.insert(tk.END, "\nObject {}:\n".format(tle[i][2:7]))
                i += 2
                for m in o.array:
                    self.output.insert(tk.END, "{}: {:.5f}, {:.5f}\n".format(
                        get_UTC_string(m.time), (m.values[0]/Constant.DEGREE + 360)%360,
                        m.values[1]/Constant.DEGREE))

        self.predict["state"] = "normal"

Application(master=tk.Tk()).mainloop()
