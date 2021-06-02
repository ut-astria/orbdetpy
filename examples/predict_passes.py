# predict_passes.py - Predict satellite passes for ground stations.
# Copyright (C) 2020-2021 University of Texas
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
from orbdetpy import add_station, configure, Constant, DragModel, Frame, MeasurementType
from orbdetpy.conversion import get_J2000_epoch_offset, get_UTC_string, pos_to_lla
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
        r = 0
        label = """Enter one Lat/Lon/Alt and optional sensor FOV for a single ground station -or-
        Enter 3+ comma separated Lat/Lon vertices in CCW order for geographic regions. 
        """
        tk.Label(self, text=label).grid(row=r, columnspan=6)

        r += 1
        tk.Label(self, text="Latitude [deg]").grid(row=r, sticky="E")
        self.latitude = tk.Entry(self, width=19)
        self.latitude.grid(row=r, column=1, sticky="W")
        tk.Label(self, text="Longitude [deg]").grid(row=r, column=2, sticky="E")
        self.longitude = tk.Entry(self, width=19)
        self.longitude.grid(row=r, column=3, sticky="W")
        tk.Label(self, text="Altitude [m]").grid(row=r, column=4, sticky="E")
        self.altitude = tk.Entry(self, width=12)
        self.altitude.grid(row=r, column=5, sticky="W")

        r += 1
        tk.Label(self, text="FOV azimuth [deg]").grid(row=r, sticky="E")
        self.fov_azimuth = tk.Entry(self, width=19)
        self.fov_azimuth.grid(row=r, column=1, sticky="W")
        tk.Label(self, text="FOV elevation [deg]").grid(row=r, column=2, sticky="E")
        self.fov_elevation = tk.Entry(self, width=19)
        self.fov_elevation.grid(row=r, column=3, sticky="W")
        tk.Label(self, text="FOV aperture [deg]").grid(row=r, column=4, sticky="E")
        self.fov_aperture = tk.Entry(self, width=12)
        self.fov_aperture.grid(row=r, column=5, sticky="W")

        r += 1
        tk.Label(self, text="Start UTC").grid(row=r, sticky="E")
        self.start_time = tk.Entry(self, width=19)
        self.start_time.grid(row=r, column=1, sticky="W")
        tk.Label(self, text="End UTC").grid(row=r, column=2, sticky="E")
        self.end_time = tk.Entry(self, width=19)
        self.end_time.grid(row=r, column=3, sticky="W")
        tk.Label(self, text="Step size [s]").grid(row=r, column=4, sticky="E")
        self.step_size = tk.Entry(self, width=12)
        self.step_size.grid(row=r, column=5, sticky="W")

        r += 1
        tk.Label(self, text="Two line elements").grid(row=r, sticky="E")
        self.tle = tk.Text(self, height=9, width=72)
        self.tle.grid(row=r, column=1, columnspan=6, sticky="W")

        r += 1
        self.predict = tk.Button(self, text="Predict passes", command=self.predict)
        self.predict.grid(row=r, columnspan=6)

        r += 1
        self.output_label = tk.Label(self, text="Output")
        self.output_label.grid(row=r, columnspan=6)

        r += 1
        self.output = tk.Text(self, height=15, width=90)
        self.output.grid(row=r, columnspan=6)

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
        self.output.insert(tk.END, "")

    def save_settings(self):
        cfg = {}
        for f in ["latitude", "longitude", "altitude", "fov_azimuth", "fov_elevation",
                  "fov_aperture", "tle"]:
            cfg[f] = getattr(self, f).get("0.0", "end-1c") if (f == "tle") else getattr(self, f).get()

        self.master.destroy()
        with open(self.cfg_file, "w") as fp:
            json.dump(cfg, fp, indent=1)

    def predict(self):
        self.predict["state"] = "disabled"
        start = get_J2000_epoch_offset(self.start_time.get())
        end = get_J2000_epoch_offset(self.end_time.get())
        data, tle = {}, [l for l in self.tle.get("0.0", "end-1c").splitlines()
                         if l.startswith("1") or l.startswith("2")]

        for f in ["latitude", "longitude", "altitude", "fov_azimuth", "fov_elevation", "fov_aperture", "step_size"]:
            data[f] = [float(t.strip()) for t in getattr(self, f).get().split(",")]
            if (f not in ["altitude", "step_size"]):
                data[f] = [d*Constant.DEGREE_TO_RAD for d in data[f]]

        cfg_list = []
        sim_meas = len(data["latitude"]) <= 2 or len(data["latitude"]) != len(data["longitude"])
        for i in range(0, len(tle), 2):
            cfg_list.append(configure(prop_start=start, prop_initial_TLE=tle[i:i+2], prop_end=end, prop_step=data["step_size"][0],
                                      sim_measurements=sim_meas, gravity_degree=-1, gravity_order=-1, ocean_tides_degree=-1,
                                      ocean_tides_order=-1, third_body_sun=False, third_body_moon=False, solid_tides_sun=False,
                                      solid_tides_moon=False, drag_model=DragModel.UNDEFINED, rp_sun=False))
            if (not sim_meas):
                cfg_list[-1].geo_zone_lat_lon[:] = [l for ll in zip(data["latitude"], data["longitude"]) for l in ll]
                continue

            add_station(cfg_list[-1], "Sensor", data["latitude"][0], data["longitude"][0], data["altitude"][0],
                        data["fov_azimuth"][0], data["fov_elevation"][0], data["fov_aperture"][0])
            cfg_list[-1].measurements[MeasurementType.AZIMUTH].error[:] = [0.0]
            cfg_list[-1].measurements[MeasurementType.ELEVATION].error[:] = [0.0]

        if (len(cfg_list)):
            i = 0
            self.output.delete("0.0", tk.END)
            if (sim_meas):
                self.output_label["text"] = "UTC, Azimuth [deg], Elevation [deg]"
            else:
                self.output_label["text"] = "UTC, Latitude [deg], Longitude [deg], Altitude [m]"

            for o in propagate_orbits(cfg_list):
                self.output.insert(tk.END, "\nObject {}:\n".format(tle[i][2:7]))
                i += 2
                for m in o.array:
                    if (sim_meas):
                        self.output.insert(tk.END, "{}: {:.5f}, {:.5f}\n".format(
                            get_UTC_string(m.time), (m.values[0]/Constant.DEGREE_TO_RAD + 360)%360,
                            m.values[1]/Constant.DEGREE_TO_RAD))
                    else:
                        lla = pos_to_lla(Frame.GCRF, m.time, m.true_state)
                        self.output.insert(tk.END, "{}: {:.5f}, {:.5f}, {:.2f}\n".format(
                            get_UTC_string(m.time), lla[0]/Constant.DEGREE_TO_RAD, lla[1]/Constant.DEGREE_TO_RAD, lla[2]))

        self.predict["state"] = "normal"

Application(master=tk.Tk()).mainloop()
