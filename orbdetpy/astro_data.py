# astro_data.py - Update Orekit astrodynamics data files.
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

import requests
from os import path
from orbdetpy import _data_dir

def format_weather(lines: str)->str:
    """Re-format space weather data into a more efficient form.
    """

    c1 = [0, 5,  8, 11, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 47, 51, 55, 59, 63, 67,
          71, 75, 79, 83, 87, 89, 93,  99, 101, 107, 113, 119, 125]
    c2 = [5, 8, 11, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 47, 51, 55, 59, 63, 67, 71,
          75, 79, 83, 87, 89, 93, 99, 101, 107, 113, 119, 125, 131]

    data = ""
    for line in lines.splitlines():
        if (line == "END DAILY_PREDICTED"):
            break
        if (len(line) > 0 and line[0].isnumeric()):
            data += ",".join([line[i:j] for i, j in zip(c1, c2)]) + "\n"

    return(data)

def update_data()->None:
    """Download and re-format astrodynamics data from multiple sources.
    """

    updates = [["http://www.celestrak.com/SpaceData/SW-All.txt", path.join(_data_dir, "SpaceWeather.dat"), format_weather],
               ["https://datacenter.iers.org/data/latestVersion/7_FINALS.ALL_IAU1980_V2013_017.txt",
                path.join(_data_dir, "Earth-Orientation-Parameters", "IAU-1980", "finals.all"), None],
               ["https://datacenter.iers.org/data/latestVersion/9_FINALS.ALL_IAU2000_V2013_019.txt",
                path.join(_data_dir, "Earth-Orientation-Parameters", "IAU-2000", "finals2000A.all"), None],
               ["http://maia.usno.navy.mil/ser7/tai-utc.dat", path.join(_data_dir, "tai-utc.dat"), None]]

    for u in updates:
        print("Updating {}".format(u[1]))
        try:
            resp = requests.get(u[0], timeout=1.0)
            if (resp.status_code != requests.codes.ok):
                print("Error {} in {}".format(resp.status_code, u[0]))
                continue
        except Exception as exc:
            print(exc)
            continue

        with open(u[1], "w") as f:
            f.write(u[2](resp.text) if (u[2] is not None) else resp.text)
