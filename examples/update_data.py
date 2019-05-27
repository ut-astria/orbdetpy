# update_data.py - Update Orekit astrodynamics data files.
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

from os import path
import requests

def format_weather(lines):
    c1 = [0, 5,  8, 11, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 47, 51,
          55, 59, 63, 67, 71, 75, 79, 83, 87, 89, 93,  99, 101, 107,
          113, 119, 125]
    c2 = [5, 8, 11, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 47, 51, 55,
          59, 63, 67, 71, 75, 79, 83, 87, 89, 93, 99, 101, 107, 113,
          119, 125, 131]

    fmtdata = ""
    for line in lines.splitlines():
        if (line == "END DAILY_PREDICTED"):
            break
        if (len(line) == 0 or not line[0].isnumeric()):
            continue

        fmtdata += ",".join([line[i:j] for i,j in zip(c1, c2)]) + "\n"

    return(fmtdata)

datadir = path.join(path.dirname(path.dirname(path.abspath(__file__))), "data")

updates = [["https://datacenter.iers.org/data/latestVersion/7_FINALS.ALL_IAU1980_V2013_017.txt",
            path.join(datadir, "Earth-Orientation-Parameters", "IAU-1980", "finals.all"), None],
           ["https://datacenter.iers.org/data/latestVersion/9_FINALS.ALL_IAU2000_V2013_019.txt",
            path.join(datadir, "Earth-Orientation-Parameters", "IAU-2000", "finals2000A.all"), None],
           ["http://maia.usno.navy.mil/ser7/tai-utc.dat",
            path.join(datadir, "tai-utc.dat"), None],
           ["http://www.celestrak.com/SpaceData/SW-All.txt",
            path.join(datadir, "SpaceWeather.dat"), format_weather]]

for u in updates:
    print("Updating %s" % u[1])
    resp = requests.get(u[0])
    if (resp.status_code != 200):
        print("Error %s in %s", (resp.status_code, u[0]))
        continue

    if (u[2] is not None):
        data = u[2](resp.text)
    else:
        data = resp.text

    with open(u[1], "w") as f:
        f.write(data)
