# setup.py - PyPI installation builder.
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

from setuptools import setup
from distutils.util import convert_path
from os import path, environ

# Get the version number from the version path.
version_dict = {}
version_path = convert_path("orbdetpy/version.py")
with open(version_path, "r") as ver_file:
    exec(ver_file.read(), version_dict)

# Get the long description from the README file.
here = path.dirname(path.abspath(__file__))
with open(path.join(here, "README.md"), "r") as a_file:
    long_description = a_file.read()

with open(path.join(here, "requirements.txt"), "r") as a_file:
    requirements = a_file.read()

CLASSIFIERS = ["Development Status :: 5 - Production/Stable",
               "Environment :: Console",
               "Intended Audience :: Developers",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
               "Natural Language :: English",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Programming Language :: Python :: 3",
               "Programming Language :: Python :: 3.7",
               "Programming Language :: Python :: 3.8",
               "Programming Language :: Java",
               "Topic :: Scientific/Engineering",
               "Topic :: Scientific/Engineering :: Astronomy",
               "Topic :: Utilities"]
PROJECT_URLS = {"Documentation": "https://github.com/ut-astria/orbdetpy",
                "Source": "https://github.com/ut-astria/orbdetpy",
                "Tracker": "https://github.com/ut-astria/orbdetpy/issues"}
PACKAGE_DATA = {"orbdetpy": ["*.py", "orekit-data/*", "orekit-data/**/*", "orekit-data/**/**/*",
                             "plotting/*.py", "rpc/*.py", "target/orbdetpy-server*.jar"]}

setup(name="orbdetpy",
      packages=["orbdetpy"],
      version=version_dict["__version__"],
      description="Orbit determination routines for Python",
      long_description=long_description,
      long_description_content_type="text/markdown",
      author="Shiva Iyer",
      author_email="shiva.iyer@utexas.edu",
      url="https://github.com/ut-astria/orbdetpy",
      project_urls=PROJECT_URLS,
      keywords=["orbit_determination utilities orbit space celestial_mechanics "
                "astrodynamics estimation satellite_tracking pass_prediction"],
      classifiers=CLASSIFIERS,
      package_data=PACKAGE_DATA,
      install_requires=requirements)
