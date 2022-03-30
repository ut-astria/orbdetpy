# setup.py - PyPI installation builder.
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

from os import path
from setuptools import setup

root_dir, version = path.dirname(path.abspath(path.realpath(__file__))), {}
with open(path.join(root_dir, "README.md"), "r") as fp:
    long_description = fp.read()
with open(path.join(root_dir, "requirements.txt"), "r") as fp:
    requirements = fp.read()
with open(path.join(root_dir, "orbdetpy", "version.py"), "r") as fp:
    exec(fp.read(), version)

CLASSIFIERS = ["Development Status :: 5 - Production/Stable",
               "Environment :: Console",
               "Intended Audience :: Developers",
               "Intended Audience :: Science/Research",
               "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
               "Natural Language :: English",
               "Operating System :: OS Independent",
               "Programming Language :: Python",
               "Programming Language :: Python :: 3",
               "Programming Language :: Python :: 3.8",
               "Programming Language :: Python :: 3.9",
               "Programming Language :: Java",
               "Topic :: Scientific/Engineering",
               "Topic :: Scientific/Engineering :: Astronomy",
               "Topic :: Utilities"]
PROJECT_URLS = {"Documentation": "https://ut-astria.github.io/orbdetpy", "Source": "https://github.com/ut-astria/orbdetpy",
                "Tracker": "https://github.com/ut-astria/orbdetpy/issues"}
PACKAGE_DATA = {"": ["requirements.txt", "setup.cfg", "setup.py"],
                "docs": ["*.html", "plotting/*.html", "rpc/*.html"],
                "examples": ["*.json", "*.py"],
                "orbdetpy": ["*.py", "*.sh", "plotting/*.py", "rpc/*.py", f"""target/orbdetpy-server-{version["__version__"]}.jar""",
                             "orekit-data/*", "orekit-data/**/*", "orekit-data/**/**/*"],
                "tests": ["*.py", "*.txt"]}

setup(name="orbdetpy", packages=list(PACKAGE_DATA.keys()), version=version["__version__"], description="Orbit determination routines for Python",
      long_description=long_description, long_description_content_type="text/markdown", author="Shiva Iyer",
      author_email="shiva.iyer@utexas.edu", url="https://github.com/ut-astria/orbdetpy", project_urls=PROJECT_URLS,
      keywords=["orbit_determination utilities orbit space celestial_mechanics astrodynamics estimation satellite_tracking pass_prediction"],
      classifiers=CLASSIFIERS, package_data=PACKAGE_DATA, install_requires=requirements, license="GPLv3+")
