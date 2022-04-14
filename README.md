**orbdetpy** is a Python orbit determination library.

# Features

Dynamics in orbdetpy can be configured with:

1. EGM96 gravity field up to degree and order 360.
2. Earth solid tides due to the influence of the Sun and Moon.
3. FES 2004 ocean tide model up to degree and order 100.
4. NRL MSISE-00 and exponential atmospheric drag models.
5. Solar radiation pressure.
6. Third body perturbations from the Sun and Moon.
7. Satellite box-wing models and maneuvers.

Range, range-rate, angles, and inertial state measurements are supported. Filtering can be done with an EKF, UKF, or Batch Least Squares. Dynamic Model Compensation (DMC) can be used to estimate unmodeled accelerations.

# Installation

If you have docker installed and wish to use a pre-built environment refer to the [Docker](#docker) section in this README.

1. Install Java SE 11 (11.0.10) from <https://www.oracle.com/javadownload>. Set the `JAVA_HOME` environment variable to the Java installation folder. The `java` executable must be added to the system path.

2. Install Python 3.8.0 or higher and run `pip install orbdetpy` to install orbdetpy and other package dependencies.

3. Update the astrodynamics data under `orbdetpy/orekit-data` periodically by running the following. You will need `root` privileges on some systems.

   `python -c "from orbdetpy.astro_data import update_data; update_data();"`

# Development

1. Download and extract <https://github.com/ut-astria/orbdetpy/releases/download/2.1.0/orekit-data.tar.gz> under the `orbdetpy/` sub-folder.

2. Developers will need Apache Maven 3+ to build the Java library. Build using the following from the `orbdetpy/` sub-folder, where `os_cpu_type` is `linux-x86_64`, `linux-x86_32`, `windows-x86_64`, `windows-x86_32`, `osx-x86_64`, or `osx-x86_32` depending on your CPU and OS:

   `mvn -e -Dos.detected.classifier=os_cpu_type package`

   The command-line is simpler on Intel/AMD 64-bit Linux:

   `mvn -e package`

3. Run `pip install -e ./` from the top-level folder containing `setup.py`.

# Docker

1. Build the docker image on a machine that has docker installed. Go to the root folder of this repository where the `Dockerfile` is and run, `docker build --build-arg ORBDETPY_VERSION=2.1.0 -t orbdetpy:2.1.0 .`

2. Run the *orbdetpy:2.1.0* image in a daemon state: `docker run -it --rm orbdetpy:2.1.0 bash`

3. Activate the python environment and run a test to determine a successful docker image build:

```bash
cd && . env_orbdetpy/bin/activate && python orbdetpy/examples/test_estimation.py
```

4. From here, you can either develop in orbdetpy or script and test in this pre-built environment.

# Examples

1. `fit_radec.py` : Run OD with real angles measurements. Also demonstrates the Laplace IOD method for estimating an initial state vector.

2. `interpolate_oem.py` : Command-line tool for interpolating state vectors from CCSDS OEM ephemeris files. OEM files are available for download at <http://astria.tacc.utexas.edu/AstriaGraph>.

3. `predict_passes.py` : Predict satellite passes for ground stations or geographic regions using TLEs. Current elements may be obtained from sites such as <http://www.celestrak.com>.

4. `propagate_tle.py` : Propagate TLEs given by command-line arguments.

5. `test_conversion.py` : Test reference frame and other conversion functions.

6. `test_estimation.py` : Demonstrates measurement simulation and orbit determination functions.

7. `test_interpolation.py` : Interpolate state vectors.

# Known Issues

1. Give orbdetpy network access permissions if you get warnings from the Microsoft Windows Firewall. 

2. If you use orbdetpy with Python's `multiprocessing` package, call the function `multiprocessing.set_start_method("spawn")` before other `multiprocessing` calls.
