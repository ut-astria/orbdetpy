Introduction
------------

This is orbdetpy, a library of Python and Java routines for orbit
determination. It is a thin Python wrapper for our estimation tools
and Orekit, which are both written in Java. 

Features
--------

The force model for orbit propagation currently includes:

1. EGM96 gravity field up to degree and order 360.
2. Earth solid tides due to the influence of the Sun and Moon.
3. FES 2004 ocean tide model up to degree and order 100.
4. The NRL MSISE-00 and simple exponential models for atmospheric drag.
5. Solar radiation pressure.
6. Third body perturbations from the Sun and Moon.

The measurement model supports range, range-rate, angles, and inertial
Cartesian coordinates. Filtering is done with our Unscented Kalman Filter
or Orekit's Extended Kalman Filter. Dynamic Model Compensation
(DMC) can be used with either filter to estimate additional perturbing
acclerations that result from unmodeled dynamics, maneuvers etc.

Installation
------------

1. Install the Java Development Kit 8 (1.8) from
   <http://openjdk.java.net/install/index.html>. Set the JAVA_HOME
   environment variable to point to your JDK installation. The `java`
   executable must also be in your system path.

2. Install Python 3.6+ and run `pip install orbdetpy` to install orbdetpy 
   and other package dependencies from the Python Package Index (PyPI).
   If you wish to use the `develop` or other experimental branches from
   GitHub, `git clone` them and run `pip install -e .` from the top level
   `orbdetpy` folder.

3. Source code, example programs and data files can be downloaded from 
   <https://github.com/ut-astria/orbdetpy>.

4. Update the astrodynamics data in `orbdetpy/data` periodically by running
   the following. You might need to do so as root on Unix-like systems.

   `python -c "from orbdetpy.astrodata import update_data; update_data();"`

5. Apache Maven 3+ is needed if you hack the Java code and need to
   rebuild the JAR files. Switch to the `orbdetpy/` folder and run the
   following, where `os_cpu_type` is `linux-x86_64`, `linux-x86_32`,
   `windows-x86_64`, `windows-x86_32`, `osx-x86_64`, or `osx-x86_32`,
   depending on your CPU architecture and OS.

   `mvn -e -Dos.detected.classifier=os_cpu_type package`

   If you are on Intel/AMD 64-bit Linux the command-line simplifies to:

   `mvn -e package`

Examples
--------

orbdetpy uses JSON files to store settings, measurements and estimation
results. The files in `examples/data` show how to configure measurement
simulation and orbit determination using radar or telescope data. The
file `docs/file_formats.md` documents the structure of the JSON files.

1. `run_tests.py` : Run all the use cases under `examples/data`. Simulated
   measurements, orbit fits, differences between simulated truth versus
   estimates, and 3-sigma of estimated covariances will be written to
   `output/` sub-folders.

   orbdetpy supports asynchronous calls to CPU intensive functions
   enabling parallelization in multi-processor environments. This
   program demonstrates asynchronous operations.

2. Simulate state vectors and radar measurements:

   `from orbdetpy.simulation import simulate_measurements`

   `simulate_measurements("radar_sim_cfg.json", output_file = "sim_data.json")`

   This will run the simulation configured in `radar_sim_cfg.json` and
   write simulated output to `sim_data.json`.

3. Plot simulation results:

   `from orbdetpy.plotting.simulation import plot as sim_plot`

   `sim_plot("radar_sim_cfg.json", "sim_data.json", interactive = True)`

   This will plot the simulated data generated in (2).

4. Run OD on simulated radar data:

   `from orbdetpy.estimation import determine_orbit`

   `determine_orbit("radar_od_cfg.json", "sim_data.json", output_file = "od_output.json")`

   This will run OD on the simulated radar data generated in (2)
   using the OD configuration in `radar_od_cfg.json` and write OD
   output to `od_output.json`.

5. Plot OD results:

   `from orbdetpy.plotting.estimation import plot as od_plot`

   `od_plot("radar_od_cfg.json", "sim_data.json", "od_output.json", interactive = True)`

   This will plot the OD results from (4).

6. `propagate_TLE.py` : Propagate TLEs given by command-line arguments.

Known Issues
------------

1. In Microsoft Windows, you might receive warnings from the "Windows
   Defender Firewall" when you import `orbdetpy`. You can safely allow
   `orbdetpy` access to the network because this only involves a single
   TCP port on `localhost`.
