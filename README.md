This is orbdetpy, a Python and Java library for orbit determination.
It is a thin Python wrapper for our Java estimation tools and Orekit.

Features
--------

The force model for orbit propagation can be configured to include:

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

1. Install the Java Development Kit 8 (1.8) from <http://openjdk.java.net/install/index.html>.
   Set the `JAVA_HOME` environment variable to the JDK installation
   folder. The `java` executable must be added to the system path.

2. Install Python 3.7.0 or higher and run `pip install orbdetpy` to install
   orbdetpy and other package dependencies.

3. Update the astrodynamics data under `orbdetpy/orekit-data` periodically by
   running the following. You might need to be the `root` user on some systems.

   `python -c "from orbdetpy.astro_data import update_data; update_data();"`

Examples
--------

1. `fit_radec.py` : Run OD with real angles measurements. Also demonstrates
   the Laplace IOD method for estimating an initial state vector.

2. `predict_passes.py` : Predict satellite passes for ground stations or
   geographic regions using TLEs. Current elements may be obtained from
   sites such as <http://www.celestrak.com>.

3. `propagate_tle.py` : Propagate TLEs given by command-line arguments.

4. `test_conversion.py` : Test reference frame and other conversion functions.

5. `test_estimation.py` : Demonstrates measurement simulation and orbit
   determination functions.

6. `test_interpolation.py` : Interpolate state vectors.

Development
-----------

1. Developers will need Apache Maven 3+ to build the Java library. Build
   using the following from the `orbdetpy/` sub-folder, where `os_cpu_type` is
   `linux-x86_64`, `linux-x86_32`, `windows-x86_64`, `windows-x86_32`,
   `osx-x86_64`, or `osx-x86_32` depending on your CPU and OS:

   `mvn -e -Dos.detected.classifier=os_cpu_type package`

   The command-line is simpler on Intel/AMD 64-bit Linux:

   `mvn -e package`

2. Download and extract <https://github.com/ut-astria/orbdetpy/releases/download/2.0.3/orekit-data.tar.gz>
   under the `orbdetpy/` sub-folder.

Known Issues
------------

1. You might receive warnings from the Windows Defender Firewall on Microsoft
   Windows. Grant `orbdetpy` network access permissions.

2. If you use the `multiprocessing` Python package, imports and calls into
   `orbdetpy` must not span `multiprocessing` function calls. That is, `orbdetpy`
   can be used in the parent process or the spawned child processes, but not both.
   A workaround is to run the `orbdetpy` RPC server using `orbdetpy/start_server.sh`
   in a separate terminal window before running your Python code.
