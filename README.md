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
   <http://openjdk.java.net/install/index.html>. Set the `JAVA_HOME`
   environment variable to point to your JDK installation. The `java`
   executable must also be in your system path.

2. Install Python 3.7.0 or higher and run `pip install orbdetpy` to install
   orbdetpy and other package dependencies.

3. Update the astrodynamics data under `orbdetpy/orekit-data` periodically by
   running the following. You might need to be the `root` user on some systems.

   `python -c "from orbdetpy.astro_data import update_data; update_data();"`

Examples
--------

1. `predict_passes.py` : Predict satellite passes for ground stations or
   geographic regions using TLEs. Current elements may be obtained from
   sites such as <http://www.celestrak.com>.

2. `propagate_tle.py` : Propagate TLEs given by command-line arguments.

3. `test_conversion.py` : Test reference frame and other conversion functions.

4. `test_estimation.py` : Demonstrates measurement simulation and orbit
   determination functions.

5. `test_interpolation.py` : Interpolate state vectors.

Development
-----------

Developers will need Apache Maven 3+ to build the Java library. Build
using the following from the `orbdetpy/` sub-folder, where `os_cpu_type` is
`linux-x86_64`, `linux-x86_32`, `windows-x86_64`, `windows-x86_32`,
`osx-x86_64`, or `osx-x86_32` depending on your CPU and OS:

`mvn -e -Dos.detected.classifier=os_cpu_type package`

The command-line is simpler on Intel/AMD 64-bit Linux:

`mvn -e package`

Download and extract <https://github.com/ut-astria/orbdetpy/releases/download/2.0.0/orekit-data.tar.gz>
under the `orbdetpy/` sub-folder.

Known Issues
------------

1. In Microsoft Windows, you might receive warnings from the "Windows
   Defender Firewall" when you import `orbdetpy`. You can safely allow
   `orbdetpy` access to the network because this only involves a single
   TCP port on `localhost`.
