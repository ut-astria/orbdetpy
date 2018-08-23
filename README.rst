========================================
orbdetpy - Orbit determination in Python
========================================

Introduction
------------

This is orbdetpy, a library of Python routines for orbit determination.
It is built on top of the Orekit astrodynamics framework.

orbdetpy is free software, distributed under the terms of the `GNU
General Public License <http://www.gnu.org/licenses/gpl.html>`_.

Features
--------

The force model for orbit propagation currently includes:

1) EGM96 gravity field up to degree and order 360.
2) Solid tides owing to the Sun and Moon.
3) FES 2004 ocean tide model up to degree and order 100.
4) Simple exponential model for atmospheric drag.
5) Solar radiation pressure.
6) Third body perturbations from the Sun and Moon.

The measurement model supports range and range-rate readings. Filtering
is done using Orekit's Extended Kalman Filter or our Unscented Kalman
Filter.

You can either use your own measurements or simulate observations using
the simdata.py module.

Future Work
-----------

The following tasks are under consideration. Community contributions are
always welcome.

1) The addition of other atmospheric models supported by Orekit.
2) Support for attitude dependent drag and solar radiation pressure.
3) Other measurement types including angles data.
4) A batch least squares implementation.
5) Rauch-Tung-Striebel smoother.
6) Parametric analysis i.e. the ability to pass-through certain
   measurement types.

Prerequisites
-------------

1) Python 3.6+ must be installed with the packages numpy, scipy, pyjnius,
   matplotlib.
2) Install the Java Development Kit 8+ (1.8+) from `here
   <http://openjdk.java.net/>`_. Set the JAVA_HOME environment variable
   to point to your JDK installation.
3) `Hipparchus 1.3+ <https://hipparchus.org/>`_ and `Orekit 9.2+
   <https://www.orekit.org/>`_ are needed for astrodynamics functions.
   As a convenience, the JAR files and associated data for these
   libraries are provided under the lib/ and data/ folders, respectively.

Configuration
-------------

orbdetpy uses JSON files to store settings, measurements and estimation
results. The file examples/data/example_cfg.json is a sample, which must
be changed to describe your own orbit determination problems.

Examples
--------

The following example programs can be found in the 'examples' folder.

1) testsim.py : Demonstrates the measurement simulator for the Galaxy 15
   satellite. examples/data/testsim.json is used to configure the simulation.
   Note that maneuvers can be incorporated into the force model during
   a simulation run.

2) plotsim.py : Plots the results of simulations created using testsim.py.

3) testodet.py : Demonstrates orbit determination in orbdetpy. The file
   examples/data/example_cfg.json contains settings that drive
   propagation and estimation. testsim.py was used to simulate the
   measurements in examples/data/example_input.json for Galaxy 15.

4) plotodet.py : Plots the results of fitting orbits using testodet.py.

Bug Reports
-----------

Comments, criticisms and bug reports are very welcome and may be sent to
the package maintainer by email or the project's website.

Shiva Iyer <shiva.iyer AT utexas DOT edu>
