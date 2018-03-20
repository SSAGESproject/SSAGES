.. _build_ssages_with_local_md_source:

Build SSAGES with local copy of MD source
=========================================

The standard procedure to build SSAGES is to auto-download the source code for the
simulation engine you intend to use. This is done by providing the option
``-DLAMMPS=YES`` (for LAMMPS) or ``-DGROMACS=YES`` (for Gromacs) to ``cmake``.
However, in many cases it will be necessary to build SSAGES using your local copy of
the MD engine source code. For example, if you have modified it to fit a special
need LAMMPS or Gromacs does not support natively.

If you want to build SSAGES using a local copy of the MD engine source code, modify
the cmake call to

.. code-block:: bash

    cmake .. -DLAMMPS_SRC=/path/to/LAMMPS/src

if you are using LAMMPS and to

.. code-block:: bash

    cmake .. -DGROMACS_SRC=/path/to/gromacs

if you are using Gromacs.

.. warning::

    The current implementation of SSAGES will patch the Gromacs source. Thus, if you
    compile the patched Gromacs source, it will no longer run. We are working to remedy
    this inconvenience.

SSAGES is not compatible with all versions of LAMMPS and Gromacs. The following
versions of LAMMPS have been tested extensively, but we are confident that SSAGES will
also work with most other LAMMPS versions.

* 10 Aug 2015
* 7 Dec 2015
* 14 May 2016
* 15 Jul 2016
* 30 Jul 2016

The following versions of Gromacs are supported

* 5.1.x
* 2016.x
* 2018

In contrast to LAMMPS, we are very confident that SSAGES will *not* work with other
versions of Gromacs out of the box. We are working hard to make SSAGES compatible with
new versions of Gromacs, as they are released.
