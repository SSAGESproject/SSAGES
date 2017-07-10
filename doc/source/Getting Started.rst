.. _getting-started:

Getting Started
===============

Prerequisites
-------------

Before you try to build SSAGES, make sure that you have the following packages
installed:

+------------+------------------+-----------------------------------+
| Package    | Required version | Package name in Ubuntu repository |
+============+==================+===================================+
| `openmpi`_ | 1.8 or higher    | openmpi-common, libopenmpi-dev    |
+------------+------------------+-----------------------------------+
| `gcc`_     | 4.9 or higher    | gcc-4.9                           |
+------------+------------------+-----------------------------------+
| `cmake`_   | 2.8 or higher    | cmake                             |
+------------+------------------+-----------------------------------+
| `python`_  | 2.7              | python2.7                         |
+------------+------------------+-----------------------------------+

.. _openmpi: https://www.open-mpi.org/
.. _gcc: https://gcc.gnu.org/
.. _cmake: https://cmake.org/
.. _python: https://www.python.org/

Get the source code
-------------------

There are two ways of getting the source code for SSAGES: Download a ZIP file
from github or clone the git repository. We strongly recommend the second method
as it allows you to easily stay up-to-date with the latest version.

To clone the git repository, call

.. code-block:: bash

    git clone https://github.com/MICCoM/SSAGES-public.git

Build SSAGES
------------

SSAGES supports a number of simulation engines including LAMMPS and Gromacs.
It is possible to have SSAGES auto-download the source codes for LAMMPS and Gromacs.
This is done by providing the option ``-DLAMMPS=YES`` for LAMMPS and ``-DGROMACS=YES`` for 
Gromacs to Cmake. However, in many cases it will be necessary to build SSAGES using 
your local copy of the MD engine source code. For example, if you have modified it 
to fit a special need LAMMPS or Gromacs does not support natively.

.. code-block:: bash

    mkdir build/
    cd build/
    cmake .. -DLAMMPS=YES
    make

or

.. code-block:: bash

    mkdir build/
    cd build/
    cmake .. -DGROMACS=YES
    make

This set of commands will automatically download LAMMPS/Gromacs and build it together
with SSAGES. Alternatively, you can :ref:`build SSAGES using a local copy of the
MD engine source code <build_ssages_with_local_md_source>`.

LAMMPS can be built with optional packages. To enable or disable these packages,
you can call

.. code-block:: bash

    make yes-package

or

.. code-block:: bash

    make no-package

respectively. For more information on optional packages for LAMMPS, refer to
the `LAMMPS user manual
<http://lammps.sandia.gov/doc/Section_start.html#making-lammps-with-optional-packages>`_.

Run SSAGES
----------

In order to run ssages you need to use run the executable followed by the input file.
For example:

.. code-block:: bash
    
    mpiexec -np 6 ./ssages input.json

Where the `-np` flag dictates the total number of processors you need and input.json
is the input file. For a description of the input file click :ref:`here <inputfiles>`.

More information on how to run SSAGES with a specific simulation engine can be found
here:

* :ref:`How to run SSAGES on LAMMPS <run_ssages_on_lammps>`.
* :ref:`How to run SSAGES with GROMACS <run_ssages_with_gromacs>`.

Advanced options
----------------

In case these simple steps do not meet your need, you can find advanced information
on building and running SSAGES here.

.. toctree::
   :maxdepth: 1

   Build SSAGES with local copy of MD source
   How to run SSAGES on LAMMPS
   How to run SSAGES with GROMACS
