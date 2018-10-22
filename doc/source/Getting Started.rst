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

Get the Source Code
-------------------

There are two ways of getting the source code for SSAGES: Download a ZIP file
from GitHub or clone the Git repository. We strongly recommend the second
method, as it allows you to easily stay up-to-date with the latest version.

To clone the Git repository, call

.. code-block:: bash

    git clone https://github.com/MICCoM/SSAGES-public.git

Build SSAGES
------------

SSAGES supports a number of simulation engines. If you have already downloaded
an engine's source files, the general steps for building SSAGES with this
engine are as follows. (The example here assumes the engine is called "Engine"
and uses the respective CMake flag ``-DENGINE_SRC``.)

.. code-block:: bash

    mkdir build/
    cd build/
    cmake .. -DENGINE_SRC=/path/to/engine
    make

The different supported engines have respective ``_SRC`` flags to specify
their source code path. For engine-specific steps and options, consult the
:ref:`Engines <engines>` page, which provides more in-depth directions. SSAGES
supports most CMake flags, including ``-DCMAKE_BUILD_TYPE`` or specifying the
compilers with ``-DCMAKE_C_COMPILER=`` and ``-DCMAKE_CXX_COMPILER=``.

It is possible to have SSAGES auto-download the source codes for LAMMPS and
GROMACS. This is done by providing the option ``-DLAMMPS=YES`` for LAMMPS and
``-DGROMACS=YES`` for GROMACS to CMake. However, in many cases, it will be
necessary to build SSAGES using your local copy of the MD engine source code.
For example, if you have modified it to fit a special need, LAMMPS or GROMACS
does not support natively.

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

This set of commands will automatically download LAMMPS/GROMACS and build it together
with SSAGES.

For these and other supported engines, you can build with a local copy of the
MD engine source code, as described in the :ref:`Engines <engines>` section.
This can be useful if you are building on a system that has limited or no
connectivity or if you want to add user modifications to the engine's
source code before compiling.

Run SSAGES
----------

In order to run SSAGES, you call the executable followed by the input file.
For example, with an input file called input.json, simple single-core jobs
can call

.. code-block:: bash

    ./ssages input.json

while jobs running on multiple threads can call

.. code-block:: bash

    mpiexec -np 6 ./ssages input.json

Here, the ``-np`` flag dictates the total number of processors on which the
simulation will run and input.json is the input file. For more information,
consult the :ref:`Input Files <inputfiles>` section.

Advanced Options
----------------

In case these simple steps do not meet your need, you can find engine-specific
options and advanced information on building and running SSAGES in the
:ref:`Engines <engines>` section.
