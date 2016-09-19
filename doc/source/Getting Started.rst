.. _getting_started:

Getting Started
===============

Describe how to get and build source here.

Pre-Reqisites
-------------

Openmpi
^^^^^^^

SSAGES requires a pre-built MPI environment to achieve the best parallel
computing performance. You can use a package manager to install related packages
(openmpi-bin, openmpi-doc, libopenmpi-dev, etc. if you are using ubuntu), or you
can obtain the source code from https://www.open-mpi.org and build it
appropriately in your system. Make sure you have a newer version installed, for
example, openmpi 1.8.

Boost
^^^^^

Make sure your have newer version of boost installed (for example, boost 1.60).
You can obtain the source code and built boost appropriately following
http://www.boost.org/

Gcc
^^^

Make sure you have gcc 4.9+ installed on your system. You can install a
compatible version of gcc using a package manager, for example, apt-get in
linux, or obtain and build cmake appropriately following https://gcc.gnu.org/.

Cmake
^^^^^

SSAGES requires a Cmake build system. You can install a compatible version of
Cmake using a package manager, for example, apt-get in linux, or obtain and
build cmake appropriately following https://cmake.org. 

Python
^^^^^^

Python 2.7 is required to enable some built-in tools in SSAGES, for example,
some python scripts are prepared to generate Json files as inputs in some
method. You can use a package manager to install it or obtain the source code
from https://www.python.org/downloads/ and build it appropriately in your
system. 

Get the source code
-------------------

There are two ways of getting the source code for SSAGES: Download a ZIP file
from github or clone the git repository. We strongly recommend the second method
as it allows you to easily stay up-to-date with the latest version.

To clone the git repository, call

```
git clone https://github.com/jonathankwhitmer/SSAGES.git
```

Build SSAGES
------------

SSAGES currently works with two simulation engines: LAMMPS and Gromacs (we are
striving to add more simulation engines soon). To build SSAGES with LAMMPS,
you can use

.. code-block:: bash

    mkdir build/\n`
    cd build/
    cmake .. -DLAMMPS=YES
    make

This set of commands will automatically download LAMMPS and build it together
with SSAGES.

Run SSAGES
----------

.. todo::

    Give a short, simple example of how to start a SSAGES simulation. A more
    detailed description can be found under Tutorials.
