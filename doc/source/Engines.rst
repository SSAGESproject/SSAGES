.. _engines:

Engines
=======

SSAGES supports multiple molecular dynamics engines. However, supported
features between MD engines. This may be the result of engine limitations
or work in progress. The table below summarizes the main features that
vary between supported engines.

+---------------+-----------------------+--------------+------------+
| Engine        | Supported versions    | Multi-walker | NPT virial |
+===============+=======================+==============+============+
| `LAMMPS`_     | 2010 or newer         | yes          | yes        |
+---------------+-----------------------+--------------+------------+
| `GROMACS`_    | 5.1.x, 2016.x, 2018.x | yes          | yes        |
+---------------+-----------------------+--------------+------------+
| `OpenMD`_     | 2.5                   | no           | no         |
+---------------+-----------------------+--------------+------------+
| `Qbox`_       | 1.60 or newer         | yes          | no         |
+---------------+-----------------------+--------------+------------+
| `HOOMD-blue`_ | 2.4 or newer          | yes          | no         |
+---------------+-----------------------+--------------+------------+

Special instructions on how to use SSAGES with a particular engine are
listed under the appropriate section.

LAMMPS
^^^^^^

Building
~~~~~~~~

SSAGES supports most recent versions of LAMMPS. To compile SSAGES with a
compatible version of LAMMPS, either ``-DLAMMPS=YES`` or
``-DLAMMPS_SRC=/path/to/LAMMPS`` must be specified in the CMake command.
For example,

.. code:: bash

	cmake -DLAMMPS=YES ..
	make

will automatically download LAMMPS (version 30 Jul 2016, tagged ``r15407``)
and compile SSAGES. If a user is interested in using a different version of
LAMMPS, SSAGES can download a specific stable release (must be supported) with

.. code:: bash

	cmake -DLAMMPS="22 Aug 2018" ..

If a user is interested in using an already-downloaded source or one with
personal modifications, then SSAGES can be pointed to that particular source
repository.

.. code:: bash

	cmake -DLAMMPS_SRC=/path/to/lammps ..

Because many users may take advantage of optional LAMMPS
packages, SSAGES forwards the make commands necessary to do so. To enable or
disable these packages, you can call

.. code-block:: bash

	make yes-package

or

.. code-block:: bash

	make no-package

For more information on optional packages for LAMMPS, refer to
the `LAMMPS User Manual <https://lammps.sandia.gov/doc/Packages.html>`_.

.. warning::

	Once you link SSAGES to a particular LAMMPS source, you will be
	**unable** to compile that LAMMPS source outside of SSAGES because of
	SSAGES dependencies which are introduced. Be sure to backup your
	repository accordingly.

The following stable versions of LAMMPS have been tested extensively, but we
are confident that SSAGES will also work with most other LAMMPS versions.

* 10 Aug 2015
* 7 Dec 2015
* 16 Feb 2016
* 14 May 2016
* 30 Jul 2016
* 5 Nov 2016
* 17 Nov 2016
* 31 Mar 2017
* 11 Aug 2017
* 16 Mar 2018
* 22 Aug 2018

Running
~~~~~~~

SSAGES integrates with LAMMPS though the flexible *fix* API offered
by LAMMPS. It is therefore necessary to define a SSAGES fix within
the LAMMPS input file as follows.

.. code-block:: none

	fix ssages all ssages

This directive ensures that SSAGES is able to locate the appropriate
adapter and interface with the LAMMPS library. **It is very important to
name the fix "ssages" as shown above. Otherwise, SSAGES will not work
properly**. It is highly recommended that the SSAGES fix command be placed
after all integrator fixes. Also, make sure that the fix is specified before
the run command, which will begin the advanced sampling simulation.

.. note::

	Due to the nature of how SSAGES forwards commands to LAMMPS, the use
	of ``include`` and ``label/jump`` within a LAMMPS input script is
	currently not supported.

SSAGES is compatible with typical LAMMPS workflows that include equilibration
or energy minimization steps before production. So long as the SSAGES fix is not
declared, LAMMPS will run without any modification.

The only LAMMPS-specific property required in a SSAGES input file is the ``input``
property which points to the LAMMPS input script. Details can be found on the
:ref:`input files page <inputfiles>`.

GROMACS
^^^^^^^

Building
~~~~~~~~

SSAGES supports most recent versions of GROMACS. To compile SSAGES with a
compatible version of GROMACS, either ``-DGROMACS=YES`` or
``-DGROMACS_SRC=/path/to/GROMACS`` must be specified in the CMake command.
For example,

.. code:: bash

	cmake -DGROMACS=YES ..
	make

will automatically download GROMACS 5.1.3 and compile SSAGES.
If a user is interested in using a different version of GROMACS, SSAGES can
download a specific release (must be supported) with

.. code:: bash

	cmake -DGROMACS=2018.3 ..

If a user is interested in using an already-downloaded source or one with
personal modifications, then SSAGES can be pointed to that particular source
repository.

.. code:: bash

	cmake -DGROMACS_SRC=/path/to/gromacs ..

.. warning::

	Once you link SSAGES to a particular GROMACS source, you will be
	**unable** to compile that GROMACS source outside of SSAGES because of
	SSAGES dependencies which are introduced. Be sure to backup your
	repository accordingly.

The following versions of GROMACS are supported by SSAGES, but we are very
confident that SSAGES will *not* work with other versions of GROMACS out of
the box, in contrast to LAMMPS. We are working hard to make SSAGES
compatible with new versions of GROMACS, as they are released.

* 5.1.x
* 2016.x
* 2018.x

Setup
~~~~~

After compiling GROMACS with SSAGES, you can use all of GROMACS’s available
tools to set up systems and generate input files. The executable is located at
``hooks/gromacs/gromacs/bin/gmx_mpi`` within the build directory.

.. note::

	Note that, the ``gmx_mpi`` executable in the SSAGES folder will NOT
	function normally for running regular GROMACS simulations via
	``gmx_mpi mdrun``.

As GROMACS has in-depth
`Documentation <http://manual.gromacs.org/documentation/current/user-guide/>`_
and a helpful
`Getting Started <http://manual.gromacs.org/documentation/current/user-guide/getting-started.html>`_
section, we will not dwell much on how to use these tools to generate systems.

Briefly, generating a GROMACS binary input file (``.tpr``) requires the
following three files:

1. A "box" of particle coordinates to simulate (``.gro`` file)
2. A topology that describes the forcefield and connectivity (``.top`` file,
   optionally ``.itp`` files)
3. A simulation detail file that sets parameters such as which thermostat
   and barostat to use, number and length of time steps, integrator, saving
   frequency, and many more (``.mdp`` file)

For example, one can convert a protein ``.pdb`` file from an
`online database <http://www.rcsb.org/>`_ using GROMACS tools to generate
a ``.gro`` and a ``.top`` file. To generate an input file, use the preprocessor
``gmx_mpi grompp`` command:

.. code-block:: bash

	gmx_mpi grompp -f npt.mdp -p topol.top -c conf.gro -o input.tpr

There are example ``.gro``, ``.mdp``, ``.top``, ``.tpr``, and ``.json`` inputs
available in the Examples folder.

After an energy minimization and brief NVT and NPT equilibration runs, you
should be ready to use SSAGES with your system. First, generate a ``.json``
file for your SSAGES input. If using a single walker, the ``inputfile`` should
be the same as your ``.tpr`` file name. If using multiple walkers, you should
number your input files right before the extension, include a numberless file,
and set the “inputfile” to be the same as the numberless. For example, if using
four walkers, you should set your “inputfile” to ``input.tpr`` and have the
following in your folder:

* ``input.tpr``
* ``input0.tpr``
* ``input1.tpr``
* ``input2.tpr``
* ``input3.tpr``

Finally, define your CV(s) and Methods, as detailed in the
:ref:`input files <inputfiles>` page.

Running
~~~~~~~

SSAGES forwards arguments to the GROMACS **mdrun** library. The
``args`` property must specified in the SSAGES input file as
described on the :ref:`input files <inputfiles>` page.

You can start your simulation by calling the SSAGES executable:

.. code-block:: bash

	mpiexec -np N ./ssages input.json

where `N` is the total number of MPI processes. For example, for three walkers
using 2 processors each, set :math:`N = 3*2 = 6`.

OpenMD
^^^^^^

Building
~~~~~~~~

SSAGES supports most recent versions of OpenMD. To compile SSAGES with a
compatible version of OpenMD, the location of the already-downloaded source
must be specified in the CMake command.

.. code:: bash

	cmake -DOPENMD_SRC=/path/to/OpenMD ..

.. warning::

	Once you link SSAGES to a particular OpenMD source, you will be
	**unable** to compile that OpenMD source outside of SSAGES because of
	SSAGES dependencies which are introduced. Be sure to backup your
	repository accordingly.

The following versions of OpenMD are supported by SSAGES, but we are very
confident that SSAGES will *not* work with other versions of OpenMD out of
the box, in contrast to LAMMPS.

* 2.5

Running
~~~~~~~

The only OpenMD-specific property required in a SSAGES input file is the ``input``
property which points to the OpenMD input script. Details can be found on the
:ref:`input files page <inputfiles>`.

Qbox
^^^^

Building
~~~~~~~~

SSAGES and Qbox can be run together to use advanced sampling methods in
*ab initio* molecular dynamics simulations. The coupling with Qbox is performed
in a server--driver mode, with SSAGES acting as the driver and Qbox as the
server. This means that if you have access to a version of Qbox (minimum 1.60)
you do not need to recompile SSAGES and Qbox together. However, it is necessary
to configure SSAGES to be used with Qbox, so that it will compile the correct
Hook and Driver. To do so, add the following flag during the configuration of
SSAGES:

.. code:: bash

	cmake -DQBOX=YES ..

It is important to remark that in this case, **SSAGES will not automatically
download Qbox**, it will be simply configured so to communicate with it. You
are required to have access to a Qbox executable. If you do not have access to
a pre-compiled version, then you will need to
`download and compile it yourself <http://qboxcode.org/build/>`_.

Setup
~~~~~
As for other engines, there are two input scripts necessary to run a
Qbox--SSAGES calculation composed of ``N`` walkers:

1. A JSON input file, specifying the methods and CVs that you want to use.
   Also, it specifies the Qbox input file names and the number of MD, density,
   and wavefunction steps that you want to use.
2. A number ``N`` of Qbox input files, that will be used in the first step of
   the calculation to obtain the ground state density in the first step.

The JSON file contains the same field that would usually have (CVs, methods, logger, etc.) with three additional options:

.. code:: javascript

	{
		"walkers": N,
		"input": ["md.1", "md.2", ..., "md.N"],
		"md_iterations" : 10,
		"qm_iterations" : 30,
		"wf_iterations" : 1,
	}

The keywords ``walkers`` and ``input`` are the standard SSAGES keywords to
declare the number of walkers and the starting input file of each walker. The
keywords ``md_iterations``, ``qm_iterations`` and ``wf_iterations``  are the
respectively the number of MD steps to perform, the number of `scf` to perform
per MD step, and the number of wave-function optimization per `scf` steps.
These parameters correspond to the first, second and third number in the
command ``run 20 10 0``. Consult the
`Qbox Documentation <http://qboxcode.org/doc/html/>`_ for more information.

The Qbox input file of each walker, specifies the parameters to be used in the
DFT calculations (``xc``, ``ecut``, ``T``, etc.). This file will be parsed by
Qbox **at the first time step of the simulations** to set up the calculations.
If the file contains a command such as ``run 200 10``, the 200 MD steps that
Qbox will perform **will be unbiased**. If wanted, this feature can be used to
equilibrate the system. After this first step, the command
``run 1 qm_iterations wf_iterations`` will be repeated for ``md_iterations``.
**The QBox input file must contain at least 1 MD step in order to run with SSAGES.**
Thus always include the ``run 1`` command.

An example of ``input.json`` and ``md.i`` is present in the ``Examples/User/ABF/NaCl-Qbox`` directory.

Running
~~~~~~~

As previously reported, Qbox and SSAGES communicate in a server--driver mode.
To launch Qbox in a server mode is sufficient to use the proper keyword and
specify its input and output file:

.. code:: bash

	mpirun -n X qb -server ssages_in_0 ssages_out_0

for a single walker or

.. code:: bash

	mpirun -n X qb -server ssages_in_0 ssages_out_0
	mpirun -n X qb -server ssages_in_1 ssages_out_1
	....
	mpirun -n X qb -server ssages_in_N ssages_out_N

for multiple walkers. At the moment, the name ``ssages_in_`` and
``ssages_out_`` are **mandatory** and cannot be changed. When launched in this
way, Qbox creates ``N`` files called ``ssages_in_N.lock``, and then wait for
input. When the files ``ssages_in_N.lock`` are deleted from disk, Qbox will
execute the commands contained in the files ``ssages_in_N``, write the result
of the calculation in ``ssages_out_N``, and create N ``ssages_in_N.lock``
files. Without the deletion of the ``.lock`` files, Qbox will not execute any
command and will remain idle.

After Qbox has started the server mode run (so it is idling and the ``.lock``
files are present on disk), we can launch SSAGES to drive the calculations:

.. code::

	mpirun -n N ssages input.json

After SSAGES has started, the two codes will alternate with each other in the
following way:

1. SSAGES will write the script ``md.i`` to file ``ssages_in_i``, which will
   initialize the DFT parameters of the calculations. Then, it will trigger
   Qbox execution by deleting the ``.lock`` files.
2. Qbox will perform the DFT calculation specified in ``ssages_in_i`` and write
   the output in ``ssages_out_i`` and will recreate the ``.lock`` files.
3. SSAGES will read the Qbox output, calculate the CVs and the bias, and write
   to the files ``ssages_in_i``, containing the external forces and the
   position of the atoms, as well as the command
   ``run 1 qm_iterations wf_iterations``. It will then delete the ``.lock``
   file, triggering another MD step calculation in Qbox.
4. Steps 2 and 3 will be repeated for ``md_iterations`` number of time.
5. After the last iteration, SSAGES will write an input file that will instruct
   Qbox to save a ``restart_i.xml`` file that can be used to restart the
   calculations, and terminate the Qbox instance.
6. Qbox and SSAGES will then finish the execution.

Normally, Qbox overwrites the output ``ssages_out_i`` in server mode. To
preserve the trajectory and avoid the loss of data, SSAGES will append the
``ssages_out_i`` file to a ``ssages_out_i_run_j.xml`` file. In the latter, the
``i`` index identifies the walker, while the ``j`` index identifies the number
of runs. (For example, if you restarted two times, you would have
``_run_1.xml``, ``_run_2.xml``, and ``_run_3.xml``.) We suggest using the
``restart_i.xml`` files to avoid discontinuities in the trajectories; when
restarting, create a ``md.i`` file that contains the ``load restart_i.xml``
instructions.

There are useful scripts to analyze and plot Qbox trajectories, which are available in the `Qbox tools webpage <http://qboxcode.org/tools//>`_. To run any of these scripts, first reformat ``ssages_out_i_run_j.xml`` file by running a python script ``Qbox-xml-cleaning.py`` present in ``Tools/`` directory. For example, if the original ssages-qbox output file is ``ssages_out_0_run_0.xml`` then the command line to reformat this xml file is

.. code:: bash

	python3 Qbox-xml-cleaning.py ssages_out_0_run_0.xml ssages_out_0_run_0_cleaned.xml

where first arugment the name of the original xml output file, and the second argument is the reformatted xml file. Now these files can be analyzed using the scripts in `Qbox tools webpage <http://qboxcode.org/tools//>`_. For example, to create xyz trajectory file from the reformatted output  ``ssages_out_0_run_0_cleaned.xml``, run the command

.. code:: bash

	python2 qbox_xyz.py -all ssages_out_0_run_0_cleaned.xml > out_0_run_0.xyz
	
Running on Clusters
~~~~~~~~~~~~~~~~~~~

Most likely, you are going to launch this calculation on a cluster or a
supercomputer, where you will need to prepare a submission scripts and then
launch through a job scheduler. Given the fact that SSAGES needs to start
*after* Qbox, it is better to either separate the scripts that submit the two
different calculations, or use a syntax that ensure that the submission occurs
in the right order. For example, on `Slurm <https://slurm.schedmd.com/>`_, we
can use one script:

.. code:: bash

	srun -n X  -N 1 qb -server ssages_in0 ssages_out0  &
	srun -n X  -N 1 qb -server ssages_in1 ssages_out1  &
	srun -n 2  -N 1 ssages input.json &
	wait

which ensures that the scripts are executed in the right way.

If you want to have different scripts for Qbox and SSAGES:

In the Qbox scripts, ``qb.sh``

.. code:: bash

	srun -n X  -N 1 qb -server ssages_in0 ssages_out0
	srun -n X  -N 1 qb -server ssages_in1 ssages_out1

In the SSAGES script, ``ssages.sh``

.. code:: bash

	srun -n 2  -N 1 ssages input.json

Then you will need to submit both of them with a third script, ``launch.sh``

.. code:: bash

	#!/bin/bash

	j_qb=`sbatch qb.sh | awk '{print $4}'`

	sbatch --dependency=after:${j_qb} ssages.sh


The advantage of the latter method, with three scripts, is that it will avoid
conflict between modules, which may be present in the first example, depending
on how you have compiled Qbox and SSAGES.

HOOMD-blue
^^^^^^^^^^

Building
~~~~~~~~

HOOMD-blue supports SSAGES on the current ``master`` branch as of 2018-08-01,
and will be supported in releases v2.4.0 and later. HOOMD-blue must be built
with **MPI support enabled**, using the HOOMD CMake flag ``-DMPI_ENABLED=ON``.
If building HOOMD from source, there are two ways to create a valid
installation, which contains the ``include`` files needed by SSAGES.

1. HOOMD can be installed with ``make install`` to the location specified by
   the HOOMD CMake flag ``-DCMAKE_INSTALL_PREFIX``.
2. HOOMD can be built using ``make`` (without installation) using the HOOMD
   CMake flag ``-DCOPY_HEADERS=ON``.

For more information on HOOMD CMake flags, see `the HOOMD documentation <https://hoomd-blue.readthedocs.io/en/stable/compiling.html>`_.
After installing HOOMD, specify the SSAGES CMake flag ``-DHOOMD_ROOT=/path/to/hoomd_installation/hoomd``, pointing to the ``hoomd`` directory within the HOOMD-blue installation prefix.
For example (using 8 ``make`` job slots),

.. code:: bash

        cd /path/to/hoomd
        mkdir build
        cd build/
        cmake ../ -DENABLE_MPI=on -DCMAKE_INSTALL_PREFIX=/path/to/hoomd_installation
        make install -j 8

        cd /path/to/SSAGES
        mkdir build
        cd build/
        # Note that this path adds a suffix .../hoomd
        cmake ../ -DHOOMD_ROOT=/path/to/hoomd_installation/hoomd
        make -j 8

Running
~~~~~~~

HOOMD-blue offers a "half-step hook" API to support SSAGES.
This feature is automatically configured when SSAGES launches
the HOOMD-blue simulation.

The HOOMD-blue SSAGES user script is written in Python. This script should
contain necessary ``import`` statements, configure the simulation, and set
types, interactions, integrators, log outputs and so forth. However, the
simulation user script should *not* call ``hoomd.run(steps)`` as this will
be called within SSAGES.

To set
`HOOMD-blue command-line options <http://hoomd-blue.readthedocs.io/en/stable/command-line-options.html>`_,
use the SSAGES JSON input file. Set the key ``"args"`` with a string of command
line options that will be passed to ``hoomd.context.initialize()``.
