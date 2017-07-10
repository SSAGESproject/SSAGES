.. _engines:

Engines
=======

SSAGES supports multiple molecular dynamics engines. However, supported 
features between MD engines. This may be the result of engine limitations 
or work in progress. The table below summarizes the main features that 
vary between supported engines. 

+----------+---------------------+--------------+----------------+
| Engine   | Supported versions  | Multi-walker | NPT ensemble   |
+==========+=====================+==============+================+
| LAMMPS   |    2010 or newer    |    yes       |   yes          |
+----------+---------------------+--------------+----------------+ 
| Gromacs  |    5.1.x, 2016.3    |    yes       |   yes          |
+----------+---------------------+--------------+----------------+
| OpenMD   |      2.4, 2.5       |    no        |   no           |
+----------+---------------------+--------------+----------------+
| Qbox     |   1.60 or newer     |    no        |   no           |
+----------+---------------------+--------------+----------------+

Special instructions on how to use SSAGES with a particular engine are
listed under the appropriate section. 

LAMMPS
^^^^^^

Building
~~~~~~~~

SSAGES supports most recent versions of LAMMPS. To compile SSAGES with a 
compatible version of LAMMPS, either ``-DLAMMPS=YES`` or 
``-DLAMMPS_SRC=/path/to/LAMMPS`` must be specified in the cmake command. 
For example, 

.. code:: bash 

	cmake -DLAMMPS=YES .. 
	make

will automatically download a recent version of LAMMPS (tagged ``r15407``) 
and compile SSAGES. Because many users may take advantage of optional LAMMPS
packages, SSAGES forwards the make commands necessary to do so, such 
as 

.. code:: bash 

	make yes-molecule
	make yes-user-drude

If a user is interested in using a different version of LAMMPS or one with 
personal modifications, then SSAGES can be pointed to that particular source 
repository.

.. code:: bash 

	cmake -DLAMMPS_SRC=/path/to/lammps .. 

.. warning:: 

	Once you link SSAGES to a particular LAMMPS source, you will be 
	**unable** to compile that LAMMPS source outside of SSAGES because of 
	SSAGES dependencies which are introduced. Be sure to backup your 
	repository accordingly. 

Running 
~~~~~~~

SSAGES integrates with LAMMPS though the flexible *fix* API offered 
by LAMMPS. It is therefore necessary to define a SSAGES fix within 
the LAMMPS input file as follows.

.. code::

	fix ssages all ssages

This directive ensures that SSAGES is able to locate the appropriate 
adapter and interface with the LAMMPS library. **It is very important to 
name the fix "ssages" as shown above. Otherwise, SSAGES will not work 
properly**. It is highly recommended that the SSAGES fix command be placed 
after all integrator fixes. Also, make sure tht the fix is specified before
the run command, which will begin the advanced sampling simulation. 

.. note::

	Due to the nature of how SSAGES forwards commands to LAMMPS, the use 
	of ``include`` within a LAMMPS input script is currently not supported.

SSAGES is compatible with typical LAMMPS workflows that include equilibration 
or energy minimzation steps before production. So long as the SSAGES fix is not 
declared, LAMMPS will run without any modification. 

The only LAMMPS-specific property required in a SSAGES input file is the ``input`` 
property which points to the LAMMPS input script. Details can be found on the 
:ref:`input files page <inputfiles>`.

Gromacs
^^^^^^^

Building
~~~~~~~~

SSAGES supports Gromacs versions 5.1.x and 2016.3. To compile SSAGES with a 
compatible version of Gromacs, either ``-DGROMACS=YES`` or 
``-DGROMACS_SRC=/path/to/Gromacs`` must be specified in the cmake command. 
For example, 

.. code:: bash 

	cmake -DGROMACS=YES .. 
	make

will automatically download Gromacs 5.1.3 and compile SSAGES. 
If a user is interested in using a different version of Gromacs or one with 
personal modifications, then SSAGES can be pointed to that particular source 
repository.

.. code:: bash 

	cmake -DGROMACS_SRC=/path/to/gromacs .. 

.. warning:: 

	Once you link SSAGES to a particular Gromacs source, you will be 
	**unable** to compile that Gromacs source outside of SSAGES because of 
	SSAGES dependencies which are introduced. Be sure to backup your 
	repository accordingly. 

Running 
~~~~~~~

SSAGES forwards arguments to the Gromacs **mdrun** library. The 
``args`` property must specified in the SSAGES input file as 
described on the :ref:`input files page <inputfiles>`.

OpenMD
^^^^^^^

.. note:: 

	Coming soon. 

Qbox
^^^^^^^

.. note:: 

	Coming soon. 
