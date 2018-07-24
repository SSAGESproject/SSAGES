.. _engines:

Engines
=======

SSAGES supports multiple molecular dynamics engines. However, supported 
features between MD engines. This may be the result of engine limitations 
or work in progress. The table below summarizes the main features that 
vary between supported engines. 

+----------+-----------------------+--------------+----------------+
| Engine   |  Supported versions   | Multi-walker | NPT virial     |
+==========+=======================+==============+================+
| LAMMPS   |     2010 or newer     |    yes       |   yes          |
+----------+-----------------------+--------------+----------------+ 
| GROMACS  | 5.1.x, 2016.x, 2018.x |    yes       |   yes          |
+----------+-----------------------+--------------+----------------+
| OpenMD   |          2.5          |    no        |   no           |
+----------+-----------------------+--------------+----------------+
| Qbox     |     1.60 or newer     |    yes       |   no           |
+----------+-----------------------+--------------+----------------+

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

will automatically download LAMMPS (version 30 Jul 2016, tagged ``r15407``)
and compile SSAGES. Because many users may take advantage of optional LAMMPS
packages, SSAGES forwards the make commands necessary to do so, such 
as 

.. code:: bash 

	make yes-molecule
	make yes-user-drude

If a user is interested in using a different version of LAMMPS, SSAGES can
download a specific stable release (must be supported) with

.. code:: bash

	cmake -DLAMMPS="5 Nov 2016" ..

If a user is interested in using an already-downloaded source or one with
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

GROMACS
^^^^^^^

Building
~~~~~~~~

SSAGES supports GROMACS versions 5.1.x, 2016.x, and 2018.x. To compile SSAGES
with a compatible version of GROMACS, either ``-DGROMACS=YES`` or 
``-DGROMACS_SRC=/path/to/GROMACS`` must be specified in the cmake command. 
For example, 

.. code:: bash 

	cmake -DGROMACS=YES .. 
	make

will automatically download GROMACS 5.1.3 and compile SSAGES. 
If a user is interested in using a different version of GROMACS, SSAGES can
download a specific release (must be supported) with

.. code:: bash

	cmake -DGROMACS=2016.4 ..

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

Running 
~~~~~~~

SSAGES forwards arguments to the GROMACS **mdrun** library. The 
``args`` property must specified in the SSAGES input file as 
described on the :ref:`input files page <inputfiles>`.

OpenMD
^^^^^^^

Building
~~~~~~~~

SSAGES supports OpenMD version 2.5. To compile SSAGES with a compatible 
version of OpenMD, the location of the already-downloaded source must be 
specified in the cmake command.

.. code:: bash

    cmake -DOPENMD_SRC=/path/to/OpenMD ..

.. warning::

	Once you link SSAGES to a particular OpenMD source, you will be 
	**unable** to compile that OpenMD source outside of SSAGES because of  
	SSAGES dependencies which are introduced. Be sure to backup your 
	repository accordingly.

Running
~~~~~~~

The only OpenMD-specific property required in a SSAGES input file is the ``input`` 
property which points to the OpenMD input script. Details can be found on the 
:ref:`input files page <inputfiles>`.

Qbox
^^^^

Building
~~~~~~~~

Qbox and SSAGES can be used together to use enhanced sampling methods in *ab initio* molecular dynamics simulations. The coupling with Qbox is performed in a server-driver mode, with SSAGES acting as the driver and Qbox as the server. This means that if you have access to a version of Qbox (minimum 1.60) you do not need to recompile SSAGES and Qbox together. However, it is necessary to configure SSAGES to be used with Qbox, so that it will compile the correct Hook and Driver. To do so, add the following flag during the configuration of SSAGES

.. code:: bash

        cmake -DQBOX=YES .. 

It is important to remark that in this case, **SSAGES will not automatically download Qbox**, it will be simply configured so to communicate with it. You are required to have access to a Qbox executable. If you do not have access to a precompiled version, then you will need to download and compile it yourself (http://qboxcode.org/build/)

Set-up
~~~~~~
As for other engines, there are two input scripts necessary to run a Qbox-SSAGES calculation composed of `N` walkers:

1. A JSON input file, specifying the methods and CVs that you want to use. Also, it specifies the qbox input file names and the number of MD, density, and wavefunction steps that you want to use.
2. A number N of Qbox input files, that will be used in the first step of the calculation to obtain the ground state density in the first step.

The JSON file contains the same field that would usually have (CVs, methods, logger, etc..) with three more options:

.. code:: javascript

       {
           "walkers": N,
           "input": "[md.1,md.2,..md.N]",
           "md_iterations" : 10,
           "qm_iterations" : 30,
           "wf_iterations" : 1,
       }

The keywords ``walkers`` and ``input`` are the standard SSAGES keywords to declare the number of walkers and the starting input file of each walker. The keywords ``md_iterations``, ``qm_iterations`` and ``wf_iterations``  are the respectively the number of MD steps to perform, the number of `scf` to perform per MD step, and the number of wave-function optimization per `scf` steps. These parameters correspond to the first, second and third number in the command `run 20 10 0`  (http://eslab.ucdavis.edu/software/qbox/QboxUserGuide.pdf).

The Qbox input file of each walker, specifies the parameters to be used in the DFT calculations (`xc`,`ecut`, T etc..). This file will be parsed by Qbox **at the first timestep of the simulations** to set up the calculations. If the file contains a command such as `run 200 10` the 200 MD steps that Qbox will perform **will be unbiased**. If wanted, this feature can be used to equilibrate the system. After this first step, the command `run 1 qm_iterations wf_iterations` will be repeated for `md_iterations`. 

An example of `input.json` and `md.i` is present in the `/Examples/User/ABF/NaCl-Qbox` folder.

Running
~~~~~~~

As previously reported, Qbox and SSAGES communicate in a server-driver mode. To launch Qbox in a server mode is sufficient to use the proper keyword and specify its input and output file:

.. code:: bash
       
       mpirun -n X qb -server ssages_in_0 ssages_out_0

for a single walker or 

.. code:: bash

       mpirun -n X qb -server ssages_in_0 ssages_out_0
       mpirun -n X qb -server ssages_in_1 ssages_out_1
       ....
       mpirun -n X qb -server ssages_in_N ssages_out_N

for multiple walkers. At the moment, the name `ssages_in_` and `ssages_out_` are **mandatory** and cannot be changed. When launched in this way, Qbox creates N files called `ssages_in_N.lock`, and then wait for input. When the files `ssages_in_N.lock` are deleted from disk, Qbox will execute the commands contained in the files `ssages_in_N`, write the result of the calculation in `ssages_out_N`, and create N `ssages_in_N.lock` files. Without the deletion of the `.lock` files, Qbox will not execute any command and will remain idle.

*After* that Qbox has *started* the server mode run (so it is idling and the file `.lock` are present on disk) we can launch SSAGES to drive the calculations: 


.. code::

       mpirun -n N ssages input.json

After that SSAGES started, the two codes will alternate each other in the following way:

1. SSAGES will write on file `ssages_in_i` the script `md.i`, that will initialize the DFT parameters of the calculations. Then, it will trigger Qbox execution by deleting the .lock files.
2. Qbox will perform the DFT calculation specified in `ssages_in_i`, write the output in `ssages_out_i` and will recreate the `.lock` file.
3. SSAGES will read the Qbox output, calculate the CVs values and the bias, and write the files `ssages_in_i` containing the external forces and the position of the atoms, as well as the command `run 1 qm_iterations wf_iterations`. It will then delete the .lock file, triggering another MD step calculation in Qbox.
4. Step 2 and 3 will be repeated for `md_iterations` number of time.
5. After the last iterations, SSAGES will write an input file that will instruct Qbox to save a `restart_i.xml` file that can be used to restart the calculations, and `quit` to terminate the Qbox instance.
6. Qbox and SSAGES will then finish the execution.

Normally, Qbox overwrites the output `ssages_out_i` in server mode. To preserve the trajectory and avoid the loss of data, SSAGES will append the `ssages_out_i` file to a `ssages_out_i_run_j.xml` file. In the latter, the `i` index identify the walker, while the `j` index identifies the number of runs (So if you restarted two times, you would have _run1.xml, _run2.xml, and _run3.xml). We suggest using the `restart_i.xml` files to avoid discontinuities in the trajectories: when restarting, create a `md.i` file that contains the `load restart_i.xml` instruction.

Running on clusters
~~~~~~~~~~~~~~~~~~~

Most likely, you are going to launch this calculation on a cluster or a supercomputer, where you will need to prepare a submission scripts and then launch it through a job scheduler. Given the fact that SSAGES need to start *after* Qbox, it is better to either separate the scripts that submit the two different calculations, or use a syntax that ensure that the submission occurs in the right order. For example on *slurm* we can either use one script 

.. code:: bash 

       srun -n X  -N 1 qb -server ssages_in0 ssages_out0  &
       srun -n X  -N 1 qb -server ssages_in1 ssages_out1  &
       srun -n 2  -N 1 ssages input.json &
       wait

which ensure that the script are executed in the right way. If you want to have different scripts for Qbox and SSAGES:

In the qbox scripts, `qb.sh`

.. code:: bash

       srun -n X  -N 1 qb -server ssages_in0 ssages_out0  
       srun -n X  -N 1 qb -server ssages_in1 ssages_out1  

In the SSAGES script, `ssages.sh`

.. code:: bash 

       srun -n 2  -N 1 ssages input.json 

Then you will need to submit both of them with a third script, `launch.sh`

.. code:: bash 

       #!/bin/bash
       
       j_qb=`sbatch qb.sh | awk '{print $4}'`
       
       sbatch --dependency=after:${j_qb} ssages.sh


The advantage of the latter method, with three scripts, is that it will avoid conflict between modules, which may be present in the first example, depending on how you have compiled Qbox and SSAGES.

