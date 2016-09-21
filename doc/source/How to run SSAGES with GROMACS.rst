How to use SSAGES with GROMACS
--------------------------------


You can compile SSAGES with GROMACS using cmake by either pointing to a GROMACS source via -DGROMACS_SRC=/path/to/gromacs, 
or use the auto-download feature with -DGROMACS=ON. Please note that current implementation of SSAGES will patch the GROMACS source, 
and thus the source GROMACS will not run normally.

After compiling GROMACS with SSAGES, you can use all of GROMACS’ available tools to set up systems and generate input files. 
The executable (gmx_mpi) is located in ssages build/gromacs/bin.

As GROMACS has an in-depth documentation and getting started section, we will not dwell much on how to use these tools to generate systems. 
Please see:
`here <http://manual.gromacs.org/online/getting_started.html>`_.
and
`here <http://www.gromacs.org/Documentation>`_.
for details.


Briefly, a GROMACS input file (.tpr) requires the following three to generate:
1) A ‘box’ of particles to simulate (.gro file)
2) A topology that describes the forcefield and connectivity (.top file, optionally .itp files)
3) A simulation details file that sets many parameters such as which thermostat and barostat to use if any, timesteps, integrator, saving frequency and many more (.mdp file)
For example, one can convert a protein .pdb file from an online database using GROMACS tools to generate a .gro and a .top file. To generate an input file, use the gmx_mpi grompp command:

.. code-block:: bash
	
	gmx_mpi grompp -f npt.mdp -p topol.top -c conf.gro -o input.tpr

Please note again that currently, the gmx_mpi executable in the SSAGES folder will NOT function normally for running regular GROMACS simulations via gmx_mpi mdrun. 

After an energy minimization and brief nvt and npt equilibration runs, you should be ready to use SSAGES with your system. First, generate a .json for your SSAGES input. 
If using a single walker, the “inputfile” should be the same as your .tpr file name. If using multiple walkers, you should number your input files right before the extension, 
include a numberless file, and set the “inputfile” to be the same as the numberless. For example, if using four walkers, you should set your “inputfile” to “input.tpr” 
and have the following in your folder:

* ``input.tpr`` 
* ``input0.tpr`` 
* ``input1.tpr`` 
* ``input2.tpr``
* ``input3.tpr``

The numberless input.tpr will not be used. Then, for each walker, set the “type” to “Gromacs”, and define the number of MPI walkers to use for each walker with “number processors”. 
Finally, define your CV(s) and Methods, either generally or for each walker. You can start your simulation by calling the ssages executable:

.. code-block:: bash

	mpirun -np X ./ssages input.json

where X = total number of MPI processes. For example, for three walkers with “number processors” : 2,  X = 3*2 = 6.

Normally, you can also define an observer in .json to automatically generate backups that will save both simulation snapshots as well as method-critical data. 
However, this feature is not yet implemented for GROMACS.

There are example .gro, .mdp, .top, .tpr and .json inputs available in the Examples folder.
