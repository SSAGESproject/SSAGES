.. _umbrella-sampling:

Umbrella Sampling
-----------------

Introduction
^^^^^^^^^^^^

Calculations of thermodynamic data and other properties rely on proper sampling of the configurational space. 
However, the presence of energy barriers can lead to configurations not being sampled properly or sampled at 
all. Umbrella sampling is a simulation technique that helps to overcome those barriers and improve sampling 
by applying a bias along a coordinate. The bias takes the form of a harmonic potential.  Usually, a series of 
umbrella-sampled simulations are performed and analyzed together using the weighted histogram analysis method 
(WHAM).

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

The following parameters need to be set under `"method"` in the JSON input file:

.. code-block:: javascript

    "type" : "Umbrella"
	
The following options are available for Umbrella Sampling:

centers (required)
	Array of target values for each CV

ksprings (required)
	Array of spring constants for each CV 

file name (required)
	Output file name for umbrella sampling data

log every (optional)
	Frequency of writing to output file (default: 1)

The umbrella sampling method can also be run with time dependent centers. 
This is equivalent to the "steered MD" method. To do so, ommit `"centers"` 
and use the following options instead: 

centers0 (required)
	Array of initial target values for each CV 

centers1 (required)
	Array of final target values for each CV 

timesteps (required)
	The number of timesteps over which to scale the umbrella centers 

By setting the above options, the umbrella method will linearly interpolate 
between `centers0` and `centers1` in `timesteps` iterations. 

.. _Umbrella_tutorial:

Tutorial
^^^^^^^^

This tutorial will go through running Umbrella Sampling on an atomistic model of butane using LAMMPS as the 
MD engine. Umbrella sampling will be performed on the torsional CV of the butane C atoms. The files that can 
be found in ``Examples/User/Umbrella`` are:

``Butane_SSAGES.in``
	LAMMPS input file

``Butane.data``
	LAMMPS data file describing butane molecule.

``Template_Input.json``
	Template JSON input containing information for one Umbrella Sampling simulation. 

``Umbrella_Input_Generator.py``
	Python script for creating Umbrella.json input file. The total number of simulations and the 'centers' values 
	are controlled in this file.
	
Once in the directory, the appropriate .json file needs to be generated. A .json file is already in the directory,
``Template_Input.json``, which contains the CV information and specifies the LAMMPS input files to be used. Using 

.. code-block:: bash

	python Umbrella_Input_Generator.py

will generate Umbrella.json. Umbrella.json contains the information from Template_Input.json duplicated 12 
times with varying values of ‘centers’. These values correspond to the target values of the torsional CV. 

To run SSAGES do:

.. code-block:: bash 

	mpiexec -np 12 /path/to/SSAGES/build/.ssages Umbrella.json
	
This will run 12 different umbrella sampling simulations, one per processor. 12 different output files 
will be generated, each containing the iteration, target value of the corresponding 'center' CV, 
and the value of the CV at the iteration number. 

