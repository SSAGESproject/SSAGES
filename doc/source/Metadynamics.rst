.. _metadynamics:

Metadynamics
--------------------

Introduction
^^^^^^^^^^^^

.. todo::

    Short introduction to Metadynamics.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Metadynamics is selected by defining ``type" : "Metadynamics"`` as the 
method in the JSON input file. It supports the following options:

widths 
   *array of doubles (number of CVs) long*.
   This array defines the width of the Gaussians being depositied over time
   in each CV dimension.

height 
	*double*. 
	This value defines the height of the Gaussians being deposited over time. 

hill_frequency 
	*double* 
	This value defines the frequency in iterations with which the Gaussians 
	are deposited. 

.. note::

	The Metadynamics method does not yet support CV bounds.

.. note::
	
	The Metadynamics method does not yet support restarts.

.. warning::

	The Metadynamics will run for the duration specified in the 
	input file under ``MDSteps``. It is the user's responsibility to ensure that 
	enough time is given for Metadynamics to obtain an satisfactory representation
	of the free energy surface of interest. A simple way to prevent early
	quitting is to define a very large number of ``MDSteps``` and to periodically
	check the output file for convergence. 

Example Input 
^^^^^^^^^^^^^

.. code-block:: javascript 

	"method" : {
		"type" : "Metadynamics", 
		"widths" : [0.1, 0.1],
		"height" : 0.1,
		"hill_frequency" : 1
	}

Output
^^^^^^

The output of Metadynamics is stored in a file called "hills.out". This file 
contains the Gaussians that have been deposited over time along the CVs. 
Each time Gaussian is deposited, a new line is written to the file in the following
format: 

*center1 center2 ... width1 width2 ... height* 

The centers represent the location in CV space where the Gaussian has been 
deposited. The widths represent the corresponding Gaussian widths for each 
CV dimension and the height is the Gaussian height. 

.. note:: 

	Although the widths and height of the Gaussian currently do not change in
	time, future additions to the Metadynamics method will allow for adaptive 
	Gaussians.

In the ``Examples/User/Meta`` directory, example MATLAB scripts are provided 
that sum the Gaussians and generate a free energy surface from the "hills.out"
file.



Tutorial
^^^^^^^^

Two Metadynamics examples are included in the ``Examples/User/Meta`` directory. 
In the first example, the free energy surface of two-dimensional Langevin 
particle (``Single_Atom`` folder) is sampled. This example requires LAMMPS.
The files included are described below: 

* ``in.LAMMPS_Meta_Test`` - LAMMPS input file describing the Langevin particle 
  and underlying free energy surface to be sampled.
* ``Meta.json`` - SSAGES JSON input file specifying Metadynamics and CVs to be 
  sampled. In this case the CVs are the *x* and *y* coordinates of the particle. 
* ``analysis.m`` - MATLAB script that analyzes the output of the Metadynamics 
  method. 
* ``Movie.m`` - MATLAB script that generates a movie of the free energy 
  surface estimate over time.

To run this example:

1) Either copy or create a symbolic link to the SSAGES executable in the
   examples directory. 

.. code-block:: bash 

	ln -s /path/to/SSAGES/build/ssages 

2) Run the example by issuing the command below. Please note that in this
   example, two walkers are used to explore the system more efficiently. If 
   you would like to use more walkers (1 processor per walker), simply include
   more drivers in the ``Meta.json`` input file. 

.. code-block:: bash 

	mpirun -np 2 ./ssages Meta.json 

3) After the run is complete use the provided ``analysis.m`` script to generate 
   a representation of the underlying free energy surface.


Developer
^^^^^^^^^

Hythem Sidky.

