.. _metadynamics:

Metadynamics
--------------------

Introduction
^^^^^^^^^^^^

Metadynamics defines a class of flat-histogram methods useful for
molecular dynamics simulations. Within metadynamics, a
history-dependent bias is accrued through the periodic application of
elemental Gaussian biases to the collective variables (CVs) of
interest. The form of the individual biases is

.. math:: g(\vec{\xi},\vec{s}_i) = W_i
   e^{-\frac{\left[\vec{\xi}-\vec{s}_i\right]^2}{2\sigma_i^2}}\;.

Here :math:`\vec{xi}` is the collective variable, :math:`\vec{s}_i` is
the location of the :math:`i^{th}` hill, :math:`W_i` is the weight and
:math:`\sigma_i` the width of the hill. This bias acts to push the
system away from previously visited states. As this bias accrues to
the height of nearby features in the free energy surface, the system's
trajectory will begin to explore an expanded area in CV space. After a
sufficient number of biases have been applied, the total applied bias
within a given region will begin to oscillate around the negative of
the free energy in that region.

The free energy surface (FES) at any time :math:`t` may be
reconstructed by summing over all applied biases thusly

.. math:: F(\xi,t) = -V(\vec{\xi}) = -\sum_{i<n(t)} W_i
	  e^{-\frac{\left[\vec{\xi}-\vec{s}_i\right]^2}{2\sigma_i^2}}\;.

Here, :math:`n(t)` refers to the number of biases applied before time
:math:`t`. The time-dependent FES can be used to determine whether or
not a simulation has reached a converged FES by computing block
averages of the applied bias (see [Singh2012]_).

Metadynamics can be applied to unbounded regions of CV space, or to
bounded regions. For the case of bounded, non-periodic CVs, boundary
corrections must be applied which alter the structure of the hills
:math:`g(\vec{\xi},\vec{s}_i)` (see [McGovern2013]_)

Metadynamics exists in many flavors, in which :math:`W_i` and
:math:`\sigma_i` are altered in a time- or trajectory- dependent
fashion. These include Well-Tempered Metadynamics, Transition-Tempered
Metadynamics and Metadynamics with Adaptive Gaussians. Each of these
methods has advantages and drawbacks. Currently, SSAGES includes only
standard Metadynamics with fixed-shape hills and no tempering
algorithm. Details on how to use these algorithms within SSAGES are
given in the following sections.
  
Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Metadynamics is selected by defining ``type" : "Metadynamics"`` as the 
method in the JSON input file. It supports the following options:

widths 
   *array of doubles (length: number of CVs)*.
   This array defines the width of the Gaussians being deposited over time
   in each CV dimension.

height 
	*double*. 
	This value defines the height of the Gaussians being deposited over time,
	in units of energy used by the MD engine. 

hill_frequency 
	*double* 
	This value defines the frequency in iterations with which the Gaussians 
	are deposited. 

.. note::

	The Metadynamics method does not yet support CV bounds.

.. note::
	
	The Metadynamics method does not yet support restarts.

.. warning::

	Metadynamics will run for the duration specified in the 
	input file under ``MDSteps``. It is the user's responsibility to ensure that 
	enough time is given for Metadynamics to obtain an satisfactory representation
	of the free energy surface of interest. A simple way to prevent 
	Metadynamics from terminating before convergence is to define a 
	very large number of ``MDSteps``` and to periodically
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
contains the location, width and height of each Gaussian that has been deposited
over time. Each time Gaussian is deposited, a new line is written to the file 
in the following format: 

*center1 center2 ... width1 width2 ... height* 

The centers denote the locations in CV space where each Gaussian has been deposited; 
these are listed in the order that the CVs appear in the SSAGES JSON input file. 
The widths denote the corresponding Gaussian widths for each CV dimension. 
The height is the Gaussian height, which should match the parameter height defined 
in the JSON input file.

.. note:: 

	Although the widths and height of the Gaussian currently do not change in
	time, future additions to the Metadynamics method will allow for adaptive 
	Gaussians.

Example MATLAB scripts are provided in the Examples/User/Meta directory. 
These scripts sum the Gaussians and generate a free energy surface from the "hills.out" 
file.


Tutorial
^^^^^^^^

Two Metadynamics examples are included in the ``Examples/User/Meta`` directory. 
In the first example, Metadynamics is used to sample the free energy surface of 
a two-dimensional particle undergoing Langevin dynamics. This example is found in 
the `Single_Atom` folder and requires LAMMPS. The files included are described below:

* ``in.LAMMPS_Meta_Test`` - LAMMPS input file describing the Langevin particle 
	and underlying free energy surface to be sampled. The free energy surface consists of two
	Gaussian wells at (0.98, 0.98) and (-0.98, -0.98) respectively, and one Gaussian 
	barrier at the origin.
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

References Cited
^^^^^^^^^^^^^^^^

.. [Singh2012] Singh, Sadanand, Manan Chopra, and Juan J. de
               Pablo. "Density of states–based molecular simulations."
               Annual review of chemical and biomolecular engineering
               3 (2012): 369-394.

.. [McGovern2013] McGovern, Michael, and Juan de Pablo. "A boundary
                  correction algorithm for metadynamics in multiple
                  dimensions." The Journal of chemical physics 139.8
                  (2013): 084102.

