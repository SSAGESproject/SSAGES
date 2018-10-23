.. _elastic-band:

Elastic Band
------------

Introduction
^^^^^^^^^^^^

There are many methods, several of which are included in SSAGES, to calculate
transition pathways between metastable states. One kind of pathway between
states in the *Minimum Energy Pathway* (MEP), quite simply the lowest energy
pathway a system can take between these states. An MEP has the condition that
the force everywhere along the pathway points only along the path, that is, it
has no perpendicular component. By finding the MEP, one also finds the saddle
points of the potential energy surface, as they are by definition the maxima of
the MEP. The *Nudged Elastic Band* (NEB) method is a popular and efficient
method to calculate the MEP between the initial and final state of a transition
:cite:`HENKELMAN20009978,HENKELMAN20009901`.

The method involves the evolution of a series of images connected by a spring
interaction (hence the "elastic" nature of the band). The force acting on the
images (a combination of the spring force along the band and the true force
acting perpendicular to the band) is minimized to ensure convergence to the MEP.
The **nudged** nature of NEB refers to a force projection that ensures the
spring forces do not interfere with the elastic band converging to the MEP, as
well as that the true force does not alter the distribution of images along the
band (that is, it ensures all the images do not fall into the metastable states).
This projection is accomplished by using the parallel portion of the spring
force and the perpendicular portion of the true force. In this way, the spring
forces act similarly to reparameterization schemes common to the string method.

Full mathematical background is available in :cite:`HENKELMAN20009978` and
:cite:`HENKELMAN20009901`, but a brief overview is given here. The
band is discretized as a series of N+1 images, and the force on each image is
given by:

.. math::

	F_{i} = F_{i,\parallel}^{s} - \nabla E(R_{i})_{\perp}

Where :math:`F_{i}` is the total force on the image, :math:`F_{i,\parallel}^{s}`
refers to the parallel component of the spring force on the
:math:`i^\text{th}` image, and
:math:`\nabla E(R_{i})_{\perp}` is the perpendicular component of the gradient
of the energy evaluated at each image :math:`R_{i}`. The second term on the
right hand side is the "true force" and is evaluated as:

.. math::

	\nabla E(R_{i})_{\perp} = \nabla E(R_{i}) - \nabla E(R_{i})\cdot\hat{\tau_{i}}

The term :math:`\hat{\tau_{i}}` represents the normalized local tangent at the
:math:`i^\text{th}` image, and thus this equation states simply that the perpendicular component
of the gradient is the full gradient minus the parallel portion of the gradient.
There are different schemes available in literature to evaluate the tangent
vector :cite:`HENKELMAN20009978`. The "spring force" is calculated as:

.. math::

	F_{i,\parallel}^{s} =
	k \left( \lvert R_{i+1} - R_{i} \rvert -
	         \lvert R_{i} - R_{i-1} \rvert \right) \cdot \hat{\tau_{i}}

Where :math:`k` is the spring constant, which can be different for each image of
the band. One can evolve the images with these forces according to any number
of schemes - a straightforward Verlet integration scheme is used in the SSAGES
implementation, described below.

Algorithmically, the NEB method is implemented in SSAGES as follows:

1. An initial band is defined between the two states of interest. This can be
   defined however one wishes; often it is simply a linear interpolation through
   the space of the collective variables. In fact, the ends of the band need
   not necessarily be in the basins of interest; the method should allow the
   ends to naturally fall into nearby metastable basins.

2. For each image of the band, a molecular system with atomic coordinates that
   roughly correspond to the collective variables of that image is constructed.
   A period of equilibration is performed to ensure that the underlying systems'
   CVs match their respective band images.

3. The gradient is sampled over a user-defined period of time and intervals,
   this being the only quantity with statistical variance that needs to be
   averaged over.

4. When sufficient sampling of the gradient is done, the band is updated one
   time-step forward with a simple Verlet scheme.

Steps 2--4 are iterated upon, leading to convergence of the method and the MEP.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

To construct an EB input file, the following options are available. A
complete EB JSON file will inherit some of its inputs from the String
schema (for parameters common to all string methods).
The options unique to EB are:

equilibration_steps
	The number of MD steps to simply perform umbrella sampling without
	invoking the NEB method. A sufficiently long number of steps ensures
	that the underlying molecular systems have CVs close to the CVs of their
	associated image on the band.

evolution_steps
	The number of steps to perform the NEB over; the band is updated after
	evolution steps times the number of samples total MD steps. A new value
	of the gradient is harvested every time the number of MD steps taken is
	an integer multiple of evolution steps.

kstring
	The constant used in calculating the spring force at each image. It
	can be specified uniquely for each image. Please notice its difference
	from ksprings.

From the String schema, the options are:

type
	This parameter identifies that a String-type method is being used, and
	thus should be set to ``"String"``.

flavor
	This parameter identifies the specific kind of string-type method
	being used; for EB, it should be set to ``"ElasticBand"``.

centers
	This parameter assigns a location in CV space for each individual image
	along the string/elastic band. This should be an array with size equal to
	the total number of images, with each entry consisting of an array with size
	equal to the number of CVs used for the elastic band method.

tolerance
	This is a tolerance threshold that can be set to trigger the end of
	the method; it is a percentage by which, if no node CV changes by this
	percentage, the method will end. It must be specified as an array with
	one entry for each CV desired.

max_iterations
	A complementary stopping criterion can be specified; the method will
	stop if it undergoes this many iterations of the string method.

ksprings
	A unique spring constant must be defined for each CV; its purpose is
	described above.

frequency
	The frequency of each integration step. This should almost always be
	set to 1.

.. _EB_tutorial:

Tutorial
^^^^^^^^

This tutorial will walk you step by step through the user example provided with
the SSAGES source code that runs the NEB method on the alanine dipeptide using
LAMMPS. First, be sure you have compiled SSAGES with LAMMPS. Then, navigate to
the ``Examples/User/ElasticBand/ADP`` subdirectory. Now, take a moment
to observe the ``in.ADP_Test and data.input`` files in order to familiarize
yourself with the system being simulated.

The next two files of interest are the ``Template_Input.json`` input file and the
``Input_Generator.py`` script. Both of these files can be modified in your
text editor of choice to customize the inputs, but for this tutorial, simply
observe them and leave them be. ``Template_Input.json`` contains all the information
necessary to fully specify one driver; ``Input_Generator.py`` copies this
information a number of times specified within the script (for this tutorial,
22 times) while also linearly interpolating through the start and end states
defined in the script and substituting the correct values into the "centers"
portion of the method definition. Execute this script as follows:

.. code-block:: bash

	python Input_Generator.py

You will produce a file called ``ElasticBand.json``. You can also open this file to
verify for yourself that the script did what it was supposed to do. Now, with
your JSON input and your SSAGES binary, you have everything you need to perform
a simulation. Simply run:

.. code-block:: bash

	mpiexec -np 22 ./ssages ElasticBand.json

Soon, the simulation will produce a ``node-X.log`` file for each driver, where
X is the number specifying the driver (in this case, 0-21 for our 22 drivers).
Each one will report the following information, in order: the node number, the
iteration number, and for each CV, the current value of the band CV as well as
the current value of the CV calculated from the molecular system.

Allow your system to run for the specified number of iterations (1000 for this
tutorial). The last line of every node file can be analyzed to view the last
positions of each image of the elastic band.

Developer
^^^^^^^^^

* Benjamin Sikora
