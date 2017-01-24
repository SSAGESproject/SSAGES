.. _tutorials:

Tutorials
=========

Basic User Tutorial
-------------------

In your SSAGES directory:

.. code-block:: bash

    cd Examples/User/Umbrella

To run a simulation using SSAGES, a ``.json`` file is needed. A ``.json`` file
will tell SSAGES what method it should run, what engine it will use, and what
parameters to use for the method chosen. It will also tell SSAGES how many
simulations to run (walkers) and what engine-input file it should read.  In this
particular example, the engine will be LAMMPS. A butane molecule is used as the
example. All the appropriate LAMMPS input files are provided. The LAMMPS input
files contains the necessary information to perform the simulation. In
``Butane_SSAEGES.in`` you will notice in the last line ``fix ssages all sages``.

A template ``.json`` file (``Template_Input.json``) is provided which contains
the necessary information for a single umbrella simulation. The python code
provided will use the template to generate a new .json file needed for the
simulation. The template ``.json`` file contains the name of the collective
variable (CV) we will use, “Torsional”, and the appropriate atom ids. Run the
python code:

.. code-block:: bash

    python Umbrella_Input_Generator.py 

A new ``.json`` file (``Umbrella.json``) will appear with the correct number of
entries. In this particular example, 12 different walkers are generated. Run
SSAGES:

.. code-block:: bash

    mpiexec -np 12 ./ssages Umbrella.json

where 12 is the number of processors. Since ``Umbrella.json`` contains 12
walkers, 12 processors should be used.

With that, SSAGES will perform Umbrella sampling on a butane molecule biasing
the torsional CV. Output files will be generated for each one of the walkers
containing the iteration number, the target value for the CV, and the CV value
at the iteration number. These values can then be used for further analysis. 

Method-specific tutorials
-------------------------

:ref:`Adaptive Biasing Force <ABF-tutorial>`

:ref:`Basis Function Sampling <BFS-tutorial>`

:ref:`Finite Temperature String <FTS_tutorial>`

:ref:`Forward Flux <FFS_tutorial>`

:ref:`Image Method <IM_tutorial>`

:ref:`Metadynamics <MD_tutorial>`

:ref:`Swarm <Swarm_tutorial>`

:ref:`Umbrella Sampling <Umbrella_tutorial>`
