.. swarm:

Swarm of Trajectories
---------------------

Introduction
^^^^^^^^^^^^

Like all string methods in general, the **string method with swarms of
trajectories** (often abbreviated to "swarm of trajectories" or even more simply
"SoT") is a method to identify a transition pathway in an arbitrarily
high-dimensional collective variable space between two metastable states of a
system.  This pathway (the string) is a parametrized curve discretized into a
set of images, each of which is itself a molecular system.  The classical *string
method in collective variables* evolves each image by estimating a mean force and
metric tensor at each image with restrained molecular dynamics simulations.  In
the SoT method, the string is instead evolved by launching a large number (a
swarm) of **unrestrained** trajectories from each image and estimating the average
drift of the collective variables over the swarm.  

The mathematical background of the method can be expressed in a few relatively
straightforward equations, with further detail available in the original work of
Benoit Roux and collaborators [1]_.  First, consider a path :math:`z(\alpha)`
constructed between two metastable states, such that :math:`\alpha=0` represents
the starting state and :math:`\alpha=1` is the final state.  The "most probable
transition pathway" (MPTP) is defined such that a molecular system started from
anywhere on the path will most probably evolve while staying on the path.  It is
shown in the original work that a mathematical definition for such a path is
given when the collective variables evolve according to:

.. math::

   z_{i}(\alpha) = z_{i}(\alpha') + \sum\limits_{j}\left(
   \beta D_{ij}\left[ z(0) \right] F_{j}\left[z(0)\right] +
   \frac{\partial}{\partial z_{j}}\left( D_{ij}\left[z(0)\right]\right)
   \right)\delta\tau

Where the following notation is used: :math:`z_{i}` represents the
collective variables belonging to the string, :math:`\alpha` represents
the parameter identifying that point on the string, :math:`\beta`
represents the temperature, :math:`D_{ij}` represents the diffusion
tensor, :math:`F_{j}` represents the mean force, :math:`z` represents
the collective variables constructed from the molecular system at a
given moment in time, and :math:`\delta\tau` represents the time step of
the evolution of the dynamics. The SoT method approximates this equation
using the average drift evaluated from a large number of unbiased
trajectories, each of length :math:`\delta\tau`, launched from each
image:

.. math::

   \overline{\Delta z_{i}(\delta\tau)} = \overline{z_{i}(\delta\tau) - z_{i}(0)} \equiv
   \sum\limits_{j} \left( \beta D_{ij}\left[z(0)\right] F_{j}\left[z(0)]\right] +
   \frac{\partial}{\partial z_{j}}\left( D_{ij}\left[ z(0)\right]\right)\right)\delta\tau

Like all string methods, there is an additional step beyond evolving the
collective variables - after one iteration of evolution, the images
along the path must be reparametrized such that they lie (for example)
an equal arc length apart. This step is necessary to ensure that all
images do not fall into one metastable basin or the other.

Algorithmically, the SoT method is implemented as follows:

#. An initial string is defined between the two states of interest. This
   can be defined however one wishes; often it is simply a linear
   interpolation through the space of the collective variables. In fact,
   the ends of the string need not necessarily be in the basins of
   interest; the dynamic nature of the method should allow the ends to
   naturally fall into nearby metastable basins.

#. For each image of the string, a molecular system with atomic
   coordinates that roughly correspond to the collective variables of
   that image is constructed.

#. A set of equilibrium trajectories are generated from that system by
   performing restrained sampling around the image’s collective
   variables.

#. That set of equilibrium trajectories is used as the starting point of
   a large number of short unbiased trajectories; the resulting average
   displacement of each collective variable is used to update the
   positions of the images.

#. A reparameterization scheme is enforced to ensure that, for example,
   the string images are equally distant in collective variable space.

Steps two through five are iterated upon, leading to convergence of the
method and the MPTP.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

| To construct a Swarm input file, the following options ae available. A
complete Swarm JSON file will inherit some of its inputs from the String
schema (for parameters common to all string methods) as well as the
Observer schema (for restarts). The options unique to Swarm are:
| ***initial\_steps***
| For each iteration of the method, this is the number of steps to spend
doing restrained sampling and not harvesting trajectories. This time is
important to ensure the underlying molecular system’s CV values are
close to the string CV values.
| ***harvest\_length***
| After the initial restraining is finished, a trajectory is harvested
for later use in launching an unrestrained trajectory every so often -
harvest length specifies how often this will be done. Harvest length
multiplied by number of trajectories (see below) will determine overall
how many more steps will be taken under restrained sampling.
| ***number\_of\_trajectories***
| The total number of unrestrained trajectories to be included in each
swarm.
| ***swarm\_length***
| The length of each unrestrained trajectory in the swarm. Swarm length
multiplied by number of trajectories specifies how many total steps will
be spent doing unrestrained sampling.

| From the String schema, the options are:

| ***type***
| This parameter identifies that a String-type method is being used, and
thus should be set to “String”
| ***flavor***
| This parameter identifies the specific kind of string-type method
being used; for swarm, it should be set to “Swarm”.
| ***centers***
| For each driver, the initial values of each CV should be specified as
a list under “centers”. In this way, the initial string is defined.
| ***tolerance***
| This is a tolerance threshold that can be set to trigger the end of
the method; it is a percentage by which, if no node CV changes by this
percentage, the method will end. It must be specified as an array with
one entry for each CV desired.
| ***max\_iterations***
| A complementary stopping criterion can be specified; the method will
stop if it undergoes this many iterations of the string method.
| ***ksprings***
| A unique spring constant must be defined for each CV; its purpose is
described above.
| ***frequency***
| The frequency of each integration step. This should almost always be
set to 1.

| From the Observer schema, the options are:

| ***type***
| This can currently only be set to “JSON”; specifies to write out JSON
type restart files.
| ***file name***
| This is a prefix which will be attached to the various restart files
that will be written out.
| ***frequency***
| This specifies how often a restart file is written out, in terms of MD
steps taken.

Tutorial
^^^^^^^^

| This tutorial will walk you step by step through the user example
provided with the SSAGES source code that runs the SoT method on the
alanine dipeptide using LAMMPS. First, be sure you have compiled SSAGES
with LAMMPS. Then, navigate to the SSAGES/Examples/User/Swarm/ADP
subdirectory. Now, take a moment to observe the in.ADP\_Test and
data.input files. In general, these should be the same as what you would
use for any other method, but for the SoT method, it is important to
define a larger skin distance than one normally would in the neighbor
command in LAMMPS. This is because, under the hood, each unrestrained
trajectory in the swarm is started by manually resetting the positions
of each atom in the LAMMPS simulation to the start of a new trajectory.
From the perspective of LAMMPS, this is a huge amount of distance to
move in a single time step; this move triggers neighbor list rebuilding,
but LAMMPS considers it a “dangerous build” which threatens to crash the
simulation. Thus, we increase the skin distance, which forces LAMMPS to
keep track of more pairs in the neighbor lists, and thus reduces the
number of dangerous builds. Keep this in mind for future runs of the SoT
method.

| The next two files of interest are the Template\_Input.json input file
and the Input\_Generator.py script. Both of these files can be modified
in your text editor of choice to customize the inputs, but for this
tutorial, simply observe them and leave them be. Template\_Input.json
contains all the information necessary to fully specify one driver;
Input\_Generator.py copies this information a number of times specified
within the script (for this tutorial, 12 times) while also linearly
interpolating through the start and end states defined in the script and
substituting the correct values into the “centers” portion of the method
definition. Execute this script as follows:

| python Input\_Generator.py

| You will produce a file called Swarm.json. You can also open this file
to verify for yourself that the script did what it was supposed to do.
Now, with your JSON input and your SSAGES binary, you have everything
you need to perform a simulation. Simply run:

| mpiexec -np 12 ./ssages Swarm.json

| Soon, the simulation will produce a node-X.log file for each driver,
where X is the number specifying the driver (in this case, 0-11 for our
12 drivers). Each one will report the following information, in order:
the node number, the iteration number, and for each CV, the current
value of the string CV as well as the current value of the CV calculated
from the molecular system.

| Allow your system to run for the desired number of MD steps, but keep
an eye on it - the system should exit once one driver reaches the
maximum number of MD steps, but it is possible that instead one driver
will exit and the rest will get stuck. Check in on your node files and
see if they have been updated recently - if not, the simulation has
likely finished. Once this is done, you can execute the included
plotter.py function in a directory containing the node files with the
command line argument of how many images your string had. The script
also accepts an argument to plot a free energy surface alongside the
string, but that goes beyond the scope of this tutorial. Thus, simply
execute:

| python plotter.py 12 none

And in a moment you should have a graph of your converged string. Thus
concludes this tutorial.

Developer
^^^^^^^^^

Cody Bezik.

References
^^^^^^^^^^

.. [1] Pan, A. C., Sezer, D. & Roux, B. *Finding Transition Pathways Using the
       String Method with Swarms of Trajectories*.
       J. Phys. Chem. B **112**, 3432–3440 (2008).
