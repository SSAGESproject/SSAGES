.. _basis-function-sampling:

Basis Function Sampling
-----------------------

Introduction
^^^^^^^^^^^^

The Basis Function enhanced sampling method is a variant of the Continuous
Wang-Landau Sampling method developed by Whitmer et al, which biases a PMF
through the summation of Kronecker deltas. In this method, the Kronecker delta
is approximated by projection of a locally biased histogram to a truncated set
of orthogonal basis functions.

.. math::

    \int_\Xi f_{i}(\vec{\xi})f_{j}(\vec{\xi})w(\vec{\xi})d\vec{\xi} = \delta_{i}c_{i}

By projecting a basis set, the system resolves the same properties as the
Kronecker deltas, but in a continuous and differentiable manner that lends well
towards MD simulations. Currently in SSAGES, Legendre Polynomials have been
implemented to work with BFS. These have the property where the weight
:math:`w(\xi) = 1` and are defined on the interval :math:`[-1, 1]`.

The method applies its bias in sweeps of $N$ through a histogram (:math:`H_{i}`)
that is updated at every :math:`j` microstate or timestep. This histogram is
then modified to an unbiased partition function estimate (:math:`\tilde{H_{i}}`)
by convolution with the current bias potential (:math:`\Phi_{i}`).

.. math::

    \tilde{H}_{i}(\xi) = H_{i}(\xi)e^{\beta \Phi_{i}}

In order to account for sampling history into this partition function estimate,
a simple weight function (:math:`W(t_{j})`) is added. 

.. math::

    Z_{i}(\xi) = \sum_{j} W(t_{j})\tilde{H_{j}}(\xi)

This final estimate is then projected to the truncated basis set. After this set
is evaluated, the coefficients of the basis set are evaluated. This process is
iterated until the surface converges, which is determined by the overall update
of the coefficients.

.. math::

    \beta \Phi_{i+1}(\xi) = \sum_j^N \alpha^i_j L_j(\xi)\\
    \alpha^i_j = \frac{2j + 1}{2} \int_{-1}^1 \log(Z_i(\xi))L_j(\xi)d\xi

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

These are all the options that SSAGES provides for running Basis Function
Sampling. In order to add BFS to the JSON file, the method should be labeled as
"Basis".

CV_coefficients
    The order of the polynomial to be projected for each collective variable. If
    the order of this array doesn't match the number of CVs, the system assumes
    the first number for all of the CVs

CV_restraint_spring_constants
    The strength of the springs keeping the system in bounds in a non-periodic
    system.

CV_restraint_maximums
    The upper bounds of each CV in a non-periodic system.

CV_restraint_minimums
    The lower bounds of each CV in a non-periodic system.

cycle_frequency
    The frequency of updating the projection bias.

frequency
    The frequency of each integration step. This should almost always be set to 1.

weight
    The weight of each visited histogram step. Should be kept at 1.0, but the
    option is available to make it slightly greater. The system has a higher
    chance of exploding at higher weight values.

basis_filename
    A suffix to name the output file. If not specified the output will be
    "basis.out"

coeff_filename
    A suffix to name the coefficient file.

temperature
    Only should be used if the MD engine cannot produce a good temperature
    value. (ie: LAMMPS with 1 particle)

tolerance
    Convergence criteria. The sum of the difference in subsequent updates of the
    coefficients squared must be less than this for convergence to work.

convergence_exit
    A boolean option to let the user choose if the system should exit once the
    convergence is met.

Required to Run BFS
^^^^^^^^^^^^^^^^^^^

In order to use the method properly a few things must be put in the JSON file. A
grid is required to run Basis Function Sampling. Refer to the Grid section in
order to understand options available for the grid implementation.
The only inputs required to run the method:

* cyclefrequency
* frequency
* CV coefficients

Example Input
^^^^^^^^^^^^^
.. code-block:: javascript

    "methods" : [{
        "type" : "Basis",
        "cvs" : [0],
        "cycle_frequency" : 10000,
        "frequency" : 1,
        "weight" : 1.0,
        "basis_filename" : "nacl_example",
        "coeff_filename" : "nacl_example",
        "CV_coefficients" : [ 50 ],
        "CV_restraint_spring_constants" : [ 30 ],
        "CV_restraint_maximums" : [ 10.3 ],
        "CV_restraint_minimums" : [ 1.7 ],
        "grid" : {
            "lower" : [2.0],
            "upper" : [15.0],
            "number_points" : [500],
            "periodic" : [false]
        }
    }]

Guidelines for running BFS
^^^^^^^^^^^^^^^^^^^^^^^^^^

* It is generally a good idea to use polynomials of order at least 25. 
* For higher order polynomials, the error in projection is less, but the number
  of bins must increase in order to accurately project the surface.
* A good rule of thumb for these simulations is to do at least one order of
  magnitude more bins than polynomial order.

If the system that is to be used requires a non-periodic boundary condition,
then it is typically a good idea to place the bounds approximately 0.1 - 0.2
units outside the grid boundaries.

The convergence exit option is available if the user chooses to continue running
past convergence, but a good heuristic for tolerance is around
:math:`1\mathrm{e}{-6}`.

.. _BFS-tutorial:

Tutorial
^^^^^^^^

This tutorial will provide a reference for running BFS in SSAGES. There are
multiple examples provided in the Examples/User directory of SSAGES, but this
tutorial will cover the Alanine Dipeptide example. 
In the ADP subdirectory of the ``Examples/User section`` there should be a
LAMMPS input file (titled ``in.ADP_BFS_example(0-1)``) and two JSON input files.
Both of these files will work for SSAGES, but the one titled ``ADP_BFS_2walkers.json``
makes use of multiple walkers.

For LAMMPS to run the example it must be made with RIGID and MOLECULE options.
In order to do so, 

1) Go to LAMMPS src folder (/build/hooks/lammps/lammps-download-prefix/src/lammps-download/src/ for -DLAMMPS=YES)
2) Do:

.. code-block:: bash

   make yes-RIGID
   make yes-MOLECULE

3) Go to your build folder and make.

Use the following command to run the example:

.. code-block:: bash

    mpiexec -np 1 ./ssages ADP_BFS_2walkers.json

This should prompt SSAGES to begin an alanine dipeptide run. If the run is
successful, the console will output the current sweep number on each node.
At this point the user can elect to read the output information after each sweep. 

basis.out
~~~~~~~~~

The ``basis.out`` file outputs in at least 4 columns. These columns refer to the
CV values, the ultimate projected PMF, the unprojected PMF, and the biased
histogram values. Depending on the number of CVs chosen for a simulation, the
number of CV columns will also correspond. Only the first CV column should be
labeled.

The important line for graphing purposes is the projected PMF, which is the
basis set projection from taking the log of the biased histogram. The biased
histgram is printed so that it can be read in for doing restart runs (subject to
change). For plotting the PMF, a simple plotting tool over the CV value and
projected PMF columns will result in the free energy surface of the simulation.
The free energy surface will return a crude estimate within the first few
sweeps, and then will take a longer period of time to retrieve the fully
converged surface. A reference image of the converged  alanine dipeptide example
is provided in the same directory as the LAMMPS and JSON input files.

coeff.out
~~~~~~~~~

This holds all the coefficient values after each bias projection update. This
file is entirely used for restart runs.

Developer
^^^^^^^^^

Joshua Moller
Julian Helfferich

