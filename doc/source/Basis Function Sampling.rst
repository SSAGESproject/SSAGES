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
towards MD simulations. The current version of SSAGES has support for Legendre, 
Chebyshev, and Fourier polynomials. Each of these has their defined weight function :math:`w(\xi)`
implemented specific to the method. Additionally, any combination of implemented basis sets can be
used for any system. It is advised that a periodic basis set be used with a periodic CV, but it
is not required.

The BFS method applies its bias in sweeps of $N$ through a histogram (:math:`H_{i}`)
that is updated at every :math:`j` microstate or timestep. This histogram is
then modified to an unbiased partition function estimate (:math:`\tilde{H_{i}}`)
by exponentiation with the current bias potential (:math:`\Phi_{i}`).

.. math::

    \tilde{H}_{i}(\xi) = H_{i}(\xi)e^{\beta \Phi_{i}}

A weight function has been added into this implementation (:math:`W(t_{j})`) so that
the user can define the effective strength of the applied bias. If not chosen, the weight is normalized to
the length of the interval.

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
"BFSMethod".

Basis Function Sampling requires the use of a basis set. These are defined by defining
an object of "basis_functions". These have the following properties

type
    Currently can either be Chebyshev, Fourier, or Legendre

polynomial_order
    Order of the polynomial. In the case of Chebyshev or Legendre this results in an
    order of input value + 1 as the method takes the 0th order internally. For a Fourier
    series, the order is the total number of coefficients including the sine and cosine series.

upper_bound
    Only exists for Chebyshev and Fourier series. This is the upper bound of the CV

lower_bound
    Only exists for Chebyshev and Fourier series. This is the lower bound of the CV

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
    The weight of each visited histogram step. Should be kept around the same value
    as the cycle_frequency (usually 0.1 times that) The system has a higher
    chance of exploding at higher weight values.

basis_filename
    A suffix to name the output file. If not specified the output will be
    "basis.out"

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
* basis_functions

Example Input
^^^^^^^^^^^^^
.. code-block:: javascript

  "methods" : [{
      "type" : "BFSMethod",
      "basis_functions" : [
      {
          "type" :"Fourier", 
          "polynomial_order" : 30,
          "upper_bound" : 3.14,
          "lower_bound" : -3.14
      },
      {
          "type" : "Fourier",
          "polynomial_order": 30,
          "upper_bound" : 3.14,
          "lower_bound" : -3.14
      }],
      "cvs" : [0,1],
      "cycle_frequency" : 100000,
      "basis_filename" : "example",
      "frequency" : 1,
      "temperature" : 1.0,
      "weight" : 100000.0,
      "tolerance" : 1e-6,
      "convergence_exit" : false,
      "grid" : {
          "lower" : [-3.14, -3.14],
          "upper" : [3.14,3.14],
          "number_points" : [100,100],
          "periodic" : [true, true]
      }
  }]

Guidelines for running BFS
^^^^^^^^^^^^^^^^^^^^^^^^^^

* It is generally a good idea to choose a lower order polynomial initially.
  Excessive number of polynomials may create an unwanted "ringing" effect that could result
  in much slower convergence.
* For higher order polynomials, the error in projection is less, but the number
  of bins must increase in order to accurately project the surface. This may
  also result in an undesired ringing phenomena.
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

    mpiexec -np 2 ./ssages ADP_BFS_2walkers.json

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

restart.out
~~~~~~~~~

This holds all the coefficient values after each bias projection update, as well
as the biased histogram. This file is entirely used for restart runs.

Developer
^^^^^^^^^

Joshua Moller
Julian Helfferich

