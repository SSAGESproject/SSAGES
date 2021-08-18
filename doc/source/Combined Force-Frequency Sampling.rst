.. _combined-force-frequency-sampling:

Combined Force-Frequency Sampling
---------------------------------

Introduction
^^^^^^^^^^^^

Combined Force-Frequency (CFF) sampling is a free energy sampling method which uses
artificial neural networks to generate an on-the-fly adaptive bias capable of rapidly
resolving free energy landscapes. It is a recent method proposed in :cite:`SEVGEN2020`,
and its extension and application to ab initio MD has been demonstrated in
:cite:`LEE2020`. The method's main strength resides in its ability to learn both from the
frequency of visits to distinct states and the generalized force estimates that arise in a
system as it evolves in phase space. This is accomplished by introducing a
self-integrating artificial neural network, which generates an estimate of the free energy
directly from its derivatives.

CFF algorithm proceeds in sweeps. Statistics are collected over a user specified interval
(``sweep``), which is a very flexible choice. At the end of each sweep artificial neural
networks are fit to the statistics to determine the optimal bias which is then applied in
the subsequent sweep. This proceeds until the free energy landscape has converged.

Example Input
^^^^^^^^^^^^^

.. code-block:: javascript

    "methods" : [
        {
            "type" : "CFF",
            "topology" : [12,8],
            "nsweep" : 10000,
            "temperature" : 298.15,
            "grid" : {
                "lower" : [-3.14159,-3.14159],
                "upper" : [3.14159,3.14159],
                "number_points" : [30,30],
                "periodic" : [false,false]
            },
            "lower_bounds" : [-5,-5],
            "upper_bounds" : [5,5],
            "lower_bound_restraints" : [0.0,0.0],
            "upper_bound_restraints" : [0.0,0.0],
            "timestep" : 0.001,
            "minimum_count" : 3000,
            "overwrite_output": true
        }
    ]

.. warning::

    Be sure to follow correct JSON syntax for your input, with a comma after every line
    except the last within each bracket.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

**Define CFF**

In the methods block, define the CFF method through the syntax:

.. code-block:: javascript

    "type" : "CFF",

To define neural network topology:

.. code-block:: javascript

    "topology" : [12,8],

Array of integers (length: number of hidden layers). This array defines the architecture
of the neural network. Like ANN Sampling, good heuristics is to use a network with a
single hidden layer if you are biasing on one CV, and two hidden layers if you are biasing
on two or more.

.. code-block:: javascript

    "temperature" : 298.15,

The temperature of the simulation must be specified.

To define the length of each sweep:

.. code-block:: javascript

    "nsweep" : 10000,

Typical values range from 1,000 to 10,000 depending on the size of the system. The slower
the system dynamics, the longer the sweep. This is not going to heavily affect convergence
time, and the method is generally quite robust to choice of sweep length. The main
consequence of this choice is that the neural network training time becomes relatively
expensive if the system's free energy is very cheap to evaluate.


**Define the grid**

To define the bounds:

.. code-block:: javascript

    "lower" : [-3.14159,-3.14159],
    "upper" : [3.14159,3.14159],

These are arrays of doubles whose length is the number of CVs used. This defines the
minimum and maximum values for the CVs for the range in which the method will be used in
order.

To define the number of CV bins used:

.. code-block:: javascript

    "number_points" : [30,30],

This array of integers defines the number of histogram bins in each CV dimension in order.

.. code-block:: javascript

    "periodic" : [false,false],

This array defines whether a given CV is periodic for restraint purposes. This is only
used to apply minimum image convention to CV restraints. The value can be safely set to
``false`` *even for periodic CVs* if no restraints are being used.

**Define the restraints**

.. code-block:: javascript

    "lower_bounds" : [-5,-5],
    "upper_bounds" : [5,5],

These arrays define the minimum and maximum values for the CV restraints in order.

.. code-block:: javascript

    "lower_bound_restraints" : [0,0],
    "upper_bound_restraints" : [0,0],

These arrays define the spring constant for the lower and upper bounds.

**Define time parameter**

.. code-block:: javascript

    "timestep" : 0.001,

The timestep of the simulation. Units depend on the conversion factor that follows. This
must be entered correctly, otherwise the generalized force estimate will be incorrect.

.. code-block:: javascript

    "minimum_count" : 3000,

This is the number of hits required to a bin in the general histogram before the full
biasing force is active. Below this value, the bias linearly decreases to zero at
``hits = 0``. Default value is ``200``, but user should provide a reasonable value for
their system.

**Handle outputs (Optional)**

.. code-block:: javascript

    "overwrite_output" : [true],

If this is enabled, output files are overwritten at each sweep such.  Otherwise, output
files are saved at each sweep. Default, ``true``.

Output
^^^^^^

There are six output files from this method: ``CFF.out``, ``F_out``, ``netstate.net``,
``netstate.net2``, ``CFF.out_gamma``, and ``traintime.out``.

The main output of this method is stored in ``CFF.out``. Each column corresponds
respectively to the CVs, visit frequencies (histogram), bias based on frequency-based ANN,
bias based on force-based ANN, average bias, and the average free energy estimate (which
is the negative of the average bias with a constant shift). The format is as follows:

``cv1 cv2 ... hist bias(freq_only) bias(force_only) bias(avg) free_energy(avg)``

A file called ``CFF.out_gamma`` outputs network complexity term ``gamma`` for each neural
net and the ratio of gammas from both neural nets. (See :cite:`SEVGEN2020` for more
information.) The format is as follows:

``sweep_iter gamma(freq_only) gamma(force_only) gamma_ratio``

File ``traintime.out`` contains CPU wall time (in seconds) taken for training of neural
networks during each sweep.

Files called ``netstate.dat`` and ``netstate2.dat`` contain the neural network parameters
for each of the two neural networks used to apply biases during the sampling
(``netstate.dat`` stores frequency-based parameters, while ``netstate2.dat`` stores
force-based ones). For more information about these files, see the ANN Sampling Method.

A file called ``F_out`` contains the Generalized Force vector field, which is the same
output file as the Adaptive Biasing Force method. Vectors defined on each point on a grid
that goes from ``CV_lower_bounds`` to ``CV_upper_bounds`` of each CV in its dimension, with
``CV_bins`` of grid points in each dimension. The printout is in the following format: 2*N
number of columns, where N is the number of CVs. First N columns are coordinates in CV
space, the N+1 to 2N columns are components of the Generalized Force vectors. See the ABF
Sampling Method for more information.

An example of CFF (e.g., alanine dipeptide in water) is located in ``Examples/User/CFF/ADP``.

Developers
^^^^^^^^^^

* Elizabeth M.Y. Lee
* Emre Sevgen
* Boyuan Yu


.. warning::

    Please make sure to cite :cite:`SEVGEN2020` and :cite:`LEE2020` if you use this
    method!
