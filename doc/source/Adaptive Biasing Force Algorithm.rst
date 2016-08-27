Adaptive Biasing Force Algorithm
--------------------------------

Introduction
^^^^^^^^^^^^

Adaptive Biasing Force is, at its heart, a flat histogram method. Like many
other methods that seek uniform sampling over CV space such as Metadynamics, it
adaptively biases the simulation until such diffusive sampling is achieved.
However, unlike metadynamics, ABF does not estimate the free energy surface.
Rather, it directly estimates the derivative of the free energy in CV directions
- the generalized force on that CV by the system.

In practice, this translates to histogramming coordinates in CV space with an
instantaneous estimation of the free energy derivative. This instantaneous
estimate fluctuates around the true, global free energy derivative at that
point, but the average quickly converges to the real value. Then, the free
energy derivatives can be integrated much like Thermodynamic Integration to get
the free energy surface. 

Thus, ABF gives a vector field and not a free energy surface.

An excellent write-up on the method can be found
`here <http://pubs.acs.org/doi/abs/10.1021/jp506633n>`_.

Details on the specific implementation used in SSAGES can be found
`here <http://mc.stanford.edu/cgi-bin/images/0/06/Darve_2008.pdf>`_.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Adaptive Biasing Force Method

* Calculate the generalized force on CVs at each timestep
* Bias with the negative of the estimated generalized force
* Define a CV range. Outside of the CV range, there will be no bias, and no
  histogram hits will be collected.
* Can optionally define a restraint range. Outside this range, a harmonic
  restraint of user-chosen spring constant will drive the CV back into the
  range. This range should be WIDER than the CV range by at least one bin size
  in each direction. To disable restraints, enter a spring constant k equal to
  or less than zero.

How to define the ABF Method: ``"type" : "ABF"``

CV_lower_bounds
    *array of doubles (nr of CVs) long*.
    This array defines the minimum values for the CVs for the range in which the
    method will be used in order. 

CV_upper_bounds
    *array of doubles (nr of CVs) long*.
    This array defines the minimum values for the CVs for the range in which the
    method will be used in order.

CV_bins
    *array of doubles (cr of CVs) long*.
    This array defines the number of histogram bins in each CV dimension in order.

CV_restraint_minimums
    *array of doubles (cr of CVs) long*.
    This array defines the minimum values for the CV restraints in order. 


CV_restraint_maximums
    *array of doubles (cr of CVs) long*.
    This array defines the maximum values for the CV restraints in order.

CV_restraint_spring_constants
    *array of doubles (cr of CVs) long*.
    This array defines the spring constants for the CV restraints in order.
    Enter a value equal to or less than zero to turn restraints off.

timestep
    *double*.
    The timestep of the simulation. Units depend on the conversion factor that
    follows.

minimum_count
    *integer*.
    Number of hits in a histogram required before the full bias is active for
    that bin. Below this value, the bias linearly decreases to equal 0 at hits = 0.

filename
    *string*.
    Name of the file to save Adaptive Force Vector Field information to - this
    is what’s useful

backup_frequency
    *integer*.
    Saves the histogram of generalized force every this many timesteps.

unit_conversion
    *double*.
    Unit conversion from d(momentum)/d(time) to force for the simulation. 
    For LAMMPS using units real, this is 2390.06
    (gram.angstrom/mole.femtosecond^2 -> kcal/mole.angstrom)
    For GROMACS, this is 1.

Orthogonalization
    *1 or 0*.
    Turn on Gram-Schmidt Orthogonalization for the projector vector fields. Best
    to leave on unless CVs are orthogonal (share no gradient directions)

frequency
    *1*.
    Leave at 1.

Example input
^^^^^^^^^^^^^

.. code-block:: javascript

    "method" : {
            "type" : "ABF",                
            "CV_lower_bounds" : [-3.13, -3.13],
            "CV_upper_bounds" : [3.13,3.13],
            "CV_bins" : [91,91],
            "CV_restraint_minimums" : [-5,-5],
            "CV_restraint_maximums" : [5,5],
            "CV_restraint_spring_constants" : [0,0],
            "timestep" : 0.002,
            "minimum_count" : 200,
            "filename" : "F_out",
            "backup_frequency" : 1000000,
            "unit_conversion" : 1,
            "orthogonalization" : 0,
            "frequency" : 1
    }

Output
^^^^^^

The method outputs a vector field, with vectors defined on each point on a grid
that goes from (min_CV) to (max CV) of each CV in its dimension, with (nr of
bins) of grid points in each dimension. For example, for 2 CVs defined from
(-1,1) and (-1,0) with 3 and 2 bins respectively would be a 3x2 grid (6 grid
points). There are several outputs of the method. For most people, the file
specified in “file name” should suffice. This file will contain the Adaptive
Force vector field printed out every 'Backup interval' steps and at the end of a
simulation in the following format: 2*N number of columns, where N is the number
of CVs. First N columns are coordinates in CV space, the N+1 to 2N columns are
components of the Adaptive Force vectors. An example for N=2 is:

+-----------+-----------+-------------+-------------+
| CV1 Coord | CV2 Coord | d(A)/d(CV1) | d(A)/d(CV2) |
+===========+===========+=============+=============+
| -1        | -1        | -1          | 1           |
+-----------+-----------+-------------+-------------+
| -1        | 0         | 2           | 1           |
+-----------+-----------+-------------+-------------+
| 0         | -1        | 1           | 2           |
+-----------+-----------+-------------+-------------+
| 0         | 0         | 2           | 3           |
+-----------+-----------+-------------+-------------+
| 1         | -1        | 2           | 4           |
+-----------+-----------+-------------+-------------+
| 1         | 0         | 3           | 5           |
+-----------+-----------+-------------+-------------+

Tutorial
^^^^^^^^

Find the following input files in Examples/User/ABF:

For LAMMPS:

* ``in.ADP_ABF_Example``
* ``example.input``
* ``ABF_AlaDP_1walker.json``
* ``ABF_AlaDP_8walkers.json``

1) Put all files in your ssages build folder
2) For a single walker example, do:

.. code-block:: bash

    mpirun -np 1 ./ssages -ABF_AlaDP_1walker.json
    
For 8 walkers, do:

.. code-block:: bash

    mpirun -np 8 ./ssages -ABF_AlaDP_8walkers.json

Multiple walkers initiated from different seeds using ssages’ SEED keyword will
explore different regions and will all contribute to the same adaptive force.

3) After the run is finished open F_out and copy the last grid that defined the
   Adaptive Force vector field (all numbers in four columns after the last line
   of text)
4) Paste into any new folder, run ABF_analysis.py (requires numpy, scipy and
   matplotlib)

For GROMACS:

Optional:

* ``adp.gro``
* ``topol.top``
* ``nvt1-8.mdp`` (8 files)

Required:

* ``example1-8.tpr`` (8 files)
* ``ABF_AlaDP_1walker.json``
* ``ABF_AlaDP_8walkers.json``

1) Put all files in your ssages build folder
2) For a single walker example, do:

.. code-block:: bash

    mpirun -np 1 ./ssages -ABF_AlaDP_1walker.json

For 8 walkers, do:

.. code-block:: bash

    mpirun -np 8 ./ssages -ABF_AlaDP_8walkers.json

These will run using the pre-prepared input files in .tpr format. If you wish to
prepare input files yourself using GROMACS tools:

.. code-block:: bash

    gmx grompp -f nvt1.mdp -p topol.top -c adp.gro -o example1.tpr
    gmx grompp -f nvt2.mdp -p topol.top -c adp.gro -o example2.tpr
    .
    .
    gmx grompp -f nvt8.mdp -p topol.top -c adp.gro -o example8.tpr

The only difference between .mdp files is the seed for random velocity
generation, so walkers can explore different places on the free energy
surface.

Developer
^^^^^^^^^

Emre Sevgen
