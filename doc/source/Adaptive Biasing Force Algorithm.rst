.. _adaptive-biasing-force:

Adaptive Biasing Force Algorithm
--------------------------------

Introduction
^^^^^^^^^^^^

Adaptive Biasing Force is a variant of a flat histogram method. Like many
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
`here <http://pubs.acs.org/doi/abs/10.1021/jp506633n>`__.

Details on the specific implementation used in SSAGES can be found
`here <http://mc.stanford.edu/cgi-bin/images/0/06/Darve_2008.pdf>`__.

An integrator for 1D and 2D surfaces are provided in SSAGES/Tools/ABF_integrator (requires numpy, scipy and matplotlib).
ABF_integrator.py -i <inputfile> -o <outputname> --periodic1 <True/False> --periodic2 <True/False> --interpolate <integer> --scale <float>


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Adaptive Biasing Force Method

* Calculate the generalized force on CVs at each timestep
* Bias with the negative of the estimated generalized force
* Define a CV range. Outside of the CV range, there will be no bias, and no
  histogram hits will be collected.
* Can optionally define a restraint range. Outside this range, a harmonic
  restraint of user-chosen spring constant(s) will drive the CV(s) back into the
  range. This range should be WIDER than the CV range by at least one bin size
  in each direction. To disable restraints, enter a spring constant k equal to
  or less than zero. If restraints are used on a periodic system, one can define
  the periodic boundaries, so that minimum image convention to CVs can be applied.
  (CV_periodic_boundary_upper/lower_bounds). For example, on a -pi to pi CV, if the 
  CV is restrained to -3.14 to -2.36 and the CV crosses the -3.14 boundary to 3.14,
  this will ensure the restraint is applied correctly back towards -3.14 rather than
  a large force applied to bring it from 3.14 all the way to -2.36.

How to define the ABF Method: ``"type" : "ABF"``

cvs
   *array of integers*.
   This array selects which CVs this method will operate on. Index starts from 0.

CV_lower_bounds
    *array of doubles (nr of CVs) long*.
    This array defines the minimum values for the CVs for the range in which the
    method will be used in order. 

CV_upper_bounds
    *array of doubles (nr of CVs) long*.
    This array defines the minimum values for the CVs for the range in which the
    method will be used in order.

CV_bins
    *array of doubles (nr of CVs) long*.
    This array defines the number of histogram bins in each CV dimension in order.

CV_restraint_minimums
    *array of doubles (nr of CVs) long*.
    This array defines the minimum values for the CV restraints in order. 

CV_restraint_maximums
    *array of doubles (nr of CVs) long*.
    This array defines the maximum values for the CV restraints in order.

CV_restraint_spring_constants
    *array of doubles (nr of CVs) long*.
    This array defines the spring constants for the CV restraints in order.
    Enter a value equal to or less than zero to turn restraints off.

CV_isperiodic
    *array of booleans (nr of CVs) long*.
    This array defines whether a given CV is periodic for restraint purposes.
    This is only used to apply minimum image convention to CV restraints.
    Can be safely set to false even for periodic CVs if no restraints are being used.
    If ANY CV is set to periodic, then CV_periodic_boundary_lower_bounds and 
    CV_periodic_boundary_upper_bounds must be provided for ALL CVs. 
    Values entered for non-periodic CVs are not used.

CV_periodic_boundary_lower_bounds
    *array of doubles (nr of CVs) long*.
    This array defines the lower end of the period.
    This only matters if CV_isperiodic is true for the CV.

CV_periodic_boundary_upper_bounds
    *array of doubles (nr of CVs) long*.
    This array defines the upper end of the period.
    This only matters if CV_isperiodic is true for the CV.

timestep
    *double*.
    The timestep of the simulation. Units depend on the conversion factor that
    follows. This must be entered correctly, otherwise generalized force estimate
    will be incorrect.

minimum_count
    *integer*.
    Number of hits in a histogram required before the full bias is active for
    that bin. Below this value, the bias linearly decreases to equal 0 at hits = 0.
    Default = 200, but user should provide a reasonable value for their system.

mass_weighing
    *boolean*
    Turns on/off mass weighing of the adaptive force.
    Default is off. Keep off if your system has massless sites such as in TIP4P water.

filename
    *string*.
    Default = F_out
    Name of the file to save Adaptive Force Vector Field information to - this
    is what’s useful

Fworld_filename
    *string*.
    Default = Fworld_cv
    Name of the file to backup the raw Fworld output to for restarts.
    There will be separate outputs for each CV.

Nworld_filename
    *string*.
    Default = Nworld
    Name of the file to backup the raw Nworld output to for restarts.

backup_frequency
    *integer*.
    Saves the histogram of generalized force every this many timesteps.

unit_conversion
    *double*.
    Unit conversion from d(momentum)/d(time) to force for the simulation. 
    For LAMMPS using units real, this is 2390.06
    (gram.angstrom/mole.femtosecond^2 -> kcal/mole.angstrom)
    For GROMACS, this is 1.

frequency
    *1*.
    OPTIONAL
    Leave at 1.
    

Example input
^^^^^^^^^^^^^

.. code-block:: javascript

"methods" : [{
                "type" : "ABF",
		"cvs" : [0,1],
  		"CV_lower_bounds" : [-3.14, -3.14],
                "CV_upper_bounds" : [3.14,3.14],
		"CV_bins" : [21,21],
  		"CV_restraint_minimums" : [-5,-5],
                "CV_restraint_maximums" : [5,5],
		"CV_restraint_spring_constants" : [0,0],
		"CV_isperiodic" : [false,false],
		"timestep" : 0.002,
		"minimum_count" : 50,
		"filename" : "F_out",
		"backup_frequency" : 1000,
		"unit_conversion" : 1,
		"frequency" : 1
            }]

Output
^^^^^^

The main output of the method is stored in a file specified in 'filename'. This 
file will contain the Adaptive Force vector field printed out every 
'backup_frequency' steps and at the end of a simulation. The method outputs a vector 
field, with vectors defined on each point on a grid that goes from 
(CV_lower_bounds) to (CV_upper_bounds) of each CV in its dimension, with (CV_bins) of grid points 
in each dimension. For example, for 2 CVs defined from (-1,1) and (-1,0) with 3 and
2 bins respectively would be a 3x2 grid (6 grid points). The printout is in the
following format: 2*N number of columns, where N is the number of CVs. First N columns 
are coordinates in CV space, the N+1 to 2N columns are components of the Adaptive Force 
vectors. An example for N=2 is:

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

.. _ABF-tutorial:

Tutorial
^^^^^^^^

For LAMMPS (must be built with RIGID and MOLECULE packages)
To build RIGID and MOLECULE: 

1) Go to LAMMPS src folder (/build/hooks/lammps/lammps-download-prefix/src/lammps-download/src/ for -DLAMMPS=YES)
2) Do:

.. code-block:: bash

   make yes-RIGID
   make yes-MOLECULE

3) Go to your build folder and make.

Find the following input files in Examples/User/ABF/Example_AlanineDipeptide:

* ``in.ADP_ABF_Example(0-1)`` (2 files)
* ``example.input``
* ``ADP_ABF_1walker.json``
* ``ADP_ABF_2walkers.json``

1) Put the contents of ABF_ADP_LAMMPS_Example folder in your ssages build folder
2) For a single walker example, do:

.. code-block:: bash

    ./ssages ADP_ABF_1walker.json.json
    
For 2 walkers, do:

.. code-block:: bash

    mpirun -np 2 ./ssages ADP_ABF_2walkers.json

For GROMACS:

Optional:

* ``adp.gro``
* ``topol.top``
* ``nvt.mdp``

Required:

* ``example_adp(0-1).tpr`` (2 files)
* ``ADP_ABF_1walker.json``
* ``ADP_ABF_2walkers.json``

1) Put the contents of ABF_ADP_Gromacs_Example in your ssages build folder
2) For a single walker example, do:

.. code-block:: bash

    ./ssages ABF_AlaDP_1walker.json

For 2 walkers, do:

.. code-block:: bash

    mpirun -np 2 ./ssages ABF_AlaDP_2walkers.json

These will run using the pre-prepared input files in .tpr format. If you wish to
prepare the input files yourself using GROMACS tools (if compiled with -DGROMACS=YES):

.. code-block:: bash

    /build/hooks/gromacs/gromacs/bin/gmx_mpi grompp -f nvt.mdp -p topol.top -c adp.gro -o example_adp0.tpr
    /build/hooks/gromacs/gromacs/bin/gmx_mpi grompp -f nvt.mdp -p topol.top -c adp.gro -o example_adp1.tpr

Be sure to change the seed in .mdp files for random velocity generation, 
so walkers can explore different places on the free energy surface.

Multiple walkers initiated from different seeds will
explore different regions and will all contribute to the same adaptive force.

After the run is finished, you can check that your output matches the sample
outputs given in the examples folders:

1) Copy ABF_integrator.py (requires numpy, scipy and matplotlib) into your build folder.
2) Run the integrator:

.. code-block:: bash

    python ABF_integrator.py --periodic1 True --periodic2 True --interpolate 200

3) This will output a contour map, a gradient field and a heatmap. Compare these to the sample outputs.	

Developers
^^^^^^^^^^

Emre Sevgen
Hythem Sidky
