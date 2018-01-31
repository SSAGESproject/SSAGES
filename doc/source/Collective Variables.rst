.. _cvs:

Collective Variables
====================

Collective variables (CVs) are arbitrary differentiable functions of the 3N Cartesian 
coordinates of the atoms in a simulation. They usually represent some chemically 
meaningful pathway along which advanced sampling can be performed. Listed below 
are the collective variables currently supported in SSAGES. In addition to 
specific properties for each CV, a name property can be defined for any CV 
which can be used to reference the CV from a method or other supported location. 

.. code-block:: javascript

	"name" : "mycvname"

Specified names for CVs must be unique. It is possible to omit a name and reference 
a CV by its index as well. 

Angle
-----

Example 
^^^^^^^

.. code-block:: javascript 
	
	{
		"type" : "Angle",
		"atom_ids" : [0, 1, 2]
	}

Description 
^^^^^^^^^^^

This CV calculates the bend angle, in radians, formed between three selected atoms :math:`i,j,k`,

.. math::

	\xi = \cos^{-1}\left(\frac{\mathbf{r}_{ij} \cdot \mathbf{r}_{kj}}{\Vert \mathbf{r}_{ij} \Vert \Vert \mathbf{r}_{kj} \Vert} \right).

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript 
	
	"type"

Property ``type`` must be set to string ``"Angle"``.


.. code-block:: javascript 
	
	"atom_ids"

Property ``atom_ids`` must contain three integers consisting of the atom ID forming the angle of interest.

Box Volume
----------
.. warning:: 

	Non-orthorhombic boxes are currently not supported. 

.. note:: 

	Currently supported only in Gromacs and LAMMPS.

Example 
^^^^^^^

.. code-block:: javascript 
	
	{
		"type" : "BoxVolume"
	}

Description 
^^^^^^^^^^^

This CV calculates the box volume as the determinant of the Parrinello-Rahman matrix :math:`\mathbf{H}`,

.. math::

	\xi = \det\left( H_{ij} \right)

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript 
	
	"type"

Property ``type`` must be set to string ``"BoxVolume"``.

Gyration Tensor
---------------

Example 
^^^^^^^

.. code-block:: javascript 

	"type" : "GyrationTensor", 
	"atom_ids" : [1, 2, 3, 4, 5, 6, 7, 8, 9, 10], 
	"component" : "shapeaniso"

Description 
^^^^^^^^^^^

This CV calculates components of the *mass-weighted** gyration tensor defined as, 

.. math::

	S_{mn} = \frac{1}{N}\sum_{i=1}^{N}{r_m^i r_n^i}

where :math:`r_m` is the coordinate of the :math:`m^{\mathrm{th}}` atom in the interial 
frame. With eigenvalues of :math:`\lambda_x^2, \lambda_y^2, \lambda_z^2`, possible components to 
use as a CV include: 

Radius of Gyration (Squared)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. math::

	R_g^2 = \lambda_x^2 + \lambda_y^2 + \lambda_z^2

Principal Moment
~~~~~~~~~~~~~~~~

.. math::

	\lambda_i^2,\ i \in \{x,y,z\}

Asphericity
~~~~~~~~~~~

.. math::

	b = \lambda_z^2 - \frac{1}{2}\left(\lambda_x^2 + \lambda_y^2 \right)

Acylindricity
~~~~~~~~~~~~~

.. math::

	c = \lambda_y^2 - \lambda_x^2

Shape Anisotropy
~~~~~~~~~~~~~~~~

.. math::

	\kappa^2 = \frac{3}{2}\frac{\lambda_x^4+\lambda_y^4+\lambda_z^4}{\left(\lambda_x^2+\lambda_y^2+\lambda_z^2\right)^2}-\frac{1}{2}

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript 
	
	"type"

Property ``type`` must be set to string ``"GyrationTensor"``.

.. code-block:: javascript 

	"atom_ids" 

Property ``atom_ids`` must be an array of integers containing the atom IDs which will enter the calculation. 

.. code-block:: javascript 

	"component" 

Property ``component`` must be a string defining the gyration tensor component of interest. 
Valid options are ``"Rg"``, ``"principal1"``, ``"principal2"``, ``"principal3"``, ``"asphericity"``, 
``"acylindricity"``, or ``"shapeaniso"``.

Particle Coordinate
-------------------

Example 
^^^^^^^

.. code-block:: javascript 

	{
		"type" : "ParticleCoordinate", 
		"atom_ids" : [1, 5, 6, 10],
		"dimension" : "x"
	}

Description 
^^^^^^^^^^^

This CV calculates the :math:`x`, :math:`y` or :math:`z` component center of mass of a
group of atoms. 
 
.. math::

	\xi = \frac{1}{\sum_i{m^i}}\sum_{i=1}^{N}{r_\alpha^i}\ \ \ \alpha \in {x,y,z}

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript 
	
	"type"

Property ``type`` must be set to string ``"ParticleCoordinate"``.

.. code-block:: javascript 

	"atom_ids" 

Property ``atom_ids`` must be an array of integers containing the atom IDs which will enter the calculation. 

.. code-block:: javascript 

	"dimension" 

Property ``dimension`` must be a string defining the Cartesian component of interest ``"x"``, ``"y"``, or ``"z"``.

Pairwise
--------

Example
^^^^^^^

.. code-block:: javascript 
	
	{
		"type" : "Pairwise", 
		"group1" : [1, 5],
		"group2" : [2, 3, 4, 6, 7, 8],
		"kernel" : {
			"type" : "gaussian",
			"mu" : 1.0,
			"sigma" : 0.2
		}
	}

Description 
^^^^^^^^^^^

This CV calculates a variety of pairwise properties. The functions (kernels) used are continous analogs for otherwise discontinuous CVs. If parameters are chosen judiciously, these kernels can be used in place of some standard, discontinuous CVs. A Gaussian kernel can emulate a count of nearest neighbors; a switching function kernel can emulate a coordination number.

.. math::

	\xi = \sum_{i \in A}\sum_{i \in B}{f_{ij}}

where :math:`f_{ij}` is a pairwise function for atoms :math:`i` and :math:`j`. are at a distance of the center of the Gaussian, :math:`r_{ij}=\mu`, and decreases to zero as the distance deviates away from :math:`\mu`.


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript 
	
	"type"

Property ``type`` must be set to string ``"Pairwise"``.

.. code-block:: javascript 
	
	"group1"

Property ``group1`` must be an array of integers containing the atom IDs in the first set.

.. code-block:: javascript 
	
	"group2"

Property ``group2`` must be an array of integers containing the atom IDs in the second set.

.. note:: 

	Atoms can exist in both ``group1`` and ``group2`` simultaneously. Contacts are automatically 
	skipped if :math:`i = j`.

.. code-block:: javascript 
	
	"kernel"

Property ``kernel`` must be an object defining the properties of the pairwise kernel function and its associated properties.

Pairwise Kernels
~~~~~~~~~~~~~~~~

Gaussian Function
*****************

The Gaussian function is defined as:

.. math::

	g_{ij} = e^{-\frac{\left(r_{ij} - \mu\right)^2}{2\sigma^2}}.

Properties
++++++++++

.. code-block:: javascript
	
	"mu"

Property ``mu`` is required and must be numeric.

.. code-block:: javascript 
	
	"sigma"

Property ``sigma`` is required and must be numeric.

Rational Switching Function
***************************

The rational switching function is defined as: 

.. math::

	s_{ij} = \frac{1-\left(\frac{r_{ij} - d_0}{r_0}\right)^n}{1-\left(\frac{r_{ij} - d_0}{r_0}\right)^m}.

Properties
++++++++++

.. code-block:: javascript 
	
	"type"

Property ``type`` must be set to string ``"rationalswitch"``.

.. code-block:: javascript 
	
	"d0"

Property ``d0`` is required and must be numeric. 

.. code-block:: javascript 
	
	"r0"

Property ``r0`` is required and must be numeric. 

.. code-block:: javascript 
	
	"n"

Property ``n`` is required and must be an integer. 

.. code-block:: javascript 
	
	"m"

Property ``m`` is required and must be an integer. 

Particle Position
-------------------

Example 
^^^^^^^

.. code-block:: javascript 

	{
		"type" : "ParticlePosition", 
		"atom_ids" : [1, 5, 6, 10],
		"fix" : [true, false, true],
		"position" : [3.51, 6.66, 2.14]
	}

Description 
^^^^^^^^^^^

This CV calculates the distance of the center of mass of a group of atoms 
from a particular point in Cartesian space. 

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript 
	
	"type"

Property ``type`` must be set to string ``"ParticlePosition"``.

.. code-block:: javascript 

	"atom_ids" 

Property ``atom_ids`` must be an array of integers containing the atom IDs which 
will enter the calculation. 

.. code-block:: javascript 

	"position" 

Property ``position`` must be a 3-dimensional array of numbers defining the reference 
point in the simulation box. 

.. code-block:: javascript 

	"fix" 

Property ``fix`` must be a 3-dimensional array of booleans specifying the components 
of the distance vector to include in the calculation.

Particle Separation
-------------------

Example
^^^^^^^

.. code-block:: javascript

    {
        "type" : "ParticleSeparation",
        "group1" : [1],
        "group2" : [5, 6, 10]
    }

Description
^^^^^^^^^^^

This CV calculates the distance between the centers of mass of two groups of
atoms.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

    "type"

Property ``type`` must be set to string ``"ParticleSeparation"``.

.. code-block:: javascript

    "group1"

Property ``group1`` must be an array of integers containing the atom ID(s) which
make up the first group of atoms. The CV will calculate the distance between
the center of mass of this group and the group defined by property ``group2``.

.. code-block:: javascript

    "group2"

Property ``group2`` must be an array of integers containing the atom ID(s) which
make up the second group of atoms. The CV will calculate the distance between
the center of mass of this group and the group defined by property ``group1``.

Optional
~~~~~~~~

.. code-block:: javascript

    "dimension"

Property ``dimension`` is a 3-dimensional array of booleans specifying which
Cartesian components to include in the calculation. If left unspecified, all
three xyz components will be used.

Polymer Rouse Modes
-------------------

Example
^^^^^^^

.. code-block:: javascript

    {
        "type": "RouseMode",
        "mode": 1,
        "groups":  [
                    [ 1, 2, 3, 4, 5],
                    [ 6, 7, 8, 9,10],
                    [11,12,13,14,15],
                    [16,17,18,19,20],
                    [21,22,23,24,25],
                    [26,27,28,29,30],
                    [31,32,33,34,35],
                    [36,37,38,39,40],
                    [41,42,43,44,45],
                    [46,47,48,49,50]
                   ]
    }

Description
^^^^^^^^^^^

This CV calculates the magnitude of a given Rouse mode for a set of atoms as

.. math::

    X_p = \sqrt{\mathbf{X}_p\cdot\mathbf{X}_p},

with the :math: `p` th Rouse mode defined as

.. math::

    \mathbf{X}_p = \sqrt{\frac{c_p}{N}}\sum_{i=1}^N \mathbf{R}_i \cos \Bigl[\frac{p\pi}{N}\bigl(i-\frac{1}{2}\bigr) \Bigr],

where :math: `N` is the number of groups or beads comprising the polymer, :math: `\mathbf{R}_i` is the center-of-mass of the :math: `i` th bead, and :math: `c_p` is a constant equal to 1 for :math: `p=0` and equal to 2 for :math: `p=1,\cdots,N-1`.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
^^^^^^^^

.. code-block:: javascript

    "type"

Property ``mode`` must be set to string ``"RouseMode"``.

.. code-block:: javascript

    "groups"

Property ``groups`` is an array of arrays containing the atom IDs (as integers) that comprise the discretized polymer beads. The number of groups provided implicitly defines :math: `N`, the number of polymer beads.

.. code-block:: javascript

    "mode"

Property ``mode`` is an integer indicating the index of the desired Rouse mode. Valid values range from 0 up to one less than the number of groups, or `0,\cdots, N-1`. 

Torsional Angle
---------------

Example 
^^^^^^^

.. code-block:: javascript 

	{
		"type" : "Torsional", 
		"atom_ids" : [1, 5, 6, 10]
	}

Description 
^^^^^^^^^^^

This CV calculates the dihedral angle, in radians, formed by four atoms :math:`i,j,k,l`.
It is computed as, 

.. math:: 

	\xi = \tan^{-1}\left( \frac{\left[(r_{lk} \times r_{jk}) \times (r_{ij} \times r_{jk}) \right] \cdot \frac{r_{jk}}{\Vert r_{jk}\Vert}}{(r_{lk} \times r_{jk}) \cdot (r_{ij} \times r_{jk}) } \right).

Specifically, the function ``atan2`` is used for the inverse tangent calculation to yield a four-quadrant angle.


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript 
	
	"type"

Property ``type`` must be set to string ``"Torsional"``.

.. code-block:: javascript 

	"atom_ids" 

Property ``atom_ids`` must be an array of 4 integers containing the atom IDs which 
form the dihedral. 

Alpha Helix RMSD
----------------

Example
^^^^^^^

.. code-block:: javascript

	{
            "type" : "AlphaRMSD",
            "residue_ids" : [3, 21],
            "reference" : "reference_structure.pdb",
            "unitconv" : 10
	}

Description
^^^^^^^^^^^

This CV calculates alpha helix character by comparision to an "ideal" alpha
helix structure composed of 6 amino acids. This is computed by performing a
summation over all possible sequences of 6 consecutive amino acids in the
segment of interest:

.. math::

	\xi = \sum_i \frac{1 - \left(\frac{r_i}{0.1\text{ nm}}\right)^8}{1 - (\frac{r_i}{0.1\text{ nm}})^{12}}

where :math:`r_i` is the pairwise RMSD calculated between the backbone atoms in
the 6 amino acid sequence and the ideal reference structure. 5 backbone atoms
are used for each amino acid, so each pairwise RMSD is calculated between two
sets of 30 atoms. In the case of glycine, the HA1 atom is used in place of CB
backbone atom.

.. note::

	Note that this CV is basically a summation of a switching function; in the
	future the user will be able to choose custom parameters for the switching
	function.


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

	"type"

Property ``type`` must be set to string ``"AlphaRMSD"``.

.. code-block:: javascript

	"residue_ids"

Property ``residue_ids`` must be an array of two integers designating the range
of amino acids for which to calculate the CV. The indices of the amino acids
must match those from the reference structure provided in the property
``reference``. The smaller index must be listed first, and the range must span
at least 6 amino acids.

.. code-block:: javascript

    "reference"

Property ``reference`` must be a string containing the name of a reference pdb
structure. This reference pdb structure is used along with the residue range
defined in ``residue_ids`` to check for alpha helix character. For now, all
residues in the system must be numbered in increasing order, even if they belong
to separate chains. For example, if your system has two chains of 20 amino acids
each, the first amino acid in the second chain should be numbered 21. 

Optional
~~~~~~~~

.. code-block:: javascript

    "unitconv"

Property ``unitconv`` must be numeric. This factor is used to reconcile the
internal MD units for your engine and the units used in the ideal alpha helix
reference structure. If your engine uses units of nanometers, this
can be ignored. Otherwise, ``unitconv`` must be set to the equivalent number of
length units in your MD engine equal to 1 nm. For example, if your default unit
length is in angstroms, ``unitconv`` will be set to 10. 

Anti Beta RMSD
----------------

Example
^^^^^^^

.. code-block:: javascript

	{
            "type" : "AntiBetaRMSD",
            "residue_ids" : [3, 21],
            "reference" : "reference_structure.pdb",
            "unitconv" : 10,
            "mode" : 0
	}

Description
^^^^^^^^^^^

This CV calculates anti beta-sheet character by comparision to an "ideal" anti
beta-sheet structure composed of 6 amino acids. This is computed by performing a
summation over all possible sequences of 6 amino acids, consisting of two
segments of 3 consecutive amino acids each, in the region of interest.

.. math::

	\xi = \sum_i \frac{1 - \left(\frac{r_i}{0.1\text{ nm}}\right)^8}{1 - (\frac{r_i}{0.1\text{ nm}})^{12}}

where :math:`r_i` is the pairwise RMSD calculated between the backbone atoms in
the 6 amino acid sequence and the ideal reference structure. 5 backbone atoms
are used for each amino acid, so each pairwise RMSD is calculated between two
sets of 30 atoms. In the case of glycine, the HA1 atom is used in place of CB
backbone atom.

.. note::

	Note that this CV is basically a summation of a switching function; in the
	future the user will be able to choose custom parameters for the switching
	function.


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

	"type"

Property ``type`` must be set to string ``"AntiBetaRMSD"``.

.. code-block:: javascript

	"residue_ids"

Property ``residue_ids`` must be an array of two integers designating the range
of amino acids for which to calculate the CV. The indices of the amino acids
must match those from the reference structure provided in the property
``reference``. The smaller index must be listed first, and the range must span
at least 6 amino acids.

.. code-block:: javascript

    "reference"

Property ``reference`` must be a string containing the name of a reference pdb
structure. This reference pdb structure is used along with the residue range
defined in ``residue_ids`` to check for anti beta-sheet character. For now, all
residues in the system must be numbered in increasing order, even if they belong
to separate chains. For example, if your system has two chains of 20 amino acids
each, the first amino acid in the second chain should be numbered 21. 

Optional
~~~~~~~~

.. code-block:: javascript

    "unitconv"

Property ``unitconv`` must be numeric. This factor is used to reconcile the
internal MD units for your engine and the units used in the ideal anti
beta-sheet reference structure. If your engine uses units of nanometers, this
can be ignored. Otherwise, ``unitconv`` must be set to the equivalent number of
length units in your MD engine equal to 1 nm. For example, if your default unit
length is in angstroms, ``unitconv`` will be set to 10. 

.. code-block:: javascript

    "mode"

Property ``mode`` is an integer specifying whether to calculate beta-sheets
formed only between residues on the same chain (intra) or only between residues
on separate chains (inter). If ``mode`` is set to 0, both modes will be used.
A value of 1 selects for the intra mode; a value of 2 selects for inter mode.


Parallel Beta RMSD
------------------

Example
^^^^^^^

.. code-block:: javascript

	{
            "type" : "ParallelBetaRMSD",
            "residue_ids" : [3, 21],
            "reference" : "reference_structure.pdb",
            "unitconv" : 10,
            "mode" : 0
	}

Description
^^^^^^^^^^^

This CV calculates anti beta-sheet character by comparision to an "ideal"
parallel beta-sheet structure composed of 6 amino acids. This is computed by 
performing a summation over all possible sequences of 6 amino acids, consisting
of two segments of 3 consecutive amino acids each, in the region of interest.

.. math::

	\xi = \sum_i \frac{1 - \left(\frac{r_i}{0.1\text{ nm}}\right)^8}{1 - (\frac{r_i}{0.1\text{ nm}})^{12}}

where :math:`r_i` is the pairwise RMSD calculated between the backbone atoms in
the 6 amino acid sequence and the ideal reference structure. 5 backbone atoms
are used for each amino acid, so each pairwise RMSD is calculated between two
sets of 30 atoms. In the case of glycine, the HA1 atom is used in place of CB
backbone atom.

.. note::

	Note that this CV is basically a summation of a switching function; in the
	future the user will be able to choose custom parameters for the switching
	function.


Options & Parameters
^^^^^^^^^^^^^^^^^^^^

Required
~~~~~~~~

.. code-block:: javascript

	"type"

Property ``type`` must be set to string ``"ParallelBetaRMSD"``.

.. code-block:: javascript

	"residue_ids"

Property ``residue_ids`` must be an array of two integers designating the range
of amino acids for which to calculate the CV. The indices of the amino acids
must match those from the reference structure provided in the property
``reference``. The smaller index must be listed first, and the range must span
at least 6 amino acids.

.. code-block:: javascript

    "reference"

Property ``reference`` must be a string containing the name of a reference pdb
structure. This reference pdb structure is used along with the residue range
defined in ``residue_ids`` to check for parallel beta-sheet character. For now,
all residues in the system must be numbered in increasing order, even if they 
belong to separate chains. For example, if your system has two chains of 20 
amino acids each, the first amino acid in the second chain should be numbered 
21. 

Optional
~~~~~~~~

.. code-block:: javascript

    "unitconv"

Property ``unitconv`` must be numeric. This factor is used to reconcile the
internal MD units for your engine and the units used in the ideal parallel
beta-sheet reference structure. If your engine uses units of nanometers, this
can be ignored. Otherwise, ``unitconv`` must be set to the equivalent number of
length units in your MD engine equal to 1 nm. For example, if your default unit
length is in angstroms, ``unitconv`` will be set to 10. 

.. code-block:: javascript

    "mode"

Property ``mode`` is an integer specifying whether to calculate beta-sheets
formed only between residues on the same chain (intra) or only between residues
on separate chains (inter). If ``mode`` is set to 0, both modes will be used.
A value of 1 selects for the intra mode; a value of 2 selects for inter mode.



