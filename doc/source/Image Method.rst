:orphan:

.. image_method:

Image Method
------------

Introduction
^^^^^^^^^^^^

Surface charging or polarization can strongly affect the nature of interactions
between charged dielectric objects, particularly when sharp dielectric
discontinuities are involved. However, not any efficient and accurate
computation tools are publicly available especially for the description of
polarization effects in many-body systems. 

For this purpose, Image Method, an analytic perturbative approach we recently
developed for evaluating the polarization energy of a many-body collection of
charged dielectric spheres embedded in a dielectric medium becomes particularly
suitable [1]_.

The polarization-induced interactions between these spheres depend on the ratio
of dielectric constants for the spheres and the medium, and the ratio of the
distance between particles and the radii of the particles. We have shown that,
in some cases, polarization completely alters the qualitative behavior, and in
some other cases, polarization leads to stable configurations that otherwise
could not occur in its absence. 

We think it is helpful to include Image Method into SSAGES for users to include
polarization corrections properly in their systems, and meanwhile, to couple
with advanced sampling methods to accelerate their simulations. 

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

SSAGES Image method is implemented in a way that is as easy as conducting a
simulation using LAMMPS that only includes pairwise Coulombic interactions into
electrostatic interactions. To achieve this, we update the electrostatic forces
acting on all objects by adding up the polarization corrections using SSAGES
engine and then pass the modified snapshot back to LAMMPS engine at each time
step. The JSON file needed for SSAGES engine should include:

einner
    The relative dielectric permittivity of polarizable object. 

ion-type-start
    For cases that you have both polarizable objects and non-polarizable objects
    in you system, for example, in which colloids and ions are treated as
    polarizable and non-polarizable, respectively. This parameter controls where
    the non-polarizable typos start. 

atom type radius
    Radius of all types of objects. 

Guidelines
^^^^^^^^^^

It is very similar as running a simulation including electrostatic interactions
using LAMMPS. Referring to the exampled LAMMPS INPUTFILE and DATAFILE, you need
to double check you have declared the following variables that are particularly
necessary for Image Method to compute polarization corrections: 

* charges
* dielectric (relative dielectric permittivity of the surrounding continuum)

Method Output
^^^^^^^^^^^^^

There are not special outputs files generated for Image method since it only
provides an updated electrostatic forces by including polarization corrections.
Nevertheless, we provided options of dumping trajectories and printing out
force-distance data in the LAMMPS INPUTFILE examples for users to visualize how
significant the polarization effects are in some cases more conveniently. 

.. _IM_tutorial:

Tutorial
^^^^^^^^

.. todo::

    Write a tutorial. 

Developer
^^^^^^^^^

* Jiyuan Li

References
^^^^^^^^^^

.. [1] J. Qin, J. Li, V. Lee, H. Jaeger, J. J. de Pablo, and K. Freed,
       *A theory of interactions between polarizable dielectric spheres*,
       J. Coll. Int. Sci. **469**, 237 - 241 (2016)
