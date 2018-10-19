.. SSAGES documentation master file, created by
   sphinx-quickstart on Mon Jun  6 16:05:35 2016.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to SSAGES!
==================

In their simplest form, particle-based simulations are limited to generating
ensembles of configurations (in Monte Carlo [MC] simulations) or trajectories
in time (in molecular dynamics [MD] or Brownian dynamics [BD]). One can then
extract mechanical variables such as the potential energy or pressure and
perform ensemble or time averages. There are two important limitations to such
calculations:

1. For complex materials, the time scales available to standard
   MD simulations are often insufficient to sample relevant regions of phase
   space
2. In order to develop a fundamental understanding of materials,
   researchers are primarily interested in calculating the free energy, the
   entropy, and their derivatives with respect to various thermodynamic
   quantities (which lead to material properties such as elastic moduli, heat
   capacity, and various other susceptibilities).

These quantities are difficult to obtain or intractable in standard MC and MD
simulations. To overcome these limitations, MC and MD simulations must be
supplemented with advanced sampling techniques. These methods are critical for
the efficient simulation of complex assembly processes.

SSAGES (Software Suite for Advanced General Ensemble Simulations) is
designed to perform these calculations. The framework is designed to treat
molecular simulation routines as a black box, using the coordinates of the
system as evolved by an MD engine to compute collective variables which
permit a meaningful reduced-dimensionality representation of the phase space
within a system. This information is then used to define evolving reactive
pathways or to bias the statistics of a simulation for the purposes of
computing free energies. The internal structure of the code has been designed
to be simple and extensible to new sampling methods and engines. For further
details on examples and capabilities of SSAGES, peruse the documentation
for :ref:`specific methods <advanced-sampling-methods>`.

Contents:

.. toctree::
   :maxdepth: 2

   Introduction
   Getting Started
   Input Files
   Engines
   Collective Variables
   Advanced Sampling Methods
   Write your own Methods and CVs
   Contribute to SSAGES
   The SSAGES Cookbook
   Acknowledgements
   Copyright and License

Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

