Forward-Flux
------------

Forward Flux Sampling (FFS) is a specialized method to simulate “rare events” in
non-equilibrium and equilibrium systems with stochastic dynamics. Several review
articles in the literature present a comprehensive perspective on the basics,
applications, implementations, and recent advances of FFS Here, we provide a
brief general introduction to FFS, and describe the Rosenbluth-like variant of
forward flux method, which is implemented in SSAGES. We also explain various
options and variables to setup and run an efficient FFS simulation using SSAGES.

Introduction
^^^^^^^^^^^^

Rare events occur infrequently in nature mainly due to significant activation
energies necessary to take a system from an initial state (commonly referred to
as “A”) to a final state (or state “B”). The outcomes of rare events are
generally substantial and thereby it is essential to obtain a molecular-level
understanding of the mechanisms and kinetics of these events. Examples of
small-scale rare events include conventional nucleation and growth phenomenon,
folding/unfolding of large proteins, and non-spontaneous chemical reactions.
“Thermal fluctuations” commonly drives the systems from an initial state to a
final state over an energy barrier :math:`\Delta E`. The transition frequency
from state A to state B is proportional to :math:`e^{\frac{-\Delta E}{k_{B}T}}`,
where :math:`k_{B}T` is the thermal energy of the system. Accordingly, the time
required for an equilibrated system in state A to reach state B grows
exponentially (at a constant temperature) as the energy barrier :math:`\Delta E`
become larger. Eventually none or only a few transitions may occur within the
typical timescale of molecular simulations. In FFS method several intermediate
states or so-called interfaces (:math:`\lambda_{i}`) are placed along a
“reaction coordinate” or an “order parameter” between the initial state A and
the final state B (Figure 1). These intermediate states are chosen such that the
energy barrier between adjacent interfaces are readily surmountable using
typical simulations. Using the stored configurations at an interface, several
attempts are made to arrive at the next interface in the forward direction (the
order parameter must increase monotonically when going from A to B). This
incremental progress makes it more probable to observe a full transition path
from state A to state B. FFS uses positive flux expression to calculate rate
constant. The system dynamics are integrated forward in time and therefore
detailed balance is not required.

.. figure:: images/forward_flux_image1.png
    :align: center

    In Forward Flux sampling method, dimensionality of the system is reduced by
    choosing one or more “reaction coordinate” or “order parameter”. Several
    equally-spaced intermediate states are placed along the order parameter to
    link the initial state A and the final state B. Incremental progress of the
    system is recorded and analyzed to obtain relevant kinetic and thermodynamic
    properties.

Several protocols of forward flux method have been adopted in the literature to

1) generate the intermediate configurations,
2) calculate the conditional probability of reaching state B starting from
   state A, :math:`P(\lambda_{B} = \lambda_{n} | \lambda_{A} = \lambda_{0})`,
3) compute various thermodynamic properties, and
4) optimize overall efficiency of the method. The followings are the widely-used
   variants of forward flux sampling method:

* Direct FFS (DFFS)
* Branched Growth FFS (BGFFS)
* **Rosenbluth-like FFS (RBFFS)**
* Restricted Branched Growth FFS (RBGFFS)
* FFS Least-Squares Estimation (FFS-LSE)
* FF Umbrella Sampling (FF-US) 

Rosenbluth-like FFS (RBFFS) has been implemented in the current version of
SSAGES. Direct FFS (DFFS) and Branched growth FFS (BGFFS) will be included in
the future release of SSAGES.

Rate Constant and Initial Flux
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

The overall rate constant or the frequency of going from state A to state B is
computed using the following equation:

.. math::

    k_{AB} = \Phi_{A,0} \cdot P\left(\lambda_{N} \vert \lambda_{0}\right)

here, :math:`\Phi_{A,0}` is the initial forward flux or the flux at the initial
interface, and :math:`P\left(\lambda_{N} \vert \lambda_{0}\right)` is the
conditional probability of the trajectories that initiated from A and reached B
before returning to A. In practice, :math:`\Phi_{A,0}` can be obtained by
simulating a single trajectory in State A for a certain amount of time
:math:`t_{A}`, and counting the number of crossing of the initial interface
:math:`\lambda_{0}`. Alternatively, a simulation may be carried out around state
A for an unlimited period of time until :math:`N_{0}` number of accumulated
checkpoints is stored (this has been implemented in SSAGES):

.. math::

    \Phi_{A,0} = \frac{N_{0}}{t_{A}}

here, :math:`N_{0}` is the number of instances in which :math:`\lambda_{0}` is
crossed, and :math:`t_{A}` is the simulation time that the system was run around
state A. Note that

1) :math:`\lambda_{0}` can be crossed in either forward
   (:math:`\lambda_{t} < \lambda_{0}`) or backward
   (:math:`\lambda_{t} > \lambda_{0}`) directions, but only "forward crossing"
   marks a checkpoint (see Figure 2) and
2) :math:`t_{A}` should only include the simulation time around state A and
   thereby the portion of time spent around state B must be excluded, if any. 

In general, the conditional probability  is computed using the following expression:

.. math::

    P\left(\lambda_{n} \vert \lambda_{0}\right) =
    \prod\limits_{i=0}^{n-1} P\left(\lambda_{i+1} \vert \lambda_{i}\right) =
    P\left(\lambda_{1}\vert\lambda_{0}\right)\cdot P\left(\lambda_{2}\vert\lambda_{1}\right)
    \dots P\left(\lambda_{n}\vert\lambda_{n-1}\right)

:math:`P\left(\lambda_{i+1}\vert\lambda_{i}\right)` is computed by initiating a
large number of trials from the current interface and recording the number of
successful trials that reaches the next interface. The successful trials in
which the system reaches the next interface are stored in the memory and used as
checkpoints in the next interface. The failed trajectories that go all the way
back to state A are terminated. Different flavors of forward flux method use
their unique protocol to select checkpoints to initiate trials at a given
interface, compute final probabilities, create transitions paths, and analyze
additional statistics.

.. figure:: images/forward_flux_image2.png
    :align: center

    A schematic representation of computation of initial flux using a single
    trajectory initiated in state A. The simulation runs for a certain period of
    time :math:`t_{A}` and number of forward crossing is recorded. Alternatively,
    we can specify the number of necessary checkpoints :math:`N_{0}` and run a
    simulation until desired number of checkpoints are collected. In this figure,
    green circles show the checkpoints that can be used to generate transition
    paths.

Rosenbluth-like Forward Flux Sampling (RBFFS)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Rosenbluth-like Forward Flux Sampling (RBFFS) method is an adaptation of
Rosenbluth method in polymer sampling to the simulation of rare events [4]_.
The RBFFS is comparable to Branched Growth Forward Flux (BGFFS) [1]_ [2]_ but,
in contrast to BGFFS, a single checkpoint is randomly selected at a non-initial
interface instead of initiation of trials from all checkpoints at a given
interface (Figure 3). In RBFFS, first a checkpoint at :math:`\lambda_{0}` is
selected and :math:`k_{0}` trials are initiated. The successful runs that reach
:math:`\lambda_{1}` are stored and the rest that go back to A are terminated.
Next, one of the checkpoints at :math:`\lambda_{1}` is randomly chosen (in
contrast to Branched Growth where all checkpoints are involved), and
:math:`k_{1}` trials are initiated to :math:`\lambda_{2}`. Last, this procedure
is continued for the following interfaces until state B is reached or all trials
fail. This algorithm is then repeated for the remaining checkpoints at
:math:`\lambda_{0}` to generate multiple “transition paths”.

.. figure:: images/forward_flux_image3.png
    :align: center

    Rosenbluth-like Forward Flux Sampling (RBFFS) involves sequential generation
    of unbranched transition paths from all available checkpoints at the first
    interface :math:`\lambda_{0}`. A single checkpoint at the interface
    :math:`\lambda_{i > 0}`  is randomly marked and :math:`k_{i}` trials are
    initiated from that checkpoint which may reach to the next interface
    :math:`\lambda_{i+1}` (successful trials) or may return to state A (failed
    trial).

In Rosenbluth-like forward flux sampling, we choose one checkpoint from each
interface independent of the number of successes. The number of available
checkpoints at an interface are not necessarily identical for different
transition paths :math:`p`. This implies that more successful transition paths
are artificially more depleted than less successful paths. Therefore, we need to
enhance those extra-depleted paths by reweighting them during post-processing.
The weight of path :math:`p` at the interface :math:`\lambda_{i}` is given by:

.. math::

    w_{i,b} = \prod\limits_{j=0}^{i-1} \frac{S_{j,p}}{k_{j}}

where :math:`S_{j,p}` is the number of successes at the interface :math:`j` for
path :math:`p`. The conditional probability is then computed using the following
expression:

.. math::

    P\left(\lambda_{n}\vert\lambda_{0}\right) =
    \prod\limits_{i=0}^{n-1} P\left(\lambda_{i+1} \vert \lambda_{i}\right) =
    \frac{ \prod_{i=0}^{n-1}\sum_{p} w_{i,p} S_{i,p} / k_{i} }{ \sum_{p} w_{i,p} }

:math:`\Sigma` here runs over all transition paths in the simulation.

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

To run a RBFFS simulation using SSAGES, an input file in JSON format is required
along with a general input file designed for your choice of molecular dynamics
engine (MD engine). For your convenience, two files ``Template_Input.json`` and
``FF_Input_Generator.py`` are provided to assist you in generating the JSON
file. Here we describe the parameters and options that should be set in
``Template_Input.json`` file in order to successfully generate an input file and
run a RBFFS simulation.

Input and parameters related to "driver"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

type 
    + Type: string
    + Default:  “LAMMPS”
    + Functionality:  Defines the preferred MD engine for running the actual
      simulation. You are encouraged to read the documentation page of the
      corresponding MD package to learn about input files and different options
      of that package.   

num processors
    + Type: integer
    + Default: 1
    + Functionality:  Sets the number of processors that each individual drivers
      uses to run the simulation. In current version of SSAGES, drivers can only
      use one processor.

MDSteps
    + Type: integer
    + Default: 1000000000
    + Functionality:  Sets the maximum number of MD steps allowed for the FFS
      simulation on a given walker. We recommend defining a large number here to
      ensure that the simulation is completed before reaching that many steps.
      SSAGES will exit upon completion of the FFS simulation.

logfile
    + Type: string
    + Default: “none”
    + Functionality: Sets the name of engine-dependent log file that MD engine
      uses to write the simulation information including timesteps, energies, etc.

Input and parameters related to "method"
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

type
    + Type: string
    + Default: "ForwardFlux"
    + Functionality:  Specifies that “ForwardFlux” module of SSAGES will be
      activated. Don’t change this if you plan to run a forward flux sampling
      simulation.

index_file
    + Type: string
    + Default: none
    + Functionality: Stores interface information in the format:
      `Interface filename origin`
      The file-naming scheme is based on the interface the simulation is on and
      cumulative hash number. Origin is the filename of the file the trajectory
      was previously fired from to reach the current interface position.

library_file
    + Type: string 
    + Default: "library_input.dat"
    + Functionality:  Sets the name of the file that stores the checkpoints at
      the initial interface  by running a serial trajectory around state A.

results_file
    + Type: string 
    + Default: "results.dat"
    + Functionality: Specifies the name of the file in which the results of the
      forward flux simulation is stored. This file can later be helpful for
      post-processing purposes.

centers
    + Type: array
    + Default: none
    + Functionality:  Defines an array of intermediate interfaces that links the
      initial state A to the final state B. This array can either be defined in
      the ``Template_Input.json`` file or ``FF_Input_Generator.py`` file. In the
      latter case, the values of **centers** is left blank in the
      ``Template_Input.json`` file.

generate_configs
    + Type: integer 
    + Default: 1
    + Functionality: Defines the number of checkpoints/configurations that ought
      to be generated at the first interface, i.e. .  

shots
    + Type: integer
    + Default: 1
    + Functionality:  Sets the number of trials that should be initiated from
      the randomly selected checkpoints at an interface (at the initial
      interface, all checkpoint are used to generate multiple transition paths).
      In principle, this can change from interface to interface but in the
      current implementation of SSAGES, the number of trials/shots from a
      checkpoint/node is assumed to be a constant number.  

frequency
    + Type: integer
    + Default: 1
    + Functionality:  Specifies the frequency (in timesteps of MD simulation)
      that SSAGES recomputes the value of “order parameter” and writes the
      output data. 

restart_type
    + Type: string
    + Default: "new_library"
    + Functionality: Defines how a FFS simulation should be restarted. Several
      options are available:

        1) "new_library": generates a new starting library. If this option is
           defined, a new FFS simulation is setup and run.
        2) "from_library": restarts from a library of available configurations
           defined by library_file and library_point.
        3) "from_interface": restarts the simulation from an interface defined
           by the current position of the CV from configurations found in
           ``index_contents``.
        4) "none": SSAGES restarts the FFS simulation using snapshots of
           trajectories that are not necessarily checkpoints/nodes located at a
           specific interface.

      "from_library", "from_interface" and "none" are typically reserved for
      restarting from crashes only. 
      
library_point
    + Type: integer
    + Default: none
    + Functionality: Specifies the current library configuration that you are on
      from the list of configurations found in the library file defined by
      ``library_file``.

current_hash
    + Type: integer
    + Default: 1
    + Functionality: Used in the file-naming scheme. Mainly needed for restarts,
      or if specifying where the number scheme should start. Default is based on
      walker_ID*1000000, meaning walker 0 files will be
      ``dump_"interface"_0.dump``, ``dump_"interface"_1.dump``, etc.

index_contents
    + Type: string
    + Default: none
    + Functionality: Only used for restarts by SSAGES, includes the same
      contents as index_file.

successes
    + Type: array
    + Default: none
    + Functionality: Only used for restarts by SSAGES, contains successes on
      each interface for each library configuration explored so far. Contents
      are exactly those as ``results_file: walker_id library_point`` "list of
      successes at each interface". For example, a library consisting of two
      configurations and 4 interfaces using 1 walker:

      0 0 2 3 4 0

      0 1 1 0 0 0

current_shot
    + Type: integer
    + Default: none
    + Functionality: Mainly used for restarts, indicates which shot this walker
      is on.

Other required input parameters
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

CVs
    + Type: array
    + Default: none
    + Functionality: Selection of "order parameter" or "reaction coordinate".
      The current implementation of FFS in SSAGES can only take one collective
      variable. See section XXX for more details.

inputfile
    + Type: string
    + Default: none
    + Functionality: Specifies the name of engine-specific input file name. The
      user is encouraged to refer to the documentation page of the corresponding
      MD package to learn about various input options as well as the structure
      and format of input files suitable for MD engine of your choice.

Tutorial
^^^^^^^^

.. todo::

    Give a tutorial. The tutorial can be based on one of the examples for this
    method. Describe how to compile the input files and how to call SSAGES.
    Describe how to understand and visualize the results.

Developer
^^^^^^^^^

Hadi & Joshua Lequieu.

References
^^^^^^^^^^

.. [1] R. J. Allen, C. Valeriani, P. R. ten Wolde, *Forward Flux Sampling for
       Rare Event Simulations*. J Phys-Condens Mat 2009, 21 (46)
       
.. [2] F. A. Escobedo, E. E. Borrero, J. C. Araque, *Transition Path Sampling
       and Forward Flux Sampling. Applications to Biological Systems*.
       J Phys-Condens Mat 2009, 21 (33).

.. [3] R. J. Allen, D. Frenkel, P. R. ten Wolde, *Forward Flux Sampling-Type
       Schemes for Simulating Rare Events: Efficiency Analysis*.
       J. Chem. Phys. 2006, 124 (19).

.. [4] M. N. Rosenbluth, A. W. Rosenbluth, *Monte-Carlo Calculation of the
       Average Extension of Molecular Chains*.
       J. Chem. Phys. 1955, 23 (2), 356-359.
