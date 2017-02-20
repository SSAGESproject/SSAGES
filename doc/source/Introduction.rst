Introduction
============

Welcome to SSAGES, our new and shiny advanced sampling package.

Over the past several decades, molecular simulation has emerged as a powerful
tool for investigating a wide range of physical phenomena. Molecular simulation
is, in essence, a computational "microscope" whereby computers are used to "look
at" the properties of a system that are difficult to observe or measure through
traditional experimental setups. The comparison between simulations and the
corresponding experimental systems can sometimes be challenging, usually due to
factors such as the length and time scales explored. In simulation, a molecular
model must have sufficient temporal and spatial accuracy to resolve the fastest
time scales and shortest length scales within a system. Unfortunately, due to
computational constraints, this detailed resolution has limited the length of
time and number of particles that a model can simulate, typically simulating
systems that are smaller than analogous experimental setups in laboratory
environments for much shorter times than the duration of the experiments.
However in recent years, advancements in computational processing power,
including custom built computer architectures, have continued to increase the
time and length scales accessible by molecular simulation, with current
state-of-the-art simulations able to analyze systems for milliseconds (10-3s). 

Another challenge arises from the difficulty in obtaining good statistics from
molecular simulations.  Thermal fluctuations dominate motion at the nano-scale
and result in motion that appears random (i.e. Brownian), with no two molecular
trajectories being identical. As a result, statistically meaningful averages are
necessary in order to calculate thermodynamic and kinetic quantities of interest
in these systems.  An incredibly powerful thermodynamic quantity referred to as
the relative free energy of a system can be calculated in this way. The relative
free energy can characterize underlying system behavior in the presence of the
thermal-induced random noise. Performing this necessary averaging within
simulations is challenging. In essence, the requirement of averaging compounds
the issue of time scales described previously; not only must long simulations be
performed, but they must be performed a prohibitively large number of times in
order to extract sufficient statistics. It is therefore necessary to develop
efficient techniques to calculate meaningful averages from simulations.

Advanced sampling methods represent a class of simulation techniques that seek
to improve this improper averaging and accelerate the extraction of useful
properties (e.g. free energies, transition paths) from simulations.  At the
heart of all advanced sampling methods, is statistical mechanics, a field of
physics that relates microscopic phenomena (i.e. the motion of particles) to
macroscopic observables (e.g. temperature and pressure). By taking advantage of
statistical mechanics, advanced sampling methods are used to apply a systematic
bias to a simulation to speed convergence, and then mathematically remove this
bias to extract the true underlying behavior. Throughout the past decade,
advanced sampling methods have become wildly successful, and have now become an
essential component in the toolbox of molecular simulation. 

Despite the demonstrated utility of advanced sampling techinques, they have only
been adopted by a fraction of the scientists working in the field. One
explanation for this slow adoption is technical: advanced sampling methods are
complicated, and not all research groups have the expertise required in order to
implement these methods themselves. In the worst case, this leads to long stages
of code development, possibly leading to unknown implementation errors or
insufficient validation. Even in cases when advanced sampling methods are
implemented, they are typically done so with a specific problem in mind and are
custom-built for a certain model or application. This specificity necessitates
modification of the custom-built advanced sampling code when studying new
systems. This prevents the distribution of code between researches in the field.
As a result, the same methods are implemented again and again by different
members of the community. Sadly, in molecular simulation, it is quite common to
"reinvent the wheel". 

SSAGES is an answer to this problem. SSAGES (Suite for Advanced Generalized
Ensemble Simulations) is a free, open-source software package that allows users
to easily apply advanced sampling techniques to any molecular system of
interest. Simply put, SSAGES is a wrapper that converts a molecular simulation
engine (e.g. LAMMPS, NAMD) into an advanced sampling engine. SSAGES contains a
library of widely used enhanced sampling methods that can be used to calculate
everything from free energies to transition pathways. Importantly, SSAGES works
with many of the widely used simulation packages, and can simply be added on top
of the simulations a researcher is already running. SSAGES is implemented in a
highly modular way, and is easily extended to incorporate a new method or to
modify an existing one and has been rigorously tested to ensure the accuracy of
its calculations. 

In short, SSAGES makes advanced sampling methods easy. We hope that it will do
just that for your research.

