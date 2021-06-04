SSAGES v0.9 Release Notes
=========================

## v0.9.2
- Critical bug fix for patching Eigen with GROMACS
- Correction to gradient of AngleCV
- Additional checks and docs for multiwalker loggers
- Support newer LAMMPS versions
- Minor source cleanup

## v0.9.1
- Critical bug fix for neural-network based methods
- Updated Eigen dependency to v3.3.7
- Many documentation improvements
- Support for newer HOOMD versions (2.6+)
- Build process fixes

## v0.9.0
- New Combined Force–Frequency sampling method
- Addition of non-weighted internal center of mass calculation
- Documentation additions
- New ANN-based collective variable

SSAGES v0.8 Release Notes
=========================

## v0.8.6
- New RMSD Collective Variable
- Improved Qbox examples
- Standardized test suite
- Improved documentation
- Addition of Travis-CI testing
- General code cleanup and bugfixes (See commits)

## v0.8.5
- Improved Elastic Band Sampling
- Temporary removal of GROMACS 2019
- Additional Documentation
- General code cleanup and bugfixes (See commits)

## v0.8.4
- Major documentation updates!
- Improved HOOMD-blue support
- ANN restarts
- Python script for Metadynamics example
- Support for newer GROMACS versions (up to 2019.1)
- Several bugfixes (See commits)

## v0.8.3
- HOOMD-blue support!
- Support for newer GROMACS and LAMMPS versions
- CV definition checking for Methods
- Documentation updates
- Several bugfixes (See commits)

## v0.8.2
- Grid internal updates
- ABF Integrator handles interpolation in each direction independently
- ABF restarts now handled with JSON member
- GyrationTensorCV can be projected into any number of dimensions
- Documentation updates
- Several bugfixes (See commits)

## v0.8.1
- GROMACS support for all 5.1.x, 2016.x, 2018.x!
- Better CMake handling for Hooks
- Handling of LAMMPS line continuations (&)
- Correct handling of multiprocessor ABF method
- BFS method cleanup
- Minor documentation updates
- Eigen include update (3.3.4)
- googletest include update (1.8.0)
- jsoncpp include update

## v0.8.0
- Added ANN sampling!
- More documentation updates
- Updates to examples
- Added Fourier and Chebyshev basis sets to BFS
- Added 3D ABF integrator
- Improved periodicity handling on grid
- Secondary structure CV bug fixes
- Improved Qbox integration


SSAGES v0.7 Release Notes
=========================

## v0.7.5
- Major documentation update!
- Major examples update!
- New generalized pairwise CV
- New secondary structure CVs (alpha and anti/parallel beta sheet RMSD)
- New CV logging capability
- CV selection by name in methods
- Updated unit tests
- Updated ABF integrators for periodic and non-periodic CVs
- Gromacs 2016.3 support!
- Many bug fixes!

## v0.7.0
- New simplified JSON syntax
- Support for multiple simultaneous methods!
- Eliminated boost dependency!
- CV selector for methods
- Argument forwarding for Gromacs
- Updated forward flux examples
- Significant under-the-hood improvements
- Fixed Gromacs auto-download


SSAGES v0.6 Release Notes
=========================

## v0.6.0
- Support for QBox first-principles MD engine
- Support for OpenMD engine
- Coordination number CV
- Polymer Rouse modes CV
- Box volume CV
- Virial contribution (NPT support) for some CVs and methods
- Updated examples and documentation
- New backend grid
- Grid-based metadynamics
- Updated forward flux sampling
- Fixed regression with string methods
- Performance and other improvements!


SSAGES v0.5 Release Notes
=========================

## v0.5.0
- Gromacs restart support
- New gyration tensor CVs
- Updated examples and documentation
- Metadynamics optimizations
- Better engine error handling
- More! (See commit log)
