<div align="center">
  <a href="http://miccomcodes.org" target="_blank">
    <img src="http://miccomcodes.org/static/img/ssageslogo.jpg" alt="SSAGES" height="200">
  </a>
</div>

<h2 align="center">
<p align="center">
  <a href="http://miccomcodes.org/manual/index.html" target="_blank">
    <img src="https://img.shields.io/badge/docs-v0.8-blue.svg" alt="Documentation">
  </a>
  &nbsp;
  <a href="https://doi.org/10.1063/1.5008853" target="_blank">
    <img src="https://img.shields.io/badge/doi-10.1063%2F1.5008853-blue" alt="Cite SSAGES">
  </a>
</p>
</h2>

**SSAGES** (**S**oftware **S**uite for **A**dvanced **G**eneral **E**nsemble
**S**imulations) is an open-source, engine agnostic, C++11 based advanced
sampling package.  It is designed to be easy to use, extendable and extremely
versatile. It is currently pre-beta, meaning that there are many rough edges,
but we are working rapidly to expand its features and fix any bugs. Keep an eye
on this page for future updates and see below on how to contribute!

## What's New (v0.8.5)
- Improved Elastic Band Sampling
- Temporary removal of GROMACS 2019
- Additional Documentation
- General code cleanup and bugfixes (See commits)

<a id="features"></a>
## Features
**SSAGES** currently works with multiple molecular dynamics engines. It contains a variety of collective variables (CVs) and advanced sampling methods. 

### Highlights 
- Engine agnostic framework 
- Simple JSON input file syntax 
- Easy to add new CVs 
- Easy to add new methods
- Much more!

### Engines 
- GROMACS 5.1.x, 2016.x, 2018.x
- LAMMPS (Most recent versions)
- OpenMD (2.5+)
- QBox (1.63+)

### CVs
- Atom group coordinate
- Atom group position 
- Atom group separation 
- Bend angle
- Box volume
- Components of gyration tensor
- Pairwise kernel (coordination number, nearest neighbors)
- Polymer Rouse modes 
- Secondary structure (alpha, anti/parallel beta sheet) RMSD
- Torsional angle

### Methods 
- Adaptive biasing force 
- Artificial neural network sampling
- Basis function sampling 
- Metadynamics 
- Umbrella sampling 
- Finite temperature string 
- Nudged elastic band 
- Swarm of trajectories 
- Forward flux sampling 

<a id="installation"></a>
## Installation
The first step is to clone the repository locally.

```bash
$ git clone https://github.com/MICCoM/SSAGES-public.git
```
**SSAGES** uses a CMake build system. It also requires the use of a support MD engine.
For example, to compile with LAMMPS, execute the following

```bash
$ cd SSAGES
$ mkdir build && cd build
$ cmake -DLAMMPS_SRC=/path/to/lammps/src .. 
$ make
```

This will build a SSAGES executable which will reside in the build directory.

If you want to use a specific compiler (or if your default compiler is not supported),
set the C and C++ compilers with `CMAKE_C_COMPILER` and `CMAKE_CXX_COMPILER`, respectively.
For example, to use gcc/g++, replace the CMake command with

```bash
$ cmake -DLAMMPS_SRC=/path/to/lammps/src -DCMAKE_C_COMPILER=gcc -DCMAKE_CXX_COMPILER=g++ ..
```

If you want to compile and run unit and integration tests, replace the cmake command
in the example above with

```bash
$ cmake -DLAMMPS_SRC=/path/to/lammps/src -DBUILD_TESTS=ON ..
```

### MPI

A requisite underlying MPI library also required to run SSAGES. 
On recent Debian based systems using OpenMPI, the requirement 
can be installed via:

```bash 
$ sudo apt-get install libopenmpi-dev openmpi-bin
```

For more detail on the build system, please check the documentation.

## Build the documentation

After you have built SSAGES you can further build the documentation from the build/
directory. Make sure you have the following tools available:

* doxygen
* dot (part of the `graphViz` package)
* sphinx (`python-sphinx` or `python3-sphinx`)
* ReadTheDocs theme (`pip install sphinx_rtd_theme`)
* sphinxcontrib-bibtex (`pip install sphinxcontrib-bibtex`)

You can build the documentation with
```bash
$ make doc
```
You can alse build the API-references and the User
Manual separately with
```bash
$ make apiref
```
and
```bash
$ make manual
```

If you have `pdflatex` installed, you can also build
a PDF file for the documentation. To compile the
API-reference into a PDF file do
```bash
$ cd doc/API-doc/latex/
$ make
```
The pdf is called refman.pdf

Similarly, you can build a PDF version of the Manual with
```bash
$ cd doc/Manual/
$ make
```
The pdf-file will be called SSAGES.pdf

### View the documentation

Once you have built the documentation you will find it
in the `doc/API-doc/` and `doc/Manual` directories. To
view the documentation in a browser call for example
```bash
$ firefox doc/API-doc/html/index.html
```
or
```bash
$ firefox doc/Manual/index.html
```

## What's Old (v0.8.4)
- Major documentation updates!
- Improved HOOMD-blue support
- ANN restarts
- Python script for Metadynamics example
- Support for newer GROMACS versions (up to 2019.1)
- Several bugfixes (See commits)

## What's Old (v0.8.3)
- HOOMD-blue support!
- Support for newer GROMACS and LAMMPS versions
- CV definition checking for Methods
- Documentation updates
- Several bugfixes (See commits)

## What's Old (v0.8.2)
- Grid internal updates
- ABF Integrator handles interpolation in each direction independently
- ABF restarts now handled with JSON member
- GyrationTensorCV can be projected into any number of dimensions
- Documentation updates
- Several bugfixes (See commits)

## What's Old (v0.8.1)
- GROMACS support for all 5.1.x, 2016.x, 2018.x!
- Better CMake handling for Hooks
- Handling of LAMMPS line continuations (&)
- Correct handling of multiprocessor ABF method
- BFS method cleanup
- Minor documentation updates
- Eigen include update (3.3.4)
- googletest include update (1.8.0)
- jsoncpp include update

## What's Old (v0.8.0)
- Added ANN sampling!
- More documentation updates
- Updates to examples
- Added Fourier and Chebyshev basis sets to BFS
- Added 3D ABF integrator
- Improved periodicity handling on grid
- Secondary structure CV bug fixes
- Improved Qbox integration

## What's Old (v0.7.5)
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

## What's Old (v0.7.0)
- New simplified JSON syntax
- Support for multiple simultaneous methods! 
- Eliminated boost dependency! 
- CV selector for methods 
- Argument forwarding for Gromacs 
- Updated forward flux examples 
- Significant under-the-hood improvements
- Fixed Gromacs auto-download

## What's Old (v0.6.0)
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

## What's Old (v0.5.0)
- Gromacs restart support
- New gyration tensor CVs
- Updated examples and documentation
- Metadynamics optimizations
- Better engine error handling
- More! (See commit log)

## Known issues 
**SSAGES** is currently in pre-beta. That means there may be known issues that are not yet resolved. Major issues are listed here. 

- Restarts are not fully functioning for all methods.

## Contributing 
Feel free to fork this project on GitHub. Any pull-requests, feature requests or other form of contributions are welcome.
