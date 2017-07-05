SSAGES
============

**SSAGES** (**S**oftware **S**uite for **A**dvanced **G**eneralized **E**nsemble **S**imulations) is an open-source, engine agnostic, C++11 based advanced sampling package. 
It is designed to be easy to use, extendable and extremely versatile. It is currently pre-beta, meaning that there are many rough edges, but we are working rapidly 
to expand its features and fix any bugs. Keep an eye on this page for future updates and see below on how to contribute!

## What's New (v0.7.0)
- New simplified JSON syntax
- Support for multiple simultaneous methods! 
- Eliminated boost dependency! 
- CV selector for methods 
- Argument forwarding for Gromacs 
- Updated forward flux examples 
- Significant under-the-hood improvements
- Fixed Gromacs auto-download

<a id="features"></a>
## Features
**SSAGES** currently works with LAMMPS and Gromacs molecular dynamics engines. It contains a variety of collective variables (CVs) and advanced sampling methods. 

### Highlights 
- Engine agnostic framework 
- Simple JSON input file syntax 
- Easy to add new CVs 
- Easy to add new methods
- Much more!

### Engines 
- Gromacs 5.x.x
- LAMMPS (Most recent versions)
- OpenMD (2.4+)
- QBox (1.63+)

### CVs
- Atom group coordinate
- Atom group position 
- Atom group separation 
- Bend angle
- Torsional angle
- Components of gyration tensor
- Polymer Rouse modes 
- Coordination number
- Box volume

### Methods 
- Metadynamics 
- Basis function sampling 
- Adaptive biasing force 
- Nudged elastic band 
- Finite temperature string 
- Forward flux sampling 
- Swarm of trajectories 
- Umbrella sampling 

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

## What's Old (v0.4.3)
- Added time dependent umbrella centers (steered MD)
- Fixed biases resizing bug in ABF method

## Known issues 
**SSAGES** is currently in pre-beta. That means there may be known issues that are not yet resolved. Major issues are listed here. 

- Examples have not been completely updated to reflect the new JSON syntax


## Contributing 
Feel free to fork this project on GitHub. Any pull-requests, feature requests or other form of contributions are welcome.
