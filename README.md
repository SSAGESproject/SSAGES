SSAGES
============

**SSAGES** (**S**oftware **S**uite for **A**dvanced **G**eneralized **E**nsemble **S**imulations) is an open-source, engine agnostic, C++11 based advanced sampling package. 
It is designed to be easy to use, extendable and extremely versatile. It is currently pre-beta, meaning that there are many rough edges, but we are working rapidly 
to expand its features and fix any bugs. Keep an eye on this page for future updates and see below on how to contribute!

## What's New (v0.4.1)
- String method fix for Gromacs

<a id="features"></a>
## Features
**SSAGES** currently works with LAMMPS and Gromacs molecular dynamics engine. It contains a variety of collective variables (CVs) and advanced sampling methods. 

###Highlights 
- Engine agnostic framework 
- Simple JSON input file syntax 
- Easy to add new CVs 
- Easy to add new methods
- Much more!

### CVs
- Particle coordinate
- Particle position 
- Particle separation 
- Bend angle
- Torsional angle
- Radius of gyration 

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
$ git clone https://github.com/WhitmerGroup/SSAGES-public.git
```
**SSAGES** uses a CMake build system. It also requires either the LAMMPS or Gromacs source codes.
To compile, execute the following

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

###Boost and MPI

SSAGES uses Boost MPI to provide a convenient MPI programming interface. The requirement is for Boost >= 1.58 with the MPI and serialization modules. 
A requisite underlying MPI library is also required. On recent Debian based systems using OpenMPI, the requirement can be installed via:

```bash 
$ sudo apt-get install libopenmpi-dev openmpi-bin libboost-all-dev
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

## Known issues 
**SSAGES** is currently in pre-beta. That means there may be known issues that are not yet resolved. Major issues are listed here. 

- Restart capabilities are currently unimplemented in Gromacs 

## Contributing 
Feel free to fork this project on GitHub. Any pull-requests, feature requests or other form of contributions are welcome.
