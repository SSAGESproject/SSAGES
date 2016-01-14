SSAGES
============

<a id="installation"></a>
## Installation
The first step is to clone the repository locally.

```bash
$ git clone https://github.com/jonathankwhitmer/SSAGES.git
```
**SSAGES** uses a CMake build system. It also requires the LAMMPS source code.
To compile, execute the following

```bash
$ cd SSAGES
$ mkdir build && cd build
$ cmake -DLAMMPS_SRC=/path/to/lammps/src .. 
$ make
```

This will build and MPI compatible LAMMPS binary which will reside in the LAMMPS source 
directory. 

The build system will is a work in progress. 
