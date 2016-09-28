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

This will build a MPI compatible LAMMPS binary which will reside in the LAMMPS source
directory.

If you want to compile and run unit and integration tests, replace the cmake command
in the example above with

```bash
$ cmake -DLAMMPS_SRC=/path/to/lammps/src -DBUILD_TESTS=ON ..
```

The build system is still a work in progress.

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
