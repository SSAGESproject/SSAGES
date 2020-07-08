SSAGES Documentation
====================

SSAGES's documentation can be found in the [online
manual](https://ssagesproject.github.io/docs/index.html). Alternatively, after you
have built SSAGES you can further build the documentation from the `build/`
directory.

## Requirements

Make sure you have the following additional tools installed:

- [Doxygen]              – generate documentation from annotated C++ source code.
- [Graphviz]             – `dot` tool used to draw graph visualizations.
- [Sphinx]               – documentation builder.
- [sphinx-rtd-theme]     – Read the Docs Sphinx theme.
- [sphinxcontrib-bibtex] – Sphinx extensions for BibTeX style citations.

[Doxygen]:              https://www.doxygen.nl/index.html
[Graphviz]:             https://graphviz.org
[Sphinx]:               https://sphinx-doc.org
[sphinx-rtd-theme]:     https://sphinx-rtd-theme.readthedocs.io
[sphinxcontrib-bibtex]: https://sphinxcontrib-bibtex.readthedocs.io

On Debian-based distributions (e.g. Ubuntu), you can easily install them with
the following commands:
```
$ sudo apt install doxygen graphviz python3-sphinx python3-pip
$ pip install sphinx_rtd_theme sphinxcontrib-bibtex
```

## Building

You can build the documentation with
```
$ make doc
```
You can also build the API-references and the User Manual separately with
```
$ make apiref
```
and
```
$ make manual
```

If you have `pdflatex` installed, you can also build a PDF file for the
documentation. To compile the API-reference into a PDF file do
```
$ cd doc/API-doc/latex/
$ make
```
The PDF will be called `refman.pdf`

Similarly, you can build a PDF version of the Manual with
```
$ cd doc/Manual/
$ make
```
The PDF will be called `SSAGES.pdf`

## Viewing the documentation

Once you have built the documentation you will find it in the `doc/API-doc/`
and `doc/Manual` directories. To view the documentation in a browser just open
the files `doc/Manual/index.html` or `doc/API-doc/html/index.html`.
