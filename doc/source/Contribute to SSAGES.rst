Contribute to SSAGES
====================

The SSAGES project is built by an inclusive and welcoming group of physicists,
chemists, and chemical engineers working on complex Molecular Dynamics (MD)
simulations employing advanced sampling techniques. Advanced sampling is an
exciting and rapidly developing field. Similarly, this project is designed to
facilitate the usage and implementation of a wide array of sampling methods.
We welcome you heartily to join us as we embark on this great adventure.

There are many ways to contribute to SSAGES, and you do not necessarily need
programming skills to be part of this project (even though they surely help).
But, if you decide to work on the code base, you will be happy to find that
SSAGES is designed to be easy to use and is just as easy to extend. We put a
high priority on maintaining a readable and clearly structured code base as
well as an inclusive community welcoming new ideas and contributions.

Here is a short summary of ideas how you can become part of SSAGES:

Reporting, Triaging, and Fixing Bugs
    No software is without errors, inconsistencies, and strange behaviors. Even
    with zero programming knowledge, you can help tremendously by reporting
    bugs or confirming issued bugs.
    :ref:`Read more... <report_bugs>`

Improving the SSAGES documentation
    We strive for SSAGES to have a detailed, yet comprehensive, documentation
    on what it does and how it does it. This should include concise
    introductions to methods, quick to learn tutorials, complete coverage of
    the intricacies of each method, and helpful pointers in case you run into
    errors. While the documentation is already expansive, improvements on it
    never go unappreciated.
    :ref:`Read more... <improve_documentation>`

Including your Method and CV in SSAGES
    You have developed a new sampling method or collective variable and
    want to make it available to the community via SSAGES? Great!
    :ref:`Read more... <add_your_method>`

Working on the core SSAGES system
    If you would like to climb into the heart of SSAGES and get your hands
    dirty, this task is for you! :ref:`Read more... <work_on_core_ssages>`

.. _report_bugs:

Reporting Bugs and Requesting Features
--------------------------------------

SSAGES is an open-source project maintained by a small team of scientific
developers. Because of this, we appreciate any help that the community can
provide!

If you have been using SSAGES and have come across an error that might affect
other users or the core code base, please report it on the `GitHub Issues`_
page within the repository. When reporting bugs, it is helpful to include as
much detail as possible. This can include, but is not limited to, input files
necessary for replicating the bug, the output you expected to see versus what
you got from SSAGES, and any relevant error messages. Before submitting a bug,
browse the current issues to make sure that it has not already been reported.

If you feel that there is a feature missing from SSAGES, you can either submit
a feature request (also via `GitHub Issues`_) or you can add the feature
yourself and submit a pull request on the `GitHub Pull Requests`_ page. When
submitting a feature request, please include a reference to relevant literature
for new methods/CVs/etc., as well as a generalized use case (or more) that can
motivate adding in the requested feature.

.. _GitHub Issues: https://github.com/MICCoM/SSAGES-public/issues
.. _GitHub Pull Requests: https://github.com/MICCoM/SSAGES-public/pulls

.. _improve_documentation:

Improving the Documentation
---------------------------

| Great documentation and great code produces great software.
| ---SSAGE advice
|

Improvements on the documentation are always highly appreciated.
The SSAGES documentation is split into two parts: The User Manual (which you
are reading right now) and the API documentation. While the User Manual uses
the `Sphinx`_ documentation and contains all information necessary to use the
program, the API docs are built on `Doxygen`_ to describe the usage of the
underlying classes and functions for everyone willing to extend and improve
SSAGES.

.. _Sphinx: http://sphinx-doc.org
.. _Doxygen: http://www.doxygen.org

Here are a few ideas on how you can help:

* Fix typos: Even though we are regularly checking, there are certainly still
  a few hidden somewhere.
* Check links: Verify all internal and external links are working.
* Bring up to date: Make sure that the documentation is current, i.e. that it
  reflects the usage of the latest version.
* Add examples: An example on how to use a method, avoid a common problem, etc.
  is more helpful than a hundred pages of dry descriptions.
* Write a tutorial: If there are aspects of the code that are missing a
  tutorial, and one might be useful, consider adding a short tutorial with a
  simple toy system.

Building the Documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Before you can work on the documentation, you first have to build it. The
documentation is part of the SSAGES source code. It is assumed that you have
already downloaded and built the source code as described in the
:ref:`Getting Started <getting-started>` section. You will find a
collection of ``.rst`` files comprising the User Manual under ``doc/source/``
where the file ending ``.rst`` stands for ReStructured Text. The API
documentation, on the other hand, resides directly in the header files, right
next to the classes and functions they describe.

Assuming you have already built SSAGES, building the documentation is as easy as
typing

.. code:: bash

	``make doc``

within your build directory. In order to make the documentation, you must have
the following programs installed:

* Sphinx (with PyPI via ``pip install Sphinx`` for example)
* Doxygen
* dot (in Ubuntu this is part of the graphViz package)
* Sphinx "`Read the docs`_" theme (via ``pip install sphinx_rtd_theme``)

.. _Read the docs: https://github.com/snide/sphinx_rtd_theme

Once you have successfully built the documentation, you will find the User Manual
under ``doc/Manual/`` and the API documentation under ``doc/API-doc/html/``
(relative to your build directory - do not confuse it with the ``doc/`` folder
in the main directory of the project). To view it in your favorite web
browser (using FireFox as an example) just type

``firefox doc/Manual/index.html``

for the User Manual or

``firefox doc/API-doc/html/index.html``

for the API documentation.

Writing Documentation
^^^^^^^^^^^^^^^^^^^^^

Here are a few pointers on how to write helpful documentation, before we dive
into the details of **Sphinx** and **Doxygen** for the User Manual and the API
documentation:

* Write documentation "along the way". Do not code first and write the
  documentation later.
* Use helpful error messages. These are considered part of the documentation and
  probably are the part that is read most frequently.
* Do everything you can to structure the text. Let's face it: most people will
  just skim the documentation. Feel encouraged to use any and all techniques that
  help to spot the relevant information, for example:

  * Format your text **bold**, *italic*, ``code``, etc.
  * Use effective headers
  * Write in short paragraphs
  * Use lists, code blocks, tables, etc.

  .. note::

    These Note blocks are extremely helpful for example.

  .. warning::

    Warnings work great, too!

  .. seealso::

    Here you can find more examples for helpful Sphinx markup:
    
    http://www.sphinx-doc.org/en/stable/markup/para.html

* Use examples, a lot of them.
* In the initial stages: Don't be a perfectionist. Missing documentation is
  the worst kind of documentation.

| "It is better to have written and coded than to have never written at all."
| ---SSAGE advice
|

Documenting with Sphinx
~~~~~~~~~~~~~~~~~~~~~~~

The **Sphinx** documentation system uses reStructuredText which is loosely
based on the Markdown format. Examples for documentations written with Sphinx
include:

* `LAMMPS`_
* `HOOMD-blue`_
* Virtually all of the `Python`_ Documentation

The following tutorials are extremely helpful:

* http://www.sphinx-doc.org/en/stable/rest.html
* http://docutils.sourceforge.net/docs/user/rst/quickref.html
* http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html

.. _LAMMPS: http://lammps.sandia.gov/doc/Manual.html
.. _HOOMD-blue: http://hoomd-blue.readthedocs.io/en/stable/index.html#
.. _Python: https://docs.python.org/3/

One of the great things of Sphinx is that most documentations have a "view page
source" link at the top of the page, where you can take a look at the Sphinx
source code. Thus, the best way to learn Sphinx is to click on this link right
now and look at the source code of this page. But here is a short summary of
the most important commands:

* Markup: You can use \*italic*, \**bold**, and \``code`` for *italic*, **bold**
  and ``code``.
* Headers: Underline your headers with at least three ``===`` for titles,
  ``---`` for subtitles, ``^^^`` for subsubtitles and ``~~~`` for paragraphs.
* Lists: Bulleted lists are indicated by lines beginning with ``*``.

.. note::

    These highlighted blocks can be created with ``.. note::``. The content of
    this block needs to be indented. You can also use ``warning`` and
    ``seealso``. Even more can be found
    `here <http://www.sphinx-doc.org/en/stable/markup/para.html>`_.

Documenting with Doxygen
~~~~~~~~~~~~~~~~~~~~~~~~

**Doxygen** follows a very different philosophy compared to Sphinx and is more
steered towards API documentation, exactly what we use it for in SSAGES.
Instead of maintaining the documentation separately from the source code, the
classes and functions are documented in the same place where they are declared:
the header files. Doxygen then reads the source code and automatically builds
the documentation. Examples for documentation created with Doxygen include:

* `PLUMED`_
* `Root`_

.. _PLUMED: https://plumed.github.io/doc-v2.3/user-doc/html/index.html
.. _Root: https://root.cern.ch/doc/master/index.html

The mainpage of the Doxygen documentation is written in a separate header file,
in our case ``doc/mainpage.h``. A good introduction to the Doxygen syntax can
be found at

* http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html

The basic rule is that Doxygen comments start with ``//!`` or ``/*!`` and
document the class, namespace or function that directly follows it. Let's start
with a short example:

.. code-block:: cpp

    //! Function taking the square of a value
    /*!
     * \param val Input value
     * \returns Square of the input value
     *
     * This function calculates the square of a given value.
     */
    double square(double val)
    {
        return val*val;
    }

This example documents the function ``square()`` which simply calculates the
square of a number. The first line, starting with ``//!``, is the brief
description and should not be longer than one line. The second comment block,
starting with ``/*!`` is the full description. Here, two special commands
are used:

\\param
    This command documents one parameter of the function.

\\returns
    This command documents the return value of the function.

There are many special Doxygen commands. They all start with a backslash and
the most important, apart from the two mentioned above, are:

\\tparam
    Used to document a template parameter.

\\ingroup
    This class is part of a group, such as Methods or Core. The groups are
    defined in ``doc/mainpage.h``.

There are also helpful boxes that highlight a given aspect of the function, such as:

\\attention
    Puts the following text in a raised box. A blank line ends the attention box.

\\note
    Starts a highlighted block. A blank line ends the note block.

\\remark
    Starts a paragraph where remarks may be entered.

\\see
    Paragraph for "See also".

\\deprecated
    The documented class or function is deprecated and only kept for backwards
    compatibility.

\\todo
    Leave a "to do" note with this command.

You can also highlight your text:

\\em
    For an *italic* word. To apply to multiple words, use <em> *Italic text* </em>.

\\b
    For a **bold** word. To apply to multiple words, use <b> **Bold text** </b>.

\\c
    For a ``code`` word (in typewriter font). To apply to multiple words, use
    <tt> ``Typewriter`` </tt>.

\\code
    Starts a ``code`` block. The block ends with **\\endcode**.

\\li
    A line starting with **\\li** is an entry in a bullet list.

Another benefit of Doxygen is that you can use some LaTeX syntax. For example:

\\f$
    Starts and ends an inline math equation, similar to $ in Latex.

\\f[ and \\f]
    Start and end a display-style LaTeX equation.

\\cite <label>
    Cite a reference. The references are listed in ``doc/references.bib`` and
    follow BibTex syntax.

Doxygen is very clever in producing automatic links. For example, there
exists a class ``Method`` in SSAGES. Thus, Doxygen automatically creates a
link to the documentation of this class where the word "Method" appears. This,
however, does not work for the plural, "Methods". Instead, you can write
``\link Method Methods \endlink``. On the other hand, if you want to prevent
Doxygen from creating an autolink, put a ``%`` in front of the word.

What to Document
^^^^^^^^^^^^^^^^

We strive for a comprehensive documentation of all the methods available in
SSAGES, as well as the core features. Thus, for each method, the documentation
should include:

* An introduction to the method, what it does, and how it does it.
* A short tutorial based on one of the working examples. The reader should be
  able to complete the tutorial in ~30 min and should leave with a sense of
  accomplishment, e.g. a nice energy profile or a picture of a folded protein.
* A detailed description on how to use the method, the parameters, constraints,
  requirements, etc.

.. _add_your_method:

Adding your Method to SSAGES
----------------------------

.. seealso::

    See :ref:`here <Write-your-own-method>` for an introduction to how to
    develop your own method.

So, you have developed a new sampling method or new collective variable (CV)?
Great! SSAGES is about collaboration and integrating your new CV or method is
a priority. But before we do that, make sure you check the following boxes:

* Your code needs to compile and run (obviously).
* If you have implemented a new method, this method should have been published
  in a peer reviewed journal and the publication should be cited in the
  documentation of the method (see next point). If you have implemented a CV,
  please give a small example of usage. In which case(s) does the new CV come
  in handy?
* Your method needs to come with the necessary documentation. For others to
  be able to use your method, you will have to explain how it works. You can
  take a look at the section :ref:`"Improving the Documentation"
  <improve_documentation>` for a starter on how to write good documentation.
* Please provide an example system. This could be the folding of an
  Alanine Dipeptide molecule, a NaCl system, or just a toy model with a simple
  energy landscape. As long as the system is small and the method can easily
  complete within a few hours, it will be fine.

Once these boxes have been checked, our team of friendly developers will review
your source code to help you meet the standards of the SSAGES code.

.. _work_on_core_ssages:

Working on the Core Classes
---------------------------

.. todo::

    Describe SSAGES code hierarchy and development workflow.
