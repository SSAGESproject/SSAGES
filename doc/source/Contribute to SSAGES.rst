Contribute to SSAGES
====================

The SSAGES project is built on an inclusive and welcoming group of physicits,
chemists, and chemical engineers working on complex Molecular Dynamics
simulations employing Metadynamic techniques. Metadynamics is an exciting and
fast developing field and similarly this project is designed to facilitate the
usage and implementation of a wide array of Metadynamics methods. And we
welcome you heartily to join us and to embark with us on this great adventure.

There are many ways to contribute to SSAGES and you do not necessarily need
programming skills to be part of this project (even though they surely help).
But, if you decide to work on the code base, you will be happy to find that
SSAGES is designed to be easy to use and is just as easy to extend. We put a
high priority on maintaining a readable and clearly structured code base as well
as an inclusive community welcoming new ideas and contributions.

Here is a short summary of ideas how you can become part of SSAGES. Further
down, each point is explained in further detail.

Reporting, Triaging and Fixing Bugs
    No software is without errors, inconsistencies and strange behaviors. Even
    without any programming knowledge you can help tremendously by reporting on
    these bugs or confirming bugs other people have reported.
    :ref:`Read more... <report_bugs>`

Improve the SSAGES documentation
    We highly value a detailed yet comprehensive documentation on what SSAGES
    does and how it does it. This should include concise introductions to the
    methods, quick to learn tutorials, complete coverage of the nooks and
    crannies of each method, and of course helpful pointers in
    case you run into errors. And while the documentation is already expansive,
    improvements on it never go underappreciated here.
    :ref:`Read more... <improve_documentation>`

Include your Method and CV in SSAGES
    You have developed a new Metadyncamics scheme or a Collective Variable and
    want to make it available to the community via SSAGES? Great!
    :ref:`Read more... <add_your_method>`

Work on the core SSAGES system
    If you like to climb down into the SSAGES engine room and get your hands
    dirty this task is for you. :ref:`Read more... <work_on_core_ssages>`

.. _report_bugs:

Report bugs and wishes
----------------------

Link to GitHub issue tracker

.. _improve_documentation:

Improve the Documentation
-------------------------

Every good software is founded on a great documentation. And improvements on the
documentation are always highly appreciated. As you might already have learned,
the SSAGES documentation is split into two parts: The User Manual (which you
are reading right now), and the API documentation. While the Manual uses the
`Sphinx`_ documentation and contains all information necessary to use the
program, the API docs are bulit on `Doxygen`_ and describe the usage of the
underlying classes and functions for everyone willing to extend and improve
SSAGES.

.. _Sphinx: http://sphinx-doc.org
.. _Doxygen: http://www.doxygen.org

Here are a few ideas on how you can help:

* Fix typos: Even though we have throroughly checked, there are certainly still
  a few hidden somewhere.
* Check if all internal and external links are working.
* Make sure that the documentation is up to date, i.e. that it reflects the
  usage of the lastest version.
* Add examples: An examples on how to use a method, avoid a common problem, etc.
  are more helpful than a hundred pages of dry descriptions.
* Write a tutorial.

Build the documentation
^^^^^^^^^^^^^^^^^^^^^^^

Before you can work on the documentation, you first have to build it. The
documentation is part of the SSAGES source code. I assume that you have already
downloaded and built the source code as described in the
:ref:`Getting Started <getting_started>` section. Then, you will find a
collection of rst files comprising the User Manual under ``doc/source/`` where
the file ending ``rst`` stands for ReStructured Text. The API documentation on
the other hand resides directly in the header files right next to the classes
and functions they describe.

Assuming you have already built SSAGES, building the documentation is as easy as
typing
``make doc``
in your build directory. However, before that, make sure that you have the
following programs installed:

* sphinx (with PyPI via ``pip install Sphinx`` for example)
* doxygen
* dot (in Ubuntu this is part of the graphViz package)
* Sphinx "`Read the docs`_" theme (via ``pip install sphinx_rtd_theme``)

.. _Read the docs: https://github.com/snide/sphinx_rtd_theme

Once you have successfully built the documentation you will find the User Manual
under ``doc/Manual/`` and the API documentation under ``doc/API-doc/html/``
(relative to your build directory - do not confuse it with the ``doc/`` folder
in the main directory of the project). Thus, to view it in your favorite web
browser (I use Firefox as an example here) just type

``firefox doc/Manual/index.html``

for the User Manual or

``firefox doc/API-doc/html/index.html``

for the API documentation.

How to write documentation
^^^^^^^^^^^^^^^^^^^^^^^^^^

Here are a few pointers on how to write helpful documentation, before we dive
into the details of **Sphinx** and **Doxygen** for the User Manual and the API
documentation:

* Write documentation "along the way". Do not code first and write the
  documentation later.
* Use helpful error messages. These are considered part of the documentation and
  probably are the part that is read most frequently.
* Do everything you can to structure the text. Let's face it: Most people will
  just skim the documentation. Feel encouraged to use all techniques that
  help to spot the relevant information:

  * Format your text **bold**, *italic*, ``code``, etc.
  * Write in short paragraphs, use headers
  * Use lists, code blocks, tables, etc.

  .. note::

    These Note blocks are extremely helpful for example.

  .. warning::

    Warnings work great, too!

  .. seealso::

    Here you can find more examples for helpful Sphinx markup:
    http://www.sphinx-doc.org/en/stable/markup/para.html

* Use examples, a lot of them
* In the initial stages: Don't be a perfectionist. Missing documentation is the
  worst kind of documentation. Thus, better write more average documentation
  than little documentation worth of a literature prize.

How to write Sphinx
~~~~~~~~~~~~~~~~~~~

The **Sphinx** documentation system uses ReStructured text which is loosely
based on the markdown format. Examples for documentations written with Sphinx
include:

* `LAMMPS`_
* `HOOMD`_
* Virtually all of the `Python`_ Documentation

I found the following tutorials extremely helpful:

* http://www.sphinx-doc.org/en/stable/rest.html
* http://docutils.sourceforge.net/docs/user/rst/quickref.html
* http://openalea.gforge.inria.fr/doc/openalea/doc/_build/html/source/sphinx/rest_syntax.html

.. _LAMMPS: http://lammps.sandia.gov/doc/Manual.html
.. _HOOMD: http://hoomd-blue.readthedocs.io/en/stable/index.html#
.. _Python: https://docs.python.org/3/

One of the great things of Sphinx is that most documentations have a "view page
source" link where you can take a look at the Sphinx source code. Thus, the best
way to learn Sphinx is to click on this link right now and look at the source
code of this page. But here is a short summary of the most important commands:

* Markup: You can use \*italic*, \**bold**, and \``code`` for *italic*, **bold**
  and ``code``.
* Headers. Underline your headers with at least three ``===`` for titles,
  ``---`` for subtitles, ``^^^`` for subsubtitles and ``~~~`` for paragraphs.
* Bullet lists are indicated by lines beginning with ``*``.

.. note::

    These highlighted blocks can be created with ``.. note::``. The content of
    this block needs to be indented. You can also use ``warning`` and
    ``seealso``. Even more can be found
    `here <http://www.sphinx-doc.org/en/stable/markup/para.html>`_.

How to write doxygen
~~~~~~~~~~~~~~~~~~~~

**Doxygen** follows a very different philosphy compared to Sphinx and is more
steered towards API documentation, exactly what we use it for here at SSAGES.
Instead of maintaining the documentation separate from the source code, the
classes and functions are documented in the same place where they are declared:
The header files. Doxygen then reads the source code and autmatically builds
the documentation. Examples for documentation created with doxygen:

* `Plumed`_
* `Root`_

.. _Plumed: http://plumed.github.io/doc-v2.2/user-doc/html/index.html
.. _Root: https://root.cern.ch/doc/master/index.html

The mainpage of the doxygen documentation is written in a separate header file,
in our case ``doc/mainpage.h``. A good introduction to the doxygen syntax can
be found at

* http://www.stack.nl/~dimitri/doxygen/manual/docblocks.html

The basic rule is that doxygen comments start with ``//!`` or ``/*!`` and
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
    This command documents one parameter of the function

\\returns
    This command documents the return value of the function

There are many special doxygen commands. They all start with a backslash and
the most important, apart from the two mentioned above, are:

\\tparam
    Used to document a template parameter.

\\ingroup
    This class is part of a group, such as Methods or Core. The groups are
    defined in ``doc/mainpage.h``.

Helpful are also boxes highlighting a given aspect of the function, such as:

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
    Leave a ToDo note with this command.

You can also highlight your text:

\\em
    For *italic* word. To highlight more text use <em> *Highlighted text* </em>.

\\b
    For **bold** text. To highlight more text use <b> **Bold text** </b>.

\\c
    For typewriter font. To have more text in typewriter font, use
    <tt>Typewriter Font</tt>.

\\code
    Starts a ``code`` block. The block ends with **\\endcode**.

\\li
    A line starting with **\\li** is an entry in a bullet list.

Another big benefit of doxygen is that you can use a lot of LaTeX syntax. For
example:

\\f$
    Starts and ends an inline math equation, similar to $ in Latex.

\\f[ and \\f]
    Start and end a display-style LaTeX equation.

\\cite <label>
    Cite a reference. The references are listed in ``doc/references.bib`` and
    follow the BibTex syntax.

Doxygen also is very clever in producing automatic links. For example, there
exists a class ``Method`` in SSAGES. Thus, doxygen automatically creates a
link to the documentation of this class where the word "Method" appears. This
does, however, not work for the plural, "Methods". Instead, you can write
``\link Method Methods \endlink``. On the other hand, if you want to prevent
doxygen from creating an autolink, put a ``%`` in front of the word.

What to document
^^^^^^^^^^^^^^^^

We are aiming for a comprehensive documentation of all the methods available in
SSAGES as well as the core features. Thus, for each method the documentation
should include

* An introduction into the method, what it does and how it does it.
* A short tutorial based on one of the working examples. The reader should be
  able to complete the tutorial in ~30min and should leave with a sense of
  accomplishment, e.g. a nice energy profile or a picture of a folded protein.
* A detailed description on how to use the method, the parameters, constraints,
  requirements, etc.

.. _add_your_method:

Add your method to SSAGES
-------------------------

.. seealso::

    See :ref:`here <Write-your-own-method>` for an introduction to how to
    develop your own method.

You have developed a new Metadynamics method or a new collective variable (CV)?
Great! We would love to add it to SSAGES. But before we do this, let us check
the following boxes:

* Your code needs to compile and run (obviously).

* If you have implemented a new method, this method should have been published
  in a peer reviewed journal and the publication should be cited in the
  documentation of the method (see next point). If you have implemented a CV,
  please give a small example of usage. In which case does the new CV come in
  handy?

* Your method needs to come with the necessary documentation. For others to
  be able to use your method, you will have to explain how to do it. You can
  take a look at the section :ref:`"How to improve the documentation"
  <improve_documentation>` for a starter on how to write good documentation.

* Please provide an example system. This could be the folding of an
  Alanine Dipeptide molecule, a NaCl system or just a toy model with a simple
  energy landscape. As long as the system is small and the method can easily
  complete within a few hours, it will be fine.

Once these boxes have been checked, our team of friendly code-reviewers will
take a look at your source code and help you meet the high standard of the
SSAGES code.

.. _work_on_core_ssages:

Work on the core classes
------------------------

Describe SSAGES development
