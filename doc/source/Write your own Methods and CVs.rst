.. _Write-your-own-method:

Write Your Own CVs and Methods
==============================

One of the basic design goals of SSAGES is that it should be easily extensible.
To this end, it provides intuitive and simple tools to implement new collective
variables (CVs) and new advanced sampling methods. This section covers the
basic steps to implement a new CV or a new Method. Let us first start with the
implementation of a new CV. The techniques to implement a new Method are
covered :ref:`below <write_new_method>`.

.. _write_new_CV:

How to Write a New CV
---------------------

Each CV consists of two components: a header file and a schema file. The header
file contains the source code for the calculation of the CV and the schema file
describes the properties of the CV in a simple JSON format. Finally, you will
have to make SSAGES aware of the new CV by incorporating it into the core
classes.

The CV Header File
^^^^^^^^^^^^^^^^^^

Each CV in SSAGES is implemented as a child of the class ``CollectiveVariable``.
The header file should be placed in the directory ``src/CVs/`` and has to
(re)implement the following functions:

:code:`void Initialize(const Snapshot&)` (optional)
    Called during the pre-simulation phase. It is typically used
    to allocate or reserve memory.

:code:`void Evaluate(const Snapshot&)`
    Evaluation of the CV based on a simulation snapshot. This function should
    calculate both the value and gradient of the CV. The gradient should be a
    vector of length :math:`n`, where :math:`n` is the number of atoms
    in the Snapshot. Each element in the vector is the derivative of the CV with
    respect to the corresponding atom's coordinates. This method is called
    in the post-integration phase of every iteration.

:code:`double GetValue() const`
    Return the current value of the CV, as calculated in
    ``Evaluate(const Snapshot&)``.

:code:`double GetPeriodicValue() const`
    Return the current value of the CV, taking periodic boundary conditions into
    account. An example would be an angular CV which is bound to the region
    :math:`(-\pi,\pi]`. In this case, ``GetValue()`` could return *any* angle,
    while ``GetPeriodicValue()`` should return the angle mapped back into the
    region :math:`(-\pi,\pi]`. If the CV does not use periodic boundaries, this
    function should return the same value as ``GetValue()``.

:code:`const std::vector<Vector3>& GetGradient() const`
    Return the gradient of the CV (see ``Evaluate(const Snapshot&)`` for how the
    gradient is defined).

:code:`const std::array<double, 2>& GetBoundaries() const`
    Return a two-element array containing the lower and the upper boundary for
    the CV.

:code:`double GetDifference(double Location) const`
    Return the distance of the current value of the CV from a specified
    location, taking periodic boundary conditions into account. If the CV does
    not use periodic boundary conditions, the return value should simply be
    ``GetValue() - Location``.

The CV Schema File
^^^^^^^^^^^^^^^^^^

Together with the header file that contains the source code of the CV, you will
have to provide a schema file to make the CV accessible to the SSAGES input
files. The schema file should be placed in the directory ``schema/CVs/``. It
has to be written in the JSON format and should contain the following items:

type
    The value of ``type`` should be set to ``object``.

varname
    The name of your new CV, of the form ``ExampleCV``.

properties
    The properties contain the ``type`` which is the internal name of the CV and
    a set of other properties that have to be supplied to the constructor of the
    CV.

required
    A list containing the required properties. Optional parameters to the CV
    constructor are not listed here.

additionalProperties
    A boolean value allowing unlisted JSON members to be included in input files.

Integrate the New CV into SSAGES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have provided the header and schema files, there are two more
steps in order to make SSAGES aware of the newly included CV.

To include your new CV, you have to edit the file
``src/CVs/CollectiveVariable.cpp``, and

1. ``#include`` your CV header file at the top of the file.
2. Add a new ``else if`` clause in ``BuildCV()``. The if-test checks for the
   type of CV selected and returns the respective ``Build(json, path)``
   function from the new CV.

.. _write_new_method:

How to Write a New Method
-------------------------

Each method consists of three components: a header file, a cpp file, and a
schema file. The header and cpp files contain the source code for the method
and the schema file describes the properties of the method in a simple JSON
format. Finally, you will have to make SSAGES aware of the new method by
incorporating it into the core classes.

The Method Header File
^^^^^^^^^^^^^^^^^^^^^^

Each method in SSAGES is implemented as a child of the class ``Methods``.
The header file should be placed in the directory ``src/methods`` and has to
(re)implement the following functions:

:code:`void PreSimulation(Snapshot* snapshot, const CVList& cvs)`
    Called before the method actually runs. It is typically used
    to allocate or reserve memory.

:code:`void PostIntegration(Snapshot* snapshot, const CVList& cvs)`
    Called after each MD integration step. This is where the heart of your
    method should go. By using Snapshot and the CVs, this function modifies
    the forces, positions, velocities, etc. appropriated by the new method.

:code:`void PostSimulation(Snapshot* snapshot, const CVList& cvs)`
    Called at the end of the simulation run. Use it to close files
    your method opened, to write out data that the method is storing, etc.

The Method Schema File
^^^^^^^^^^^^^^^^^^^^^^

Together with the source code of the method, you will
have to provide a schema file to make the method accessible to the SSAGES input
files. The schema file should be placed in the directory ``schema/Methods/``. It
has to be written in the JSON format and should contain the following items:

type
    The value of ``type`` should be set to ``object``.

varname
    The name of your new method, of the form ``ExampleMethod``.

properties
    The properties contain the ``type`` which is the internal name of the
    method and a set of other properties that have to be supplied to the
    constructor of the method.

required
    A list containing the required properties. Optional parameters to the method
    constructor are not listed here.

additionalProperties
    A boolean value allowing unlisted JSON members to be included in input files.

Integrate the New Method into SSAGES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have provided the header, source, and schema files, there are two more
steps in order to make SSAGES aware of the newly included method.

To include your new method, you have to edit the file
``src/methods/Method.cpp``, and

1. ``#include`` your method header file at the top of the file.
2. Add a new ``else if`` clause in ``BuildMethod()``. The if-test checks for
   the type of method selected and calls the respective
   ``Build(json, world, comm, path)`` function from the new method. A pointer
   to the newly created object should be stored in the variable named ``method``.

Finally, add the method ``.cpp`` file to CMakeLists.txt as a source.
