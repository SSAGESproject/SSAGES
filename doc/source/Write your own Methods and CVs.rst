.. _Write-your-own-method:

Write your own Methods and CVs
==============================

One of the basic design goals of SSAGES is that it should be easily extensible.
To this end, it provides intuitive and simple tools to implement new collective
variables (CVs) and new metadynamic methods. This section covers the basic steps
to implement a new CV and a new Method. Let us start first with the
implementation of a new CV. The techniques to implement a new Method are
covered :ref:`below <write_new_method>`.

.. _write_new_CV:

How to write a new CV
---------------------

Each CV consists of two components: A header file and a schema file. The header
file contains the source code for the calculation of the CV and the schema file
describes the properties of the CV in a simple JSON format. Finally, you will
have to make SSAGES aware of the new CV.

The CV header file
^^^^^^^^^^^^^^^^^^

Each CV in SSAGES is implemented as a child of the class ``CollectiveVariable``.
The header file should be placed in the directory ``src/CVs`` and has to
(re)implement the following functions:

:code:`void Initialize(const Snapshot&)` (optional)
    This method is called during the pre-simulation phase. It is typically used
    to allocate or reserve memory.

:code:`void Evaluate(const Snapshot&)`
    Evaluation of the CV based on a simulation snapshot. Together with the
    value, this function should also calculate the gradient of the CV. The
    gradient should be a vector of length `n`, where `n` is the number of atoms
    in the Snapshot. Each element in the vector is the derivative of the CV with
    respect to the corresponding atom's coordinates. This method is called
    in the post-integration phase of every iteration.

:code:`double GetValue() const`
    Return the current value of the CV.

:code:`double GetPeriodicValue() const`
    Return the current value of the CV, taking periodic boundary conditions into
    account. An example would be an angular CV which is bound to the region
    :math:`(-\pi,\pi]`. In this case, ``GetValue()`` could return any angle,
    while ``GetPeriodicValue()`` should return the angle mapped back into the
    region :math:`(-\pi,\pi]`. If the CV does not use periodic boundaries, this
    function should return the same value as ``GetValue()``.

:code:`const std::vector<Vector33>& GetGradient() const`
    Return the gradient of the CV (see ``Evaluate(const Snapshot&)`` for how the
    gradient is defined).

:code:`const std::array<double, 2>& GetBoundaries() const`
    Return a two-element array containing the lower and the upper boundary for
    the CV.

:code:`double GetDifference(const double Location) const`
    Return the distance of the current value of the CV from a specified
    location, taking periodic boundary conditions into account. If the CV does
    not use periodic boundary conditions, the return value should simply be
    ``GetValue() - Location``.

The CV schema file
^^^^^^^^^^^^^^^^^^

Together with the header file that contains the source code of the CV, you will
have to provide a schema file to make the CV accessible to the SSAGES input
files. The schema file should be placed in the directory ``schema/CVs/``. It
has to be written in the JSON format and should contain the following items:

type
    The value of *type* should be set to ``object``.

varname
    The name of your new CV.

properties
    The properties contain the *type* which is the internal name of the CV and
    a set of other properties that have to be supplied to the constructor of the
    CV.

required
    A list containing the required properties. Optional parameters to the CV
    constructor are not listed here.

additionalProperties
    Optional properties.

Integrate the new CV into SSAGES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have provided the header and the schema file, there is one more
steps to do in order to make SSAGES aware of the newly included CV.

.. note::

    We are currently working on a method to automate this step. Revisit this
    section in future releases. Chances are, that you no longer have to worry
    about this step.

To include your new CV, you have to edit the file
``src/CVs/CollectiveVariable.cpp``, and

1. ``#include`` your CV header file at the top of the file.
2. Add a new ``else if`` clause in ``BuildCV()``. The if-test checks for the
   CV type set as an enum in the list of properties. Within the if-clause you
   should parse and validate the JSON schema, read the required properties and
   create the CollectiveVariable. A pointer to the newly created object should
   be stored in the variable named ``cv``.

.. _write_new_method:

How to write a new method
-------------------------

Each method consists of three components: A header file, a cpp file, and a schema file. The header
file and cpp file contains the source code for the method and the schema file
describes the properties of the method in a simple JSON format. Finally, you will
have to make SSAGES aware of the new method.

The method header file
^^^^^^^^^^^^^^^^^^^^^^

Each method in SSAGES is implemented as a child of the class ``Methods``.
The header file should be placed in the directory ``src/methods`` and has to
(re)implement the following functions:

:code:`void PreSimulation(Snapshot* snapshot, const CVList& cvs)`
    Setup done before the method actually runs. This function will be called
    vefore the simulation is started.

:code:`void PostIntegration(Snapshot* snapshot, const CVList& cvs)`
    This is where the heart of your method should go. By using snapshot and 
    the cvs, modify the forces, positions, velocities, etc. appropriated by 
    the new method. This function will be called after each integration MD step.

:code: `void PostSimulation(Snapshot* snapshot, const CVList& cvs)`
    This function is called at the end of the simulation run. Use it to close files
    your method opened, to write out data that you have been storing, etc.

The method schema file
^^^^^^^^^^^^^^^^^^^^^^

Together with the source code of the method, you will
have to provide a schema file to make the CV accessible to the SSAGES input
files. The schema file should be placed in the directory ``schema/methods/``. It
has to be written in the JSON format and should contain the following items:

type
    The value of *type* should be set to ``object``.

varname
    The name of your new method.

properties
    The properties contain the *type* which is the internal name of the method and
    a set of other properties that have to be supplied to the constructor of the
    method.

required
    A list containing the required properties. Optional parameters to the method
    constructor are not listed here.

additionalProperties
    Optional properties.

Integrate the new method into SSAGES
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Once you have provided the header and the schema file, there is one more
steps to do in order to make SSAGES aware of the newly included method.

.. note::

    We are currently working on a method to automate this step. Revisit this
    section in future releases. Chances are, that you no longer have to worry
    about this step.

To include your new method, you have to edit the file
``src/methods/Methods.cpp``, and

1. ``#include`` your method header file at the top of the file.
2. Add a new ``else if`` clause in ``BuildMethod()``. The if-test checks for the
   method type set as an enum in the list of properties. Within the if-clause you
   should parse and validate the JSON schema, read the required properties and
   create the method. A pointer to the newly created object should
   be stored in the variable named ``method``.
