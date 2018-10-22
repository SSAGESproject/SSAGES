.. _inputfiles:

Input Files
============

A SSAGES input file contains multiple sections that define the Collective
Variables (CVs), Methods, and various other components that go into an advanced
sampling simulation. There is a brief primer below on JSON, the format used by
SSAGES for input files. The remaining topics describe the basic syntax and
requirements for each section of an input file. Detailed information for
particular methods or collective variables can be found in their respective
locations in this manual.

JSON
----

SSAGES is run using input files written in the JSON_ file format. JSON is a
lightweight text format which is easy for humans to read and write and for
machines to parse and generate. Almost every programming language offers some
level of native JSON support, making it particularly convenient to script or
automate input file generation using, say, Python. If you've never used JSON
before, don't worry. Throughout the documentation we make no assumptions about
the user's knowledge of JSON and provide clear easy-to-follow examples.

A SSAGES input file is a valid JSON document. Here, we will define a bit of
terminology relating to JSON. Take the following JSON structure as an example,
obtained from Wikipedia_:

.. code-block:: javascript

	{
		"firstName": "John",
		"lastName": "Smith",
		"age": 25,
		"address": {
			"streetAddress": "21 2nd Street",
			"city": "New York",
			"state": "NY",
			"postalCode": "10021"
		},
		"phoneNumber": [
			{
				"type": "home",
				"number": "212 555-1234"
			},
			{
				"type": "fax",
				"number": "646 555-4567"
			}
		],
		"gender": {
			"type": "male"
		}
	}

The first pair of curly brackets define the *root* section, we will signify
this using ``#``. An item in the hierarchy, such as the street address, can be
referenced like this:  ``#/address/streetAddress``.

Square brackets ``[]`` in JSON refer to *arrays*, while curly brackets refer
to  *objects*. They can be thought of as Python lists and dictionaries
respectively. That would make ``#/phoneNumber`` an array of phone number
objects, each containing a type and a number. The fax number can be referenced
by ``#/phoneNumber/1/number``, where ``1`` is the array index beginning from
zero.

Items in a JSON object (Python dictionary) are unique. In the example above,
``#/age`` can only be defined once - it is a key in the root tree.  Defining
``#/age`` again will not throw an error, but instead the last definition will
override any previous definitions. This is actually very powerful behavior.
It means a user can import a general template JSON file and override whatever
parameters they wish. The exact behavior of the merging process is described in
detail in the user guide.

Types matter in JSON. Notice how ``#/age`` is specified by a number that is not
surrounded in quotes. This is a number, more specifically an integer. On the
other hand, ``#/address/postalCode`` is a string, even though the contents of
the string are all numbers. Certain fields in a SSAGES input file may be
required to be a string, integer, or number. The user should be aware of this
and take care to format their input file appropriately.

Simulation Properties
---------------------

A SSAGES build is compiled with support for a particular MD engine, and the
requirements for each engine vary slightly. For detailed information on
specific engines and their options check the :ref:`Engines <engines>` section.
The following parameters are needed to define a simulation in the JSON root.

.. warning::

	The properties specified below are case-sensitive. Please be sure to check
	that you have defined it according to the documentation.

Input
~~~~~

The ``"input"`` property specifies the name of the input file used by the
simulation engine.

.. code-block:: javascript

	"input": "in.system"

.. code-block:: javascript

	"input": ["in.system1","in.system2","in.system3"]

The first syntax is used if there is a single input file. For multi-walker
simulations, it is possible to use a single file for all walkers (though this
may not be recommended depending on the method) or specify a separate input
file for each walker.

.. note::

	This property not used by GROMACS (see ``"args"`` property).

Args
~~~~

.. warning::

	This property is *exclusively* for GROMACS and HOOMD-blue.

The ``"args"`` property specifies additional command line arguments to be
passed to the engine.

.. code-block:: javascript

	"args": ["-v", "-deffnm", "runfile"]

.. code-block:: javascript

	"args": "-v -deffnm runfile"

For GROMACS, a standard simulation can be invoked using
``gmx mdrun -deffnm runfile`` to execute a ``runfile.tpr`` binary, the
equivalent arguments must be specified in the ``"args"`` property. This
provides the user with the flexibility of calling command-line arguments in the
same fashion as the standard **mdrun** utility. The only exception is in the
case of multi-walker simulations. If a user wishes to use the multi-walker
capabilities, then ``"args"`` is invoked in the same fashion as a single-walker
simulation. **Do not specify the** ``-multi`` **option. This will be done
automatically.** If ``-deffnm`` is called, GROMACS expects the ``.tpr`` files
for each walker to  be named according to the walker ID starting from zero. In
the example above, if there were three walkers, then GROMACS will look for the
files "runfile0.tpr", "runfile1.tpr", and "runfile2.tpr".

Walkers
~~~~~~~

The ``"walkers"`` property specifies the number of walkers (independent instances
of the simulation engine) to run with SSAGES.

.. code-block:: javascript

	"walkers": 5

Many advanced sampling methods support multi-walker simulations which improve
the convergence of many algorithms. Typically, each walker has an independent
system configuration in a separate input file. It is important to note that
when specifying more than a single walker, the number of processors passed to
``mpiexec`` must be divisible by the number of walkers requested. Otherwise,
SSAGES will terminate with an error.

.. note::

	It is not possible to allocate a different number of processors to each
	walker, at this time.

Collective Variables
~~~~~~~~~~~~~~~~~~~~

The ``"CVs"`` property specifies the collective variables on which SSAGES
will perform its advanced sampling.

.. code-block:: javascript

	"CVs":
	[
		{
			"type": "Torsional",
			"name": "mytorsion_1",
			"atom_ids": [5,7,9,15]
		},
		{
			"type": "ParticleCoordinate",
			"atom_ids": [1],
			"dimension": "x"
		}
	]

Collective variables are specified in an array, where each element is a CV
object. Collective variables can be assigned names or referenced
by index, beginning with zero.

Methods
~~~~~~~

The ``"methods"`` property specifies the advanced sampling algorithms to which
SSAGES will apply to the system.

.. code-block:: javascript

	"methods":
	[
		{
			"type": "Umbrella",
			"ksprings": [100],
			"output_file": "ulog.dat",
			"output_frequency": 10,
			"centers": [1.0],
			"cvs": ["mytorsion_1"]
		},
		{
			"type": "Metadynamics",
			"widths": [0.3],
			"height": 1.0,
			"hill_frequency": 500,
			"lower_bounds": [0.2],
			"upper_bounds": [1.4],
			"lower_bound_restraints": [100],
			"upper_bound_restraints": [100],
			"cvs": [1]
		}
	]

Methods are specified in an array, since it is possible to run multiple methods
simultaneously. This is useful if a user is interested in performing advanced
sampling on a system subject to some restraint, typically applied via an
umbrella. Each method can selectively operate on a subset of CVs by referencing
them either by name or index, as shown above.

Logger
~~~~~~

The ``"logger"`` property specifies an output file to track any or all CVs
as the simulation proceeds.

.. code-block:: javascript

	"logger": {
		"frequency": 100,
		"output_file": "cvs.dat",
		"cvs": [0, 3]
	}

The logger is useful in tracking the evolution of the CVs over the course of an
advanced sampling calculation. Logging CVs can allow for post-simulation
reweighting, or indicate if there are sampling problems in the system being
studied. The frequency of logging the CVs can be specified and each walker in a
multi-walker simulation will have a separate output file. A user can choose to
selectively log individual CVs as well.

Putting It All Together
~~~~~~~~~~~~~~~~~~~~~~~

Combining the previous sections into a single input file yields the following
(purely hypothetical) example input for a LAMMPS simulation.

.. code-block:: javascript

	{
		"walkers": 2,
		"input": ["in.first", "in.second"],
		"CVs":
		[
			{
				"type": "Torsional",
				"name": "mytorsion_1",
				"atom_ids": [5,7,9,15]
			},
			{
				"type": "ParticleCoordinate",
				"atom_ids": [1],
				"dimension": "x"
			}
		],
		"methods":
		[
			{
				"type": "Umbrella",
				"ksprings": [100],
				"output_file": ["walker1.dat", "walker2.dat"],
				"output_frequency": 10,
				"centers": [1.0],
				"cvs": ["mytorsion_1"]
			},
			{
				"type": "Metadynamics",
				"widths": [0.3],
				"height": 1.0,
				"hill_frequency": 500,
				"lower_bounds": [0.2],
				"upper_bounds": [1.4],
				"lower_bound_restraints": [100],
				"upper_bound_restraints": [100],
				"cvs": [1]
			}
		]
	}

To execute this input file, assigning two processors per walker, one would call the command below.

.. code-block:: bash

	mpirun -np 4 ./ssages inputfile.json

.. _JSON: http://json.org
.. _Wikipedia: https://en.wikipedia.org/wiki/JSON#JSON_sample