.. _tutorials:

Tutorials
=========

JSON 
----

SSAGES is run using input files written in the JSON_ file format. JSON is a 
lightweight text format which is easy for humans to read and write and for 
machines to parse and generate. Almost every programming language offers some 
level of native JSON support, making it particularly convenient to script or 
automate input file generation using, say, Python. If you've never used JSON before, 
don't worry. Throughout the documentation we make no assumptions about 
the user's knowledge of JSON and provide clear easy-to-follow examples.

A SSAGES input file is a valid JSON document. Here we will define a bit of 
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


The first pair of curly brackets define the *root* section, we will signify this using 
``#``. An item in the hierarchy, such as the street address, can be referenced like this: 
``#/address/streetAddress``. 

Before moving on to SSAGES specific JSON, we'll mention a few more things about JSON.

- Square brackets ``[]`` in JSON refer to *arrays*, while curly brackets refer to 
  *objects*. They can be thought of as Python lists and dictionaries respectively. 
  That would make ``#/phoneNumber`` an array of phone number objects, each containing 
  a type and a number. The fax number can be referenced by ``#/phoneNumber/1/number``,
  where ``1`` is the array index beginning from zero.

- Items in a JSON object (Python dictionary) are unique. In the example above, 
  ``#/age`` can only be defined once - it is a key in the root tree. 
  Defining ``#/age`` again will not throw an error, but instead the last definition 
  will override any previous definitions. This is actually very powerful behavior. 
  It means a user can import a general template JSON file and override whatever
  parameters they wish. The exact behavior of the merging process is described in 
  detail in the user guide.

- Types matter in JSON. Notice how ``#/age`` is specified by a number that is not
  surrounded in quotes. This is a number, more specifically an integer. 
  ``#/address/postalCode`` on the other hand is a string, even though the contents of 
  the string are all numbers. Some fields in SSAGES may require the input to be a string,
  integer or number and the user should be aware of this difference.


Basic User Tutorial
-------------------

In your SSAGES directory:

.. code-block:: bash

    cd Examples/User/Umbrella

To run a simulation using SSAGES, a ``.json`` file is needed. A ``.json`` file
will tell SSAGES what method it should run, what engine it will use, and what
parameters to use for the method chosen. It will also tell SSAGES how many
simulations to run (walkers) and what engine-input file it should read.  In this
particular example, the engine will be LAMMPS. A butane molecule is used as the
example. All the appropriate LAMMPS input files are provided. The LAMMPS input
files contains the necessary information to perform the simulation. In
``Butane_SSAEGES.in`` you will notice in the last line ``fix ssages all sages``.

A template ``.json`` file (``Template_Input.json``) is provided which contains
the necessary information for a single umbrella simulation. The python code
provided will use the template to generate a new .json file needed for the
simulation. The template ``.json`` file contains the name of the collective
variable (CV) we will use, “Torsional”, and the appropriate atom ids. Run the
python code:

.. code-block:: bash

    python Umbrella_Input_Generator.py 

A new ``.json`` file (``Umbrella.json``) will appear with the correct number of
entries. In this particular example, 12 different walkers are generated. Run
SSAGES:

.. code-block:: bash

    mpiexec -np 12 ./ssages Umbrella.json

where 12 is the number of processors. Since ``Umbrella.json`` contains 12
walkers, 12 processors should be used.

With that, SSAGES will perform Umbrella sampling on a butane molecule biasing
the torsional CV. Output files will be generated for each one of the walkers
containing the iteration number, the target value for the CV, and the CV value
at the iteration number. These values can then be used for further analysis. 

Method-specific tutorials
-------------------------

:ref:`Adaptive Biasing Force <ABF-tutorial>`

:ref:`Basis Function Sampling <BFS-tutorial>`

:ref:`Finite Temperature String <FTS_tutorial>`

:ref:`Forward Flux <FFS_tutorial>`

:ref:`Image Method <IM_tutorial>`

:ref:`Metadynamics <metadynamics-tutorial>`

:ref:`Swarm <Swarm_tutorial>`

:ref:`Umbrella Sampling <Umbrella_tutorial>`

.. _JSON: http://json.org
.. _Wikipedia: https://en.wikipedia.org/wiki/JSON#JSON_sample