.. _artificial-neural-network-sampling: 

Artificial Neural Network Sampling
----------------------------------

Introduction 
^^^^^^^^^^^^

Artificial neural network (ANN) sampling is a free energy sampling method which uses 
ANNs to generate an on-the-fly adaptive bias capable of rapidly resolving free energy 
landscapes. It is a recent method proposed by Sidky & Whitmer which is demonstrated 
to be robust to user inputs and requires a minimal number of parameters. It is quite 
simple to use as described below. 

Like Basis Function Sampling, the algorithm proceeds in sweeps. Statistics are collected
over a user specified interval, which is a very flexible choice. At the end of each sweep 
an ANN is fit to the statistics to determine the optimal bias which is then applied in the 
subsequent sweep. This proceeeds until the free energy landscape has converged. 

Detailed information on the method and implementation can be found in
the publication [1]_. 

Options & Parameters
^^^^^^^^^^^^^^^^^^^^

ANN sampling is selected by defining ``"type" : "ANN"`` as the 
method in the JSON input file. It supports the following options:

topology
	*Array of integers (length: number of hidden layers)*. 
	This array defines the architecture of the neural network. ANN sampling 
	is quite robust to the choice of network architure. Nonetheless, there are 
	very good heuristics that will ensure you have a network that trains quickly
	and performs well. The key rule is to use a network with a single hidden layer 
	if you are biasing on one CV, and two hidden layers if you are biasing on two or 
	more. **Please take a look at the example section for more details on network
	architectures**.

grid 
	This is a grid object which defines how CV space will be discretized. 

nsweep 
	*integer* 
	Defines the length of a sweep in iterations. Typical values range from 1,000 to 10,000 depending 
	on the size of the system. The slower the system dynamics, the longer the sweep. This is not going to heavily 
	affect convergence time, and the method is generally quite robust to choice of sweep length. The main 
	consequence of this choice is that the ANN training time become relatively expensive if the system is 
	very cheap to evaluate. 

weight
	*double* 
	(Default = 1) 
	This defines how much relative weight is assigned to the statistics collected during a sweep. 
	The default value works fine in all cases. However, if you know the free energy barriers are very large 
	or very small, you can adjust this to speed up convergence.

temperature
	*double* 
	The temperature of the simulation must be specified. 

output_file
	*string* 
	(Default = ann.out) 
	This specifies the name of the output file containing the bias and reweighted histogram. 

overwrite_output
	*boolean* 
	(Default = true)
	Setting this to false will append the sweep number to the output file name so that a series of files 
	will be created over time. This is useful if a user is interested in plotting the evolution of the bias 
	over time. 

lower_bounds
    *array of doubles (nr of CVs) long*.
    This array defines the minimum values for the CV restraints.

upper_bounds
    *array of doubles (nr of CVs) long*.
    This array defines the maximum values for the CV restraints.

lower_bound_restraints
    *array of doubles (nr of CVs) long*.
    This array defines the spring constants for the lower bound CV restraints. If you do not 
    wish to restrain the CVs to a particular interval, you can set this to zero.

upper_bound_restraints
    *array of doubles (nr of CVs) long*.
    This array defines the spring constants for the upper bound CV restraints. If you do not 
    wish to restrain the CVs to a particular interval, you can set this to zero.

max_iters
	*integer* 
	(Default = 1000)
	Defines the maximum number of training iterations/epochs per sweep. If you are running a very 
	large network or a large number of CVs then this can be made quite small (10 to 100) since the 
	network will continually improve as the system proceeds. Unless you have a very large network, 
	or a very cheap system, you can leave this alone.

prev_weight
	*double*
	(Default = 1) 
	This is a special feature that allows ANN sampling to "forget" previous history specified by a fraction.
	For example, if you want to retain 90% of accumulated data each sweep, set this to 0.9. This is useful if 
	you are sampling along a CV you know to be bad, and there is clear quasi-nonergodicity in the system. By 
	allowing the ANN method to forget accumulated data, you allow it to adapt more efficiently to newly accessible 
	regions of phase space. Note that unless you really know why you need to be using this, you should leave it 
	alone. 

Example Input 
^^^^^^^^^^^^^

There is an online repository containing complete examples and analysis for systems of up to 4 CVs. 
Along with the original publication, this should provide you a good idea of how to specify the network 
architecture. The repository can be found `here <https://github.com/hsidky/ann_sampling>`__. 

For a 1D CV such as Na-Cl distance, the input is as follows 

.. code-block:: javascript 

	"methods" : [
		{
			"type" : "ANN",
			"topology" : [15],
			"nsweep" : 10000,
			"temperature" : 298.15,
			"grid" : {
			"lower" : [0.24],
			"upper" : [0.80],
			"number_points" : [300],
			"periodic": [false]
			},
			"lower_bounds" : [0.23],
			"upper_bounds" : [0.81],
			"lower_bound_restraints" : [0],
			"upper_bound_restraints" : [1000]
		}
	]

for a 2D CV such as the dihedral angles of ADP, the input is as follows

.. code-block:: javascript

	"methods" : [
		{
			"type" : "ANN", 
			"topology" : [10, 6],
			"nsweep" : 5000, 
			"overwrite_output" : false,
			"temperature" : 298.15,
			"grid" : {
				"lower" : [-3.141592653589793, -3.141592653589793],
				"upper" : [3.141592653589793, 3.141592653589793], 
				"number_points" : [30, 30],
				"periodic" : [true, true]
			},
			"lower_bounds" : [-4, -4],
			"upper_bounds" : [4, 4],
			"lower_bound_restraints" : [0, 0],
			"upper_bound_restraints" : [0, 0]
		}
	]

For more examples, and higher dimensions, please check out the repository linked above. 

Output
^^^^^^

ANN sampling writes either a single output file or a series of output files over time. Each 
file contains columns corresponding to the CVs, a column containing the unbiased histogram 
estimate and a final column containing the bias. The format is as follows: 

* cv1 cv2 ... histogram bias * 

This file can be loaded and visualized easily in many scripting languages such as Python and 
MATLAB. An exmaple of how to load data in Python for a 2D CV is shown below.

.. code-block:: python 

    # Load data.
    X = np.loadtxt("ann.dat")
    xg = np.reshape(X[:,0], (61, 61))
    yg = np.reshape(X[:,1], (61, 61))
    zg = np.reshape(-X[:,3], (61, 61))
    zg = zg - np.max(zg)
    
    # Plot data.
    fig = plt.figure(figsize=(5,5))
    plt.contour(xg, yg, zg, linewidths=0.5, colors="k")
    plt.contourf(xg, yg, zg)

A file called "netstate.dat" is also written out which contains the neural network parameters. 
This network can be evaluated in Python using a ANN library such as Tensorflow or Keras.

.. code-block:: python 

	from keras.models import Sequential 
	from keras.layers import Dense, Activation

	# Import and define Keras network.
	params = [] 
	xshift = []
	xscale = []
	yshift = []
	yscale = []
	net = Sequential()
	with open("netstate.dat", "r") as f: 
		# Topology. 
		layers = int(f.readline())
		arch = [int(x) for x in f.readline().split()]
		
		# Scaling and shifting. 
		xscale = [float(x) for x in f.readline().split()]
		xshift = [float(x) for x in f.readline().split()]
		yscale = [float(x) for x in f.readline().split()]
		yshift = [float(x) for x in f.readline().split()]
		
		# Weights and biases.    
		for i in range(1, layers):
			b = []
			for j in range(arch[i]):
				b.append(float(f.readline()))
			b = np.array(b) 
			
			w = []
			for j in range(arch[i]*arch[i-1]):
				w.append(float(f.readline()))
			w = np.array(w).reshape(arch[i-1], arch[i])
			
			params.append(w)
			params.append(b)
			
			if i==1:
				net.add(Dense(arch[i], activation="tanh", input_dim=arch[i-1]))
			elif i==layers-1:
				net.add(Dense(arch[i], activation="linear"))
			else:
				net.add(Dense(arch[i], activation="tanh"))

	net.set_weights(params)

The network can then be evaluated on a high resolution grid and plotted. 

.. code-block:: python 

	# Define new high-resolution grid. 
	x = np.linspace(-np.pi, np.pi, 500, endpoint=True)
	y = np.linspace(-np.pi, np.pi, 500, endpoint=True)
	xg, yg = np.meshgrid(x, y)

	# Scale data. 
	xs = np.vstack((xg.flatten(), yg.flatten())).T
	xs = (xs - xshift)*xscale

	# Evaluate network. Unscale data.
	ys = net.predict(xs)
	ys = ys/yscale + yshift
	zg = -ys.reshape(500, 500)

	# Plot data.
	plt.figure(figsize=(12,10))
	zg = zg - np.max(zg)
	plt.contour(xg, yg, zg, linewidths=0.5, colors="k")
	plt.contourf(xg, yg, zg)
	cb = plt.colorbar()
	cb.set_label("G (kJ/mol)")
	plt.xlabel("$\phi$")
	plt.ylabel("$\psi$")

These examples and more are also found in the `online repository <https://github.com/hsidky/ann_sampling>`__.

Developer
^^^^^^^^^

Hythem Sidky.

.. warning:: 
	
	Please make sure to cite the paper [1]_ if you use this method!

References
^^^^^^^^^^

.. [1]  `H. Sidky and J. K. Whitmer`, J. Chem. Phys. **148**, 104111 (2018).
