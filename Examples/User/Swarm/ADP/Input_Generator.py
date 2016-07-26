#! /usr/bin/env python

# Script for generating JSON input file for SSAGES FTS method from a template
# input JSON file with multiple walkers (one per node on string)

import json
import numpy as np
import copy

# Open template and load in the json data.
root = {}
with open('Template_Input.json') as f:
	root = json.load(f)

#Number of processors/string nodes (make sure this matches everywhere)
num = 22

#Start and end location of CVs 1, 2, etc...
centers_1 = np.linspace(-3.1, 3.1, num)
centers_2 = np.linspace(3.1, -3.1, num)

# Add on the requested number of objects -1 because we are appending
for i in range(0,num - 1):
	root['driver'].append(copy.deepcopy(root['driver'][0]))

for i in range(num):
	# Change the log file name so each driver uses a different log file
	root['driver'][i]['logfile'] = "log"

	# Change the node's location
        root['driver'][i]['method']['centers'][0] = round(centers_1[i], 3)
        root['driver'][i]['method']['centers'][1] = round(centers_2[i], 3)

# Convert python dictionary into JSON file
with open('Swarm.json', 'w') as f:
		json.dump(root, f, indent=4, separators=(',', ': '))
