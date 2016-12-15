#! /usr/bin/env python

# This file is used to take a template input json file for forward flux
# and create a input file for ssages that uses multiple drivers.

import json 
import numpy as np
import copy

class LessPrecise(float):
    def __repr__(self):
        return str(self)

# User must set these variables        
nDrivers = 1
interfaces = np.array([-1.0, -0.95, -0.8, 0, 1])
trials = np.empty(5, dtype=int)
trials.fill(50)
# Use interfaces = np.linspace(<firstInterface>, <lastInterface>, num=<nInterfaces>) if you have many equally-spaced interfaces

# Open template and load in the json data. 
root = {} 
with open('Template_Input.json') as f:
	root = json.load(f)


#print interfaces
root['driver'][0]['method']['interfaces'] = interfaces.tolist()
root['driver'][0]['method']['trials'] = trials.tolist()


# Add on the requested number of objects -1 because we are appending
for i in range(0, nDrivers):
	root['driver'].append(copy.deepcopy(root['driver'][0]))

#for i in range(0,nDrivers):

	# Change the input filename so each driver uses a different input file file
#	root['driver'][i]['inputfile'] = "none"
	
# Because appending remove original copy
root['driver'] = root['driver'][0:-1]

# Convert python dictionary into JSON file
with open('Input-1proc.json', 'w') as f:
		json.dump(root, f, indent=4, separators=(',', ': '))
