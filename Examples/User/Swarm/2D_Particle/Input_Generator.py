#! /usr/bin/env python

# Script for generating JSON input file for SSAGES Swarm method from a template
# input JSON file with multiple walkers (one per node on string)

import json
import numpy as np
import copy

# Open template and load in the json data.
root = {}
with open('Template_Input.json') as f:
	root = json.load(f)

#Number of processors/string nodes (make sure this matches everywhere)
num = 16

#Start and end location of CVs 1, 2, etc...
centers_1 = np.linspace(-0.98, 0.98, num)
centers_2 = np.linspace(-0.98, 0.98, num)

for i in range(num):
    x1 = round(centers_1[i], 3)
    x2 = round(centers_2[i], 3)
    root['methods'][0]['centers'].append([x1, x2])

# Convert python dictionary into JSON file
with open('Swarm.json', 'w') as f:
		json.dump(root, f, indent=4, separators=(',', ': '))
