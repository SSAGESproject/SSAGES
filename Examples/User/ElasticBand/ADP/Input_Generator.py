#! /usr/bin/env python

# Script for generating JSON input file for SSAGES nudged elastic band method from a template
# input JSON file with multiple walkers (one per node on string)

import json
import numpy as np
import copy

# Open template and load in the json data.
root = {}
with open('Template_Input.json') as f:
	root = json.load(f)

# number of nodes to use for nudged elastic band
num = 22

# Create vector of Torsional angle range. -pi to pi
centers_1 = np.linspace(-1.5, 0.25, num)
centers_2 = np.linspace(2.0, -2.0, num)

for i in range(num):
    x1 = round(centers_1[i], 3)
    x2 = round(centers_2[i], 3)
    root['methods'][0]['centers'].append([x1, x2])

# Convert python dictionary into JSON file
with open('ElasticBand.json', 'w') as f:
		json.dump(root, f, indent=4, separators=(',', ': '))
