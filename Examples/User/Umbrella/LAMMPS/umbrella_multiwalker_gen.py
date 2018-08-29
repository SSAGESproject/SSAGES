# This file is used to take a template input json file for umbrella sampling
# and create a input file for ssages that uses multiple drivers
# each with different umbrella centers.
import json
import numpy as np
import copy

# Open template and load in the json data.
root = {}
with open("umbrella_input.json") as f:
	root = json.load(f)

# Define number of walkers.
nwalkers = 12

# Create vector of Torsional angle range. -pi to pi
# We remove last value because torsional CV is periodic.
centers = np.linspace(-3.14159, 3.14159, num=nwalkers+1)[:-1]

# Define centers, and filenames as arrays. We will use the same
# spring constant for all.
root["walkers"] = nwalkers
root["methods"][0]["centers"] = []
root["methods"][0]["output_file"] = []

for i,center in enumerate(centers):

	# Define center for each walker
	root["methods"][0]["centers"].append([round(center, 3)])

	# Define output file for each walker
	root["methods"][0]["output_file"].append("node-{}-data.dat".format(i))

# Convert python dictionary into JSON file
with open('multiwalker_umbrella.json', 'w') as f:
		json.dump(root, f, indent=4, separators=(',', ': '))
