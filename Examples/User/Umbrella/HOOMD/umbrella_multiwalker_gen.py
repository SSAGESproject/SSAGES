"""
Author: Bradley Dice (bdice@umich.edu)
Description: The SSAGES Umbrella Sampling example for HOOMD-blue runs a
             simulation of a butane molecule. This script, modified from the
             LAMMPS example, generates input files for SSAGES and wham.
"""

# This file is used to take a template input json file for umbrella sampling
# and create a input file for ssages that uses multiple drivers
# each with different umbrella centers.
import json
import numpy as np

# Open template and load in the json data.
with open('umbrella_input.json') as f:
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

# generate input file
for i, center in enumerate(centers):

    # Define center for each walker
    root["methods"][0]["centers"].append([round(center, 6)])

    root["methods"][0]["output_file"].append('{:03d}_umbrella.dat'.format(i))

# generate wham file
with open('wham_metafile', 'w') as wham:
    wham.write("# This metafile is meant for use with wham: "
               "http://membrane.urmc.rochester.edu/content/wham\n")
    ksprings = root["methods"][0]["ksprings"][0]
    for i, center in enumerate(centers):
        wham.write('{:03d}_umbrella.dat\t{}\t{}\n'.format(
            i, round(center, 6), ksprings))


# Convert Python dictionary into JSON file
with open('multiwalker_umbrella_input.json', 'w') as f:
    json.dump(root, f, indent=4, separators=(',', ': '))
