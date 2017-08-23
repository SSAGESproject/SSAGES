#! /usr/bin/env python

# This file is used to take a template input json file for forward flux
# and create a input file for ssages that uses multiple drivers.

import json 
import numpy as np
import copy
from random import randint

class LessPrecise(float):
    def __repr__(self):
        return str(self)

def lammps_random_number(lammps_filename, nWalkers):
    """Takes a lammps input file, finds the seed number of langevin, replaces it with a random number, and
	generates new input files
	"""
    randon_numbers = set()
    while (len(randon_numbers) < nWalkers):
        randon_numbers.add(randint(0, 100000))
    randon_numbers = list(randon_numbers)

    in_lammps_file = open(lammps_filename, 'r')

    for i in range(0, nWalkers):
        out_lammps_file = open(lammps_filename + '-' + str(i), 'w')
        for line in in_lammps_file:
            if line.find('velocity') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 4:
                        line = line.replace(columns[4], str(randon_numbers[i]))
                        line = line + '\n'
            out_lammps_file.write(line)

        in_lammps_file.seek(0)
        out_lammps_file.close()

    in_lammps_file.close()


# User must set these variables
nWalkers = 2
input_filename = "in.ADP_Test"
interfaces = np.array([-2.61, -2.44, -2.26, -2.09, -1.04])
trials = np.empty(5, dtype=int)
trials.fill(50)
# Use interfaces = np.linspace(<firstInterface>, <lastInterface>, num=<nInterfaces>) if you have many equally-spaced interfaces

#Generate the new lammps input files
lammps_random_number(input_filename, nWalkers)

# Open template and load in the json data. 
root = {}
with open('Template_Input.json') as f:
    root = json.load(f)


root["walkers"] = nWalkers
root['methods'][0]['interfaces'] = interfaces.tolist()
root['methods'][0]['trials'] = trials.tolist()

input = []
for i in range(0, nWalkers):
    input.append(input_filename+ '-' + str(i))
root['input'] = input
	
# Convert python dictionary into JSON file
with open('Input-2walkers.json', 'w') as f:
        json.dump(root, f, indent=4, separators=(',', ': '))
