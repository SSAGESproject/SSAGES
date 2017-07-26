#! /usr/bin/env python

# This file is used to take a template input json file for forward flux
# and create a input file for ssages that uses multiple drivers.

import json 
import numpy as np
from random import randint
import os

class LessPrecise(float):
    def __repr__(self):
        return str(self)

def gromacs_random_number(gromacs_filename, nWalkers):
    """Takes a gromacs input file, finds the seed number of langevin, replaces it with a random number, and
	generates new input files
	"""
    randon_numbers = set()
    while (len(randon_numbers) < nWalkers):
        randon_numbers.add(randint(0, 100000))
    randon_numbers = list(randon_numbers)

    in_gromacs_file = open(gromacs_filename, 'r')

    for i in range(0, nWalkers):
        out_gromacs_file = open(gromacs_filename + '-' + str(i) + ".mdp", 'w')
        for line in in_gromacs_file:
            if line.find('gen_seed') != -1:
                line = line.strip()
                if not line.startswith('#'):
                    columns = line.split()
                    if len(columns) >= 3:
                        line = line.replace(columns[2], str(randon_numbers[i]))
                        line = line + '\n'
            out_gromacs_file.write(line)

        in_gromacs_file.seek(0)
        out_gromacs_file.close()

    in_gromacs_file.close()


# User must set these variables
nWalkers = 2
input_filename = "nvt"
interfaces = np.array([-2.61, -2.44, -2.26, -2.09, -1.04])
trials = np.empty(5, dtype=int)
trials.fill(50)
# Use interfaces = np.linspace(<firstInterface>, <lastInterface>, num=<nInterfaces>) if you have many equally-spaced interfaces

#Generate the new gromacs input files
gromacs_random_number(input_filename, nWalkers)

# Open template and load in the json data. 
root = {}
with open('Template_Input.json') as f:
    root = json.load(f)


root["walkers"] = nWalkers
root['methods'][0]['interfaces'] = interfaces.tolist()
root['methods'][0]['trials'] = trials.tolist()

for i in range(0, nWalkers):
    #Generate the gromacs tpr files
    os.system("./gmx_mpi grompp -f nvt-" + str(i) + ".mdp -p topol.top -c adp.gro -o adp" + str(i) + '.tpr')

# Convert python dictionary into JSON file
with open('Input-2walkers.json', 'w') as f:
        json.dump(root, f, indent=4, separators=(',', ': '))
