#! /usr/bin/env python

# This file is used to take a template input json file for forward flux
# and create a input file for ssages that uses multiple drivers.

import json 
import numpy as np
import copy
from random import randint
import os

class LessPrecise(float):
    def __repr__(self):
        return str(self)

def gromacs_random_number(gromacs_filename, nDrivers):
    """Takes a gromacs input file, finds the seed number of langevin, replaces it with a random number, and
	generates new input files
	"""
    randon_numbers = set()
    while (len(randon_numbers) < nDrivers):
        randon_numbers.add(randint(0, 100000))
    randon_numbers = list(randon_numbers)

    in_gromacs_file = open(gromacs_filename, 'r')

    for i in range(0, nDrivers):
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
nDrivers = 2
gromacs_random_number("nvt", nDrivers)
interfaces = np.array([-2.61, -2.44, -2.26, -2.09, -1.04])
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
    # Change the input filename so each driver uses a different input file file
    #root['driver'][i]['inputfile'] = root['driver'][i]['inputfile'] + str(i) + '.tpr'
    #Generate the gromacs tpr files
    os.system("./gmx_mpi grompp -f nvt-" + str(i) + ".mdp -p topol.top -c adp.gro -o adp" + str(i) + '.tpr')
	
# Because appending remove original copy
root['driver'] = root['driver'][1:]

# Convert python dictionary into JSON file
with open('Input-2proc.json', 'w') as f:
        json.dump(root, f, indent=4, separators=(',', ': '))
