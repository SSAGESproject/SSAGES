# This file is used to take a template input json file for forward flux
# and create a input file for ssages that uses multiple drivers.
import json 
import numpy as np
import copy

class LessPrecise(float):
    def __repr__(self):
        return str(self)

Numdrivers = 8
# Open template and load in the json data. 
root = {} 
with open('Template_Input.json') as f:
	root = json.load(f)

# Create vector of interfaces
centers = np.linspace(1.1, -1.1, num=100)

for i,center in enumerate(centers):
	centers[i] = round(center,4)

for center in centers:

	# Record all interfaces for forward flux
	key = "center"
	c = dict()
	c.setdefault(key,[])
	center = LessPrecise(round(center, 4))
	c["center"].append(center)

	root['method']['centers'].append(c)

# Add on the requested number of objects -1 because we are appending
for i in range(0,Numdrivers):
	root['driver'].append(copy.deepcopy(root['driver'][0]))

for i in range(0,Numdrivers):

	# Change the log file name so each driver uses a different log file
	root['driver'][i]['logfile'] = "none"
	
# Because appending remove original copy
root['driver'] = root['driver'][0:-1]

# Convert python dictionary into JSON file
with open('FF.json', 'w') as f:
		json.dump(root, f, indent=4, separators=(',', ': '))
