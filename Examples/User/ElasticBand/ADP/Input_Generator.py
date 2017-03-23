# This file is used to take a template input json file for umbrella sampling
# and create a input file for ssages that uses multiple drivers
# each with different umbrella centers.
import json
import numpy as np
import copy

# Open template and load in the json data.
root = {}
with open('Template_Input.json') as f:
	root = json.load(f)

# Create vector of Torsional angle range. -pi to pi
centers1 = np.linspace(-2.0944, 1.50, num=23)
centers2 = np.linspace(2.0944, -1.50, num =23)
# Periodic torsional CV so remove the last value
centers1 = centers1[0:-1]
centers2 = centers2[0:-1]

# Add on the requested number of objects -1 because we are appending
for i in range(0,len(centers1)-1):
	root['driver'].append(copy.deepcopy(root['driver'][0]))

for i,center in enumerate(centers1):
	# Change the log file name so each driver uses a different log file
	#root['driver'][i]['logfile'] = "log_"+str(i)

	# Change the umbrella's location
	root['driver'][i]['method']['centers'][0] = round(center,4)
	root['driver'][i]['method']['centers'][1] = round(centers2[i],4)

# Convert python dictionary into JSON file
with open('ElasticBand.json', 'w') as f:
		json.dump(root, f, indent=4, separators=(',', ': '))
