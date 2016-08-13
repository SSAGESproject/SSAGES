# This file is used to take a template input json file for umbrella sampling
# and create a input file for ssages that uses multiple drivers
# each with different umbrella centers.
import json 
import numpy as np
import copy

nwalker=6
min_temp=100
max_temp=500
# Create vector of Temperatures
temps = np.linspace(min_temp, max_temp, num=nwalker)
print temps

#generate input file
inputfile = "Butane_SSAGES.in"

for i in range(0,nwalker):
    fin = open(inputfile)
    outputfile = inputfile + str(i)
    print outputfile
    fout = open(outputfile, "wt")
    for line in fin:
        fout.write( line.replace('ptemp', str(temps[i])) )
    fout.close()
    fin.close()


# Open template and load in the json data. 
root = {} 
with open('Template_Input.json') as f:
	root = json.load(f)

# Add on the requested number of objects -1 because we are appending
for i in range(0,len(temps)-1):
	root['driver'].append(copy.deepcopy(root['driver'][0]))

for i,temp in enumerate(temps):
        log = "log" + str(i)
	# Change the log file name so each driver uses a different log file
	root['driver'][i]['logfile'] = log

	# Change the inputfile name for each driver 
        outputfile = inputfile + str(i)
        root['driver'][i]['inputfile'] = outputfile

root2 = root['driver'][0:(nwalker-1)]
# Convert python dictionary into JSON file
with open('ParallelTemp.json', 'w') as f:
    json.dump(root, f, indent=4, separators=(',', ': '))
