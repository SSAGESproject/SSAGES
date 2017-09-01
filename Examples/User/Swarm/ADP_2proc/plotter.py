#!/usr/bin/python

#Execute this script in the folder with your node-00XX.log outputs to get a graph of initial and final images

#Run this script as python plotter.py "number of nodes" "none/bfs"
#For example python plotter.py 22 bfs runs this script for 22 nodes and overlays a bfs.csv surface (must be generated from bfs)
#Alternative python plotter.py 22 none simply plots the nodes without a free energy surface

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import math
import scipy.interpolate
import csv
import sys, getopt

nodes = int(sys.argv[1])
free_energy = sys.argv[2]
end_lines = list()
end_data = list()
start_lines = list()
start_data = list()

for i in range(0, nodes):
    if i <= 9:
        fileName = "node-000{0}.log".format(i)
    else:
        fileName = "node-00{0}.log".format(i)

    fileHandle = open(fileName)
    lineList = fileHandle.readlines()
    end_lines.append(lineList[-1])
    start_lines.append(lineList[0])
    fileHandle.close()

for line in end_lines:
    end_data.append(line.split(' '))

for line in start_lines:
    start_data.append(line.split(' '))

end_x = list()
end_y = list()
start_x = list()
start_y = list()

for line in end_data:
    end_x.append(line[2])
    end_y.append(line[4])

for line in start_data:
    start_x.append(line[2])
    start_y.append(line[4])


#im = mpimg.imread('ADP_Background.jpg')
#implot = plt.imshow(im)

#fig = plt.figure()

#ax1 = fig.add_subplot(111)

#Convert to degrees
converter = 1.0
end_xdegree = converter*np.array(end_x, dtype=float)
end_ydegree = converter*np.array(end_y, dtype=float)
start_xdegree = converter*np.array(start_x, dtype=float)
start_ydegree = converter*np.array(start_y, dtype=float)

plt.title("ADP Isomerization")
plt.xlabel('X')
plt.ylabel('Y')

plt.plot(end_xdegree,end_ydegree, 'ro-', label='Final Images')
plt.plot(start_xdegree, start_ydegree, 'bo-', label='Initial Images')

plt.axis([-3.1, 3.1, -3.1, 3.1])

leg = plt.legend()

#Add contour lines from BFS if command line argument specified
if(free_energy == 'bfs'):

    nx = 100;
    ny = 100;
    x, y, z = np.genfromtxt(r'basis.csv', unpack=True, delimiter=',')

    xi = x.reshape((nx,ny))
    yi = y.reshape((nx,ny))
    zi = z.reshape((nx,ny))

    plt.contour(xi, yi, zi, 45)

plt.savefig("swarm_plot.png")
#plt.show()
