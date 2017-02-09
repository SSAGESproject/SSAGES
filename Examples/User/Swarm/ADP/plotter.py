##special_bonds	lj/coul 0.0 0.0 0.5!/usr/bin/python

#Execute this script in the folder with your node-00XX.log outputs to get a graph of initial and final images

import numpy as np
import matplotlib.pyplot as plt
import math

nodes = 12
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

fig = plt.figure()

ax1 = fig.add_subplot(111)

#Convert to degrees
converter = 57.295779513
end_xdegree = converter*np.array(end_x, dtype=float)
end_ydegree = converter*np.array(end_y, dtype=float)
start_xdegree = converter*np.array(start_x, dtype=float)
start_ydegree = converter*np.array(start_y, dtype=float)

ax1.set_title("ADP Isomerization")
ax1.set_xlabel('X')
ax1.set_ylabel('Y')

ax1.plot(end_xdegree,end_ydegree, 'ro-', label='Final Images')
ax1.plot(start_xdegree, start_ydegree, 'bo-', label='Initial Images')


leg = ax1.legend()

plt.show()
