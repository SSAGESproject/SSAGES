#!/usr/bin/python

#Execute this script in the folder with your node-00XX.log outputs to get a graph of initial and final images

import numpy as np
import matplotlib.pyplot as plt
import math

nodes = 22
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

#print lines

#with open("2dstring.txt") as f:
#    data = f.read()

for line in end_lines:
    end_data.append(line.split(' '))

for line in start_lines:
    start_data.append(line.split(' '))

#print data
#print start_data

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

#print x
#print y

fig = plt.figure()

ax1 = fig.add_subplot(111)

#Convert to degrees
#converter = 57.295779513
converter = 1.0
end_xdegree = converter*np.array(end_x, dtype=float)
end_ydegree = converter*np.array(end_y, dtype=float)
start_xdegree = converter*np.array(start_x, dtype=float)
start_ydegree = converter*np.array(start_y, dtype=float)


ax1.set_title("Langevin Particle")
ax1.set_xlabel('X')
ax1.set_ylabel('Y')

ax1.plot(end_xdegree,end_ydegree, 'ro-', label='Final Images')
ax1.plot(start_xdegree, start_ydegree, 'ro-', label='Initial Images')


leg = ax1.legend()

#Add contour lines of Mueller potential
"""aa = [-1, -1, -6.5, 0.7]
bb = [0, 0, 11, 0.6]
cc = [-10, -10, -6.5, 0.7]
AA = [-200, -100, -170, 15]

XX = [1, 0, -0.5, -1]
YY = [0, 0.5, 1.5, 1]

xxx = np.linspace(-1.5, 1.2);
yyy = np.linspace(-0.2, 2);

xx, yy = np.meshgrid(xxx, yyy)

V1 = 0
for j in range(0,4):
    xx2 = np.power(xx-XX[j],2)
    xxyy = np.multiply(xx-XX[j], yy-YY[j])
    yy2 = np.power(yy-YY[j],2)
    V1 = V1 + AA[j]*np.exp(aa[j]*xx2+bb[j]*xxyy+cc[j]*yy2)

V1[V1 > 200] = 200;

CS = plt.contour(xx, yy, V1, 40)"""

plt.show()
