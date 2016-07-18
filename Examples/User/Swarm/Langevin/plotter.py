#!/usr/bin/python

import numpy as np
import matplotlib.pyplot as plt
import math

with open("2dstring.txt") as f:
    data = f.read()

data = data.split('\n')

x = [row.split(' ')[0] for row in data[:-1]]
y = [row.split(' ')[1] for row in data[:-1]]

fig = plt.figure()

ax1 = fig.add_subplot(111)

ax1.set_title("ADP Isomerization")
ax1.set_xlabel('X')
ax1.set_ylabel('Y')

ax1.plot(x,y, 'ro-', label='Images')

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
