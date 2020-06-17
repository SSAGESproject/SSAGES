#!/usr/bin/python

#Execute this script in the folder with your node-00XX.log outputs to get a graph of initial and final images

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
import math
import scipy.interpolate
import csv
import sys, getopt
import sys
import getopt
import numpy as np
import scipy as sp
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import scipy.interpolate as interp
import os
import re
from collections import deque

import numpy.matlib as matlib
import numpy.linalg as npl

import matplotlib
matplotlib.use('Agg')
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.inf)


matplotlib.use('Agg')

def intgrad2(fx,fy,nx,ny,dx,dy,intconst):

	rhs = np.ravel((fy,fx))

	#interp(

	Af=np.zeros((4*nx*ny,3))

	n=0
	#Equations in x
	for i in range(0,nx):
		#Leading edge
		Af[2*ny*i][0] = 2*ny*i/2
		#Af[2*ny*i][1] = ny*i
		Af[2*ny*i][1] = ny*i+(ny-1)
		Af[2*ny*i][2] = -0.5/dx

		Af[2*ny*i+1][0] = 2*ny*i/2
		Af[2*ny*i+1][1] = ny*i+1
		Af[2*ny*i+1][2] = 0.5/dx

		#Loop over inner space
		for j in range(1,ny-1):
			Af[2*ny*i+2*j][0] = int((2*ny*i+2*j)/2)
			Af[2*ny*i+2*j][1] = ny*i+j
			Af[2*ny*i+2*j][2] = -1/dx

			Af[2*ny*i+2*j+1][0] = int((2*ny*i+2*j)/2)
			Af[2*ny*i+2*j+1][1] = ny*i+j+1
			Af[2*ny*i+2*j+1][2] = 1/dx

		#Trailing edge
		Af[2*ny*(i+1)-2][0] = int((2*ny*(i+1)-2)/2)
		Af[2*ny*(i+1)-2][1] = ny*i+(ny-2)
		Af[2*ny*(i+1)-2][2] = -0.5/dx

		Af[2*ny*(i+1)-1][0] = int((2*ny*(i+1)-2)/2)
		#Af[2*ny*(i+1)-1][1] = ny*i+(ny-1)
		Af[2*ny*(i+1)-1][1] = ny*i
		Af[2*ny*(i+1)-1][2] = 0.5/dx


	n=2*nx*ny

	#Equations in y
	#Leading edge
	for j in range(0,ny):

		Af[2*j+n][0] = 2*j/2 + n/2
		#Af[2*j+n][1] = j
		Af[2*j+n][1] = (nx-1)*ny+j
		Af[2*j+n][2] = -0.5/dy

		Af[2*j+n+1][0] = 2*j/2 + n/2
		Af[2*j+n+1][1] = j+ny
		Af[2*j+n+1][2] = 0.5/dy

	#Loop over inner space
	for i in range(1,nx-1):
		for j in range(0,ny):

			Af[2*ny*i+2*j+n][0] = int((2*ny*i+2*j+n)/2)
			Af[2*ny*i+2*j+n][1] = j+(i)*ny
			Af[2*ny*i+2*j+n][2] = -1/dy

			Af[2*ny*i+2*j+n+1][0] = int((2*ny*i+2*j+n)/2)
			Af[2*ny*i+2*j+n+1][1] = j+(i+1)*ny
			Af[2*ny*i+2*j+n+1][2] = 1/dy
			a=2*ny*i+2*j+n+1
	n=n+2*(nx-1)*ny

	#Trailing edge
	for j in range(0,ny):
		Af[2*j+n][0] = int((2*j+n)/2)
		Af[2*j+n][1] = (nx-2)*ny+j
		Af[2*j+n][2] = -0.5/dy

		Af[2*j+n+1][0] = int((2*j+n)/2)
		#Af[2*j+n+1][1] = (nx-1)*ny+j
		Af[2*j+n+1][1] = j
		Af[2*j+n+1][2] = 0.5/dy


	#Boundary conditions
	Af[0][2]=1
	Af[1][:]=0
	rhs[0] = intconst

	#Solve
	A=sps.csc_matrix((Af[:,2],(Af[:,0],Af[:,1])),shape=(2*nx*ny,nx*ny))
	fhat=spsl.lsmr(A,rhs)
	fhat=fhat[0]

	return fhat


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
    fileHandle.close()

for line in end_lines:
    end_data.append(line.split(' '))

for i in range(0, nodes):
    if i <= 9:
        fileName = "node-000{0}.log".format(i)
    else:
        fileName = "node-00{0}.log".format(i)

    fileHandle = open(fileName)
    lineList = fileHandle.readlines()
    start_lines.append(lineList[0])
    fileHandle.close()

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

#TEST
#end_xdegree = sorted(end_xdegree, key=float)
#end_ydegree = sorted(end_ydegree, key=float, reverse=True)

plt.title("ADP Isomerization")
plt.xlabel('X')
plt.ylabel('Y')

plt.plot(end_xdegree,end_ydegree, 'ro-', label='Final Images')
plt.plot(start_xdegree, start_ydegree, 'bo-', label='Initial Images')

plt.axis([-3.1, 3.1, -3.1, 3.1])

leg = plt.legend()

#Add contour lines from BFS if command line argument specified
if(free_energy == 'bfs'):

    N = 10000 #Number of points from BFS
    x, y, z = np.genfromtxt(r'basis.csv', unpack=True, delimiter=',')

    xi = np.linspace(x.min(), x.max(), N)
    yi = np.linspace(y.min(), y.max(), N)
#print "Calculating zi..."
    zi = scipy.interpolate.griddata((x, y), z, (xi[None, :], yi[:, None]), method='cubic')
#print zi

    plt.contour(xi, yi, zi, 45)

if(free_energy == 'abf'):
	inputfile = 'F_out'
	outputname = 'G'
	f = open(outputname+'_integrated.txt','w')
	vfield = []
	with open(inputfile,"r") as infile:
		for line in reversed(infile.readlines()):
			if "The columns" in line:
				headerinfo=re.findall(r'\d+',line)
				break
			elif not (line.strip()==""):
				vfield.append(map(float,line.split()))


	headerinfo = [int(i) for i in headerinfo]
	vfield = np.flipud(np.asarray(vfield))
	shapeinfo = vfield.shape
	if (headerinfo[0] == shapeinfo[0]) and (headerinfo[3]*2 == shapeinfo[1]):
		print "Columns reported and header information internally consistent. Continuing with analysis."
	else:
		print "WARNING! Header information does not match column information."


	if shapeinfo[1] == 2:
		print "Two columns detected in input, which translates to a 1-dimensional grid of "+str(shapeinfo[0])+ "points."
		dx = vfield[2,0]-vfield[1,0]
		vfieldy = vfield[:,1]
		Asurf=[]
		for n in range (0,shapeinfo[0]):
			temp = vfieldy[0:n]
			temp2 = np.trapz(temp,dx=dx)
			Asurf.append(-temp2)

		for k in range(0,shapeinfo[0]):
			f.write("{0:4f} {1:4f} \n".format(vfield[k,0],Asurf[k]))

		plt.plot(vfield[:,0],Asurf[:])
		plt.axis([min(vfield[:,0]),max(vfield[:,0]),min(vfield[:,1]),max(vfield[:,1])])
		#plt.show()
		plot.savefig(outputname+'_integrated.png')


	elif shapeinfo[1] == 4:
		print "Four columns detected in input, which translates to a 2-dimensional grid of "+str(shapeinfo[0])+" points."
		print "CV1 number of bins is "+str(headerinfo[1])
		nx = headerinfo[1]
		print"CV2 number of bins is "+str(headerinfo[2])
		ny = headerinfo[2]
		if not shapeinfo[0] == (nx*ny):
			print("np.ndim of input file (",shapeinfo[0],") does not match given CV bin dimensions of ",nx," x ",ny," = ",nx*ny)
			exit()
		#plot1 = plt.figure()
		#plt.quiver(vfield[:,0],vfield[:,1],vfield[:,2],vfield[:,3],width = (0.0002*(vfield[nx*ny-1,0]-vfield[0,0])), headwidth=3, scale=0.4)
		#plt.axis([min(vfield[:,0]),max(vfield[:,0]),min(vfield[:,1]),max(vfield[:,1])])
		#plt.title('Gradient Field of the Free Energy')
		#plot1.savefig(outputname+'_GradientField.png')

		bound = -vfield[0,0];

		grid_x, grid_y = np.mgrid[-bound:bound:500j,-bound:bound:500j]
		dx = 2*3.14/500
		dy = dx

		nx = 500
		ny = 500

		#print(grid_x)
		#print(grid_y)

		fx=interp.griddata(vfield[:,0:2], vfield[:,2],(grid_x,grid_y), method='cubic')
		fy=interp.griddata(vfield[:,0:2], vfield[:,3],(grid_x,grid_y), method='cubic')

		#print(fx)

		#fx = vfield[:,2]
		#fy = vfield[:,3]

		#dx = vfield[ny,0] - vfield[0,0]
		#dy = vfield[1,1] - vfield[0,1]

		fhat = intgrad2(fx,fy,nx,ny,dx,dy,1)
		fhat = fhat


		#print fhat

		#X = vfield[:,0].reshape((nx,ny))
		#Y = vfield[:,1].reshape((nx,ny))
		zmin = min(fhat)
		zmax = max(fhat)
		fhat = fhat.reshape((nx,ny))

		X = grid_x
		Y = grid_y

		for n in range(0,nx):
			for k in range(0,ny):
				f.write("{0:4f} {1:4f} {2:4f} \n".format(X[n,k],Y[n,k],-fhat[n,k]))


		#plot = plt.figure()

		#plt.pcolormesh(X,Y,-fhat)
		plt.contour(X,Y,-fhat,antialiased=True,levels=np.linspace(min(-fhat[:,0]),max(-fhat[:,0]),30))
		#plt.colorbar()
		plt.axis([-3.14, 3.14, -3.14, 3.14])
		#plot3.savefig(outputname+'_EnergySurface.png')

		#plot4=plt.figure()
		#CS=plt.contour(X,Y,-fhat,antialiased=True,levels=np.linspace(min(-fhat[:,0]),max(-fhat[:,0]),30))

		#plt.show()
		#plot4.savefig(outputname+'_contourmap.png')

	else:
		print("Input file does not contain the corrent number of dimensions. This code only supports 1D and 2D currently.")
		f.close()
		exit()

	f.close()
plt.savefig("fts_plot.png")
plt.show()
