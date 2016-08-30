import sys 
import numpy as np
import scipy as sp
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import os

import numpy.matlib as matlib
import numpy.linalg as npl

from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import matplotlib.pyplot as plt
np.set_printoptions(threshold=np.inf)


def intgrad2(fx,fy,nx,ny,dx,dy,intconst):
	
	rhs = np.ravel((fy,fx))
	
	Af=np.zeros((4*nx*ny,3))
	
	n=0
	#Equations in x
	for i in range(0,nx):
		#Leading edge
		Af[2*ny*i][0] = 2*ny*i/2
		Af[2*ny*i][1] = ny*i
		Af[2*ny*i][2] = -1/dx
	
		Af[2*ny*i+1][0] = 2*ny*i/2
		Af[2*ny*i+1][1] = ny*i+1
		Af[2*ny*i+1][2] = 1/dx

		#Loop over inner space
		for j in range(1,ny-1):
			Af[2*ny*i+2*j][0] = int((2*ny*i+2*j)/2)
			Af[2*ny*i+2*j][1] = ny*i+j-1
			Af[2*ny*i+2*j][2] = -0.5/dx
	
			Af[2*ny*i+2*j+1][0] = int((2*ny*i+2*j)/2)
			Af[2*ny*i+2*j+1][1] = ny*i+j+1
			Af[2*ny*i+2*j+1][2] = 0.5/dx

		#Trailing edge
		Af[2*ny*(i+1)-2][0] = int((2*ny*(i+1)-2)/2)
		Af[2*ny*(i+1)-2][1] = ny*i+(ny-2)
		Af[2*ny*(i+1)-2][2] = -1/dx
	
		Af[2*ny*(i+1)-1][0] = int((2*ny*(i+1)-2)/2)
		Af[2*ny*(i+1)-1][1] = ny*i+(ny-1)
		Af[2*ny*(i+1)-1][2] = 1/dx
	
	
	n=2*nx*ny
	
	#Equations in y
	#Leading edge
	for j in range(0,ny):

		Af[2*j+n][0] = 2*j/2 + n/2
		Af[2*j+n][1] = j
		Af[2*j+n][2] = -1/dy
	
		Af[2*j+n+1][0] = 2*j/2 + n/2
		Af[2*j+n+1][1] = j+ny
		Af[2*j+n+1][2] = 1/dy
	
	#Loop over inner space
	for i in range(1,nx-1):
		for j in range(0,ny):
			
			Af[2*ny*i+2*j+n][0] = int((2*ny*i+2*j+n)/2)
			Af[2*ny*i+2*j+n][1] = j+(i-1)*ny
			Af[2*ny*i+2*j+n][2] = -0.5/dy
	
			Af[2*ny*i+2*j+n+1][0] = int((2*ny*i+2*j+n)/2)
			Af[2*ny*i+2*j+n+1][1] = j+(i+1)*ny
			Af[2*ny*i+2*j+n+1][2] = 0.5/dy
			a=2*ny*i+2*j+n+1
	n=n+2*(nx-1)*ny
	
	#Trailing edge
	for j in range(0,ny):
		Af[2*j+n][0] = int((2*j+n)/2)
		Af[2*j+n][1] = (nx-2)*ny+j
		Af[2*j+n][2] = -1/dy
	
		Af[2*j+n+1][0] = int((2*j+n)/2)
		Af[2*j+n+1][1] = (nx-1)*ny+j
		Af[2*j+n+1][2] = 1/dy


	#Boundary conditions
	Af[0][2]=1
	Af[1][:]=0
	rhs[0] = intconst

	#Solve
	A=sps.csc_matrix((Af[:,2],(Af[:,0],Af[:,1])),shape=(2*nx*ny,nx*ny))
	fhat=spsl.lsmr(A,rhs)
	fhat=fhat[0]
	
	return fhat



	    
	 

f = open('output.txt','w')
print("Please enter the filename in which the grid and vector component information is stored. This file should ONLY contain numbers in the following format:")
print("CV1 coord    CV2 coord  ...            d(A)/d(CV1)  d(A)/d(CV2)...")

filename = raw_input()

if os.path.isfile(filename):
	vfield = np.loadtxt(filename)
else:
	print("File not found.")
	exit()

shapeinfo = vfield.shape

if shapeinfo[1] == 2:
	print("Two columns detected in input, which translates to a 1-dimensional grid of ",shapeinfo[0]," points.")
	dx = vfield[2,0]-vfield[1,0]
	vfieldy = vfield[:,1]
	Asurf=[]
	for n in range (0,shapeinfo[0]):
		temp = vfieldy[0:n]
		temp2 = np.trapz(temp,dx=dx)
		Asurf.append(temp2)

	for k in range(0,shapeinfo[0]):
		f.write("{0:4f} {1:4f} \n".format(vfield[k,0],-Asurf[k]))
	
	plt.plot(vfield[:,0],-Asurf[:])
	plt.axis([min(vfield[:,0]),max(vfield[:,0]),min(vfield[:,1]),max(vfield[:,1])])
	plt.show()
	
	
elif shapeinfo[1] == 4:
	print("Four columns detected in input, which translates to a 2-dimensional grid of ",shapeinfo[0]," points.")
	print("Please enter CV1 number of bins:")
	nx = int(input())
	print("Please enter CV2 number of bins:")
	ny = int(input())
	if not shapeinfo[0] == (nx*ny):
		print("np.ndim of input file (",shapeinfo[0],") does not match given CV bin dimensions of ",nx," x ",ny," = ",nx*ny)
		exit()
	plot1 = plt.figure()
	plt.quiver(vfield[:,0],vfield[:,1],vfield[:,2],vfield[:,3],scale=2000)
	plt.axis([min(vfield[:,0]),max(vfield[:,0]),min(vfield[:,1]),max(vfield[:,1])])
	plt.title('Gradient Field of the Free Energy')

	fx = vfield[:,2]
	fy = vfield[:,3]

	dx = vfield[ny,0] - vfield[0,0]
	dy = vfield[1,1] - vfield[0,1]	

	fhat = intgrad2(fx,fy,nx,ny,dx,dy,1)
	plot2 = plt.figure()
	
	X = vfield[:,0].reshape((nx,ny))
	Y = vfield[:,1].reshape((nx,ny))
	zmin = min(fhat)
	zmax = max(fhat)
	fhat = fhat.reshape((nx,ny))

	for n in range(0,nx):
		for k in range(0,ny): 
			f.write("{0:4f} {1:4f} {2:4f} \n".format(X[n,k],Y[n,k],-fhat[n,k]))
	
	plt.pcolormesh(X,Y,-fhat)
	plt.colorbar()
	plt.axis([min(vfield[:,0]),max(vfield[:,0]),min(vfield[:,1]),max(vfield[:,1])])
	plt.show()

f.close()

else
	print("Input file does not contain the corrent number of dimensions. This code only support 1D and 2D currently.")
	f.close()
	exit()
