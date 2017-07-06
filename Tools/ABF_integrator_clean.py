import sys
import getopt
import numpy as np
import scipy as sp
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import scipy.interpolate as interp
import os
import re

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
		Af[2*ny*i][1] = ny*i
		#Af[2*ny*i][1] = ny*i+(ny-1)
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
		#Af[2*ny*(i+1)-1][1] = ny*i
		Af[2*ny*(i+1)-1][2] = 1/dx
	
	
	n=2*nx*ny
	
	#Equations in y
	#Leading edge
	for j in range(0,ny):

		Af[2*j+n][0] = 2*j/2 + n/2
		Af[2*j+n][1] = j
		#Af[2*j+n][1] = (nx-1)*ny+j
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
		#Af[2*j+n+1][1] = j
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


def main(argv):
	inputfile = ''
	outputname = ''
	try:
		opts, args = getopt.getopt(argv,"hi:o:",["ifile=","ofile="])
	except getopt.GetoptError:
		print 'ABF_integrator.py -i <inputfile> -o <outputname>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'ABF_integrator.py -i <inputfile> -o <outputname>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
			print 'Input file is '+inputfile
		elif opt in ("-o", "--ofile"):
			outputname = arg
			print 'Output name is '+outputname
	if inputfile == '':
		print 'Defaulting to input "F_out"'
		inputfile = 'F_out'
	elif():
		print "Input file is "+inputfile
	if outputname == '':
		print 'Defaulting to output name "G"'
		outputname = 'G'
	elif():
		print 'Output file is '+outputname

	f = open(outputname+'_integrated.txt','w')
	#print("Please enter the inputfile in which the grid and vector component information is stored. This file should ONLY contain numbers in the following format:")
	#print("CV1 coord    CV2 coord  ...            d(A)/d(CV1)  d(A)/d(CV2)...")

	#inputfile = raw_input()
	

	if not (os.path.isfile(inputfile)):
		print'File not found.'
		print 'Hint: usage is "ABF_integrator.py -i <inputfile> -o <outputname>"'
		exit()		

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
	if (headerinfo[0] == shapeinfo[0]) and (headerinfo[-1]*2 == shapeinfo[1]):
		print "Columns reported and header information internally consistent. Continuing with analysis."
	else:
		print "WARNING! Header information does not match column information."


	if shapeinfo[1] == 2:
		print "Two columns detected in input, which translates to a 1-dimensional grid of "+str(shapeinfo[0])+ "points."

		A=np.zeros((vfield.shape[0],vfield.shape[0]))
		dx = vfield[2,0]-vfield[1,0]
	
		A[0][0] = -1
		A[0][1] = 1
		A[vfield.shape[0]-1][vfield.shape[0]-1] = 1
		A[vfield.shape[0]-1][vfield.shape[0]-2] = -1
		for i in range (1,vfield.shape[0]-1):
			A[i][i-1]=-1
			A[i][i]=1
		A=A/(dx)
		A[1][0] = 0
		A[1][1] = 1
		A[1][2] = 0
		print A
		b = vfield[:,1]
		b[1]=0		
		
		Asurf,c,d,e = npl.lstsq(A,b)
		print Asurf
		Asurf = -Asurf
		
		dx = vfield[2,0]-vfield[1,0]
		vfieldy = vfield[:,1]

		for k in range(0,shapeinfo[0]):
			f.write("{0:4f} {1:4f} \n".format(vfield[k,0],Asurf[k]))
		minimum = min(Asurf[:])
		plt.plot(vfield[:,0],Asurf-minimum)

		#plt.show()
		plt.savefig(outputname+'_integrated.png')
	
	
	elif shapeinfo[1] == 4:
		print "Four columns detected in input, which translates to a 2-dimensional grid of "+str(shapeinfo[0])+" points."
		print "CV1 number of bins is "+str(headerinfo[1])
		nx = headerinfo[1]
		print"CV2 number of bins is "+str(headerinfo[2])
		ny = headerinfo[2]
		if not shapeinfo[0] == (nx*ny):
			print("np.ndim of input file (",shapeinfo[0],") does not match given CV bin dimensions of ",nx," x ",ny," = ",nx*ny)
			exit()
		plot1 = plt.figure()
		plt.quiver(vfield[:,0],vfield[:,1],vfield[:,2],vfield[:,3],width = (0.0002*(vfield[nx*ny-1,0]-vfield[0,0])), headwidth=3, scale=0.4)
		plt.axis([min(vfield[:,0]),max(vfield[:,0]),min(vfield[:,1]),max(vfield[:,1])])
		plt.title('Gradient Field of the Free Energy')
		plot1.savefig(outputname+'_GradientField.png')

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


		plot3 = plt.figure()
	
		plt.pcolormesh(X,Y,-fhat)
		plt.colorbar()
		plt.axis([min(vfield[:,0]),max(vfield[:,0]),min(vfield[:,1]),max(vfield[:,1])])
		plot3.savefig(outputname+'_EnergySurface.png')
	
		plot4=plt.figure()
		CS=plt.contour(X,Y,-fhat,antialiased=True,levels=np.linspace(min(-fhat[:,0]),max(-fhat[:,0]),30))

		#plt.show()
		plot4.savefig(outputname+'_contourmap.png')

	else:
		print("Input file does not contain the corrent number of dimensions. This code only supports 1D and 2D currently.")
		f.close()
		exit()

	f.close()

if __name__ == "__main__":
   main(sys.argv[1:])
	 




