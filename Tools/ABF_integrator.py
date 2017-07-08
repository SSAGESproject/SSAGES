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


def intgrad2(fx,fy,nx,ny,dx,dy,intconst,per1,per2):
	
	#fx, fy flipped here, be careful if editing.
	rhs = np.ravel((fy,fx))
	
	Af=np.zeros((4*nx*ny,3))
	
	n=0
	#Equations in x
	for i in range(0,nx):
		#Leading edge
		Af[2*ny*i][0] = 2*ny*i/2
		if(per2):
			Af[2*ny*i][1] = ny*i+(ny-1)
		else:
			Af[2*ny*i][1] = ny*i
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
		if(per2):
			Af[2*ny*(i+1)-1][1] = ny*i
		else:
			Af[2*ny*(i+1)-1][1] = ny*i+(ny-1)
		Af[2*ny*(i+1)-1][2] = 0.5/dx
	
	
	n=2*nx*ny
	
	#Equations in y
	#Leading edge
	for j in range(0,ny):

		Af[2*j+n][0] = 2*j/2 + n/2
		
		if(per1):
			Af[2*j+n][1] = (nx-1)*ny+j
		else:
			Af[2*j+n][1] = j
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
		if(per1):
			Af[2*j+n+1][1] = j
		else:
			Af[2*j+n+1][1] = (nx-1)*ny+j
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


def main(argv):
	inputfile = ''
	outputname = ''
	periodic = [False,False]
	interpolate = 500
	scale = 1
	try:
		opts, args = getopt.getopt(argv,"hi:o:p:l:n:s",["ifile=","ofile=","periodic1=","periodic2=","interpolate=","scale="])
	except getopt.GetoptError:
		print 'ABF_integrator.py -i <inputfile> -o <outputname> --periodic1 <True/False> --periodic2 <True/False> --interpolate <integer> --scale <float>'
		sys.exit(2)
	for opt, arg in opts:
		if opt == '-h':
			print 'ABF_integrator.py -i <inputfile> -o <outputname> --periodic1 <True/False> --periodic2 <True/False> --interpolate <integer> --scale <float>'
			sys.exit()
		elif opt in ("-i", "--ifile"):
			inputfile = arg
		elif opt in ("-o", "--ofile"):
			outputname = arg
		elif opt in ("-p", "--periodic1"):
			if arg == 'True':
				periodic[0] = True
			elif arg == 'False':
				periodic[0] = False
			elif(1):
				print "--periodic1 option not understood. Please provide \"True\" or \"False\""
				sys.exit()
		elif opt in ("-l", "--periodic2"):
			if arg == 'True':
				periodic[1] = True
			elif arg == 'False':
				periodic[1] = False
			elif(1):
				print "--periodic2 option not understood. Please provide \"True\" or \"False\""
				sys.exit()
		elif opt in ("-n", "--interpolate"):
			interpolate = int(arg)
		elif opt in ("-s", "--scale"):
			scale = float(arg)
	if inputfile == '':
		print 'Defaulting to input \"F_out\"'
		inputfile = 'F_out'
	elif(1):
		print "Input file is \""+inputfile+"\""
	if outputname == '':
		print 'Defaulting to output name \"G\"'
		outputname = 'G'
	elif(1):
		print 'Output file is \"'+outputname+"\""
	if periodic[0] == False:
		print 'Periodicity in CV1 is set to false. You may change this by passing \"--periodic1 True\"'
	elif(1):
		print 'Periodicity in CV1 is set to true.'
	if periodic[1] == False:
		print 'Periodicity in CV2 is set to false. You may change this by passing \"--periodic2 True\"'
	elif(1):
		print 'Periodicity in CV2 is set to true.'
	if interpolate == 500:
		print 'Interpolation is set to its default value of 500. You may change this by passing \"--interpolate <integer>\" or turn it off with \"--interpolate 0\"'
	elif interpolate == 0:
		print 'Interpolation is turned off'
	elif(1):
		print 'Interpolation grid is set to '+str(interpolate)

	f = open(outputname+'_integrated.txt','w')	

	if not (os.path.isfile(inputfile)):
		print'File not found.'
		print 'Hint: usage is \"ABF_integrator.py -i <inputfile> -o <outputname> --periodic1 <True/False> --periodic2 <True/False> --interpolate <integer> --scale <float>\"'
		exit()		

	vfield = []
	with open(inputfile,"r") as infile:
		for line in reversed(infile.readlines()):
			if "The columns" in line:
				headerinfo=re.findall(r'\d+',line)
				break
			elif not (line.strip()==""):
				vfield.append(map(float,line.split()))
	
	vfield = np.flipud(np.asarray(vfield))
	shapeinfo = vfield.shape
	headerfound = 0
	try:
		headerinfo = [int(i) for i in headerinfo]
		headerfound = 1
		if (headerinfo[0] == shapeinfo[0]) and (headerinfo[-1]*2 == shapeinfo[1]):
			print "Columns reported and header information internally consistent. Continuing with analysis."
		else:
			print "WARNING! Header information does not match column information. Continuing with analysis, but there may be errors."
	except NameError:
		print "Header not found. Proceeding with analysis, but the file read in is likely not created by ABF directly. Proceed with caution."


	if shapeinfo[1] == 2:
		print "Two columns detected in input, which translates to a 1-dimensional grid of "+str(shapeinfo[0])+ " points."

		if interpolate == 0:
			A=np.zeros((vfield.shape[0],vfield.shape[0]))
			dx = vfield[2,0]-vfield[1,0]
			X = vfield[:,0]
			b = vfield[:,1]

		else:
			A=np.zeros((interpolate,interpolate))
			dx = (vfield[-1,0]-vfield[0,0])/interpolate
			X = np.linspace(vfield[0,0],vfield[-1,0],interpolate)
			b = np.interp(X,vfield[:,0],vfield[:,1])

		for i in range (1,A.shape[0]-1):
			A[i][i-1]=-1
			A[i][i]=1
	
		if(periodic[0]):
			A[0][-1] = -2
			A[0][0] = 2
			A[-1][-1] = 1
			A[-1][-2] = -1
		else:
			A[0][0] = -1
			A[0][1] = 1
			A[-1][-1] = 1
			A[-1][-2] = -1

			A[1][1] = 1
			A[1][0] = 0
			b[1] = 0			

		A=A/dx			
		
		Asurf,c,d,e = npl.lstsq(A,b)
		Asurf = -Asurf*scale
		
		dx = vfield[2,0]-vfield[1,0]
		vfieldy = vfield[:,1]

		for k in range(0,A.shape[0]):
			f.write("{0:4f} {1:4f} \n".format(X[k],Asurf[k]))
		minimum = min(Asurf[:])
		plt.plot(X,Asurf-minimum)

		#plt.show()
		plt.savefig(outputname+'_integrated.png')
	
	
	elif shapeinfo[1] == 4:
		print "Four columns detected in input, which translates to a 2-dimensional grid of "+str(shapeinfo[0])+" points."
		if(headerfound):
			print "CV1 number of bins is "+str(headerinfo[1])
			nx = headerinfo[1]
			print"CV2 number of bins is "+str(headerinfo[2])
			ny = headerinfo[2]
		
		else:
			print "Since header is missing, cannot read in CV dimensions."
			nx = input("Please enter CV1 number of bins:")
			ny = input("Please enter CV2 number of bins:")	

		if not shapeinfo[0] == (nx*ny):
			print("np.ndim of input file (",shapeinfo[0],") does not match given CV bin dimensions of ",nx," x ",ny," = ",nx*ny)
			exit()
		

		boundx = [vfield[0,0],vfield[-1,0]]
		boundy = [vfield[0,1],vfield[-1,1]]

		plot1 = plt.figure()
		plt.quiver(vfield[:,0],vfield[:,1],vfield[:,2],vfield[:,3],width = (0.0002*(vfield[nx*ny-1,0]-vfield[0,0])), headwidth=3)
		plt.axis([boundx[0],boundx[1],boundy[0],boundy[1]])
		plt.title('Gradient Field of the Free Energy')
		plot1.savefig(outputname+'_GradientField.png')

		if interpolate != 0:

			grid_x, grid_y = np.mgrid[boundx[0]:boundx[1]:interpolate*1j,boundy[0]:boundy[1]:interpolate*1j]
		
			dx = boundx[1]-boundx[0]/interpolate
			dy = boundy[1]-boundy[0]/interpolate
	
			nx = interpolate
			ny = interpolate

			fx=interp.griddata(vfield[:,0:2], vfield[:,2],(grid_x,grid_y), method='cubic')
			fy=interp.griddata(vfield[:,0:2], vfield[:,3],(grid_x,grid_y), method='cubic')
		else:			
		
			fx = vfield[:,2]
			fy = vfield[:,3]		

			dx = vfield[ny,0] - vfield[0,0]
			dy = vfield[1,1] - vfield[0,1]

			grid_x = vfield[:,0].reshape((nx,ny))
			grid_y = vfield[:,1].reshape((nx,ny))
		
		fhat = intgrad2(fx,fy,nx,ny,dx,dy,1,periodic[0],periodic[1])
		fhat = fhat*scale

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
		CS=plt.contour(X,Y,-fhat,antialiased=True,levels=np.linspace(np.amin(-fhat),np.amax(-fhat),30))

		plot4.savefig(outputname+'_contourmap.png')

	else:
		print("Input file does not contain the corrent number of dimensions. This code only supports 1D and 2D currently.")
		f.close()
		exit()

	f.close()

if __name__ == "__main__":
   main(sys.argv[1:])
	 




