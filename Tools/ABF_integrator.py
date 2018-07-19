import sys
import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl
import scipy.interpolate as interp
import os
import re
import numpy.linalg as npl
import matplotlib.pyplot as plt
import argparse

from intgrad2 import intgrad2
from intgrad3 import intgrad3fw
from intgrad3 import intgrad3bw
from intgrad3 import intgrad3double

def main():

	parser = argparse.ArgumentParser(description='This script post-process the F_out file that SSAGES print out, and create a FES from the estimation of the thermodynamic force')
	parser.add_argument('-i','--input',default='F_out',help='name of the file containing the force')
	parser.add_argument('-o','--output',default='G',help='file containing the free energy surface')
	parser.add_argument('-p','--periodicity',nargs='*',default=[False,False,False],help='periodicity of the CVs (True is periodic, False is not)')
	parser.add_argument('-n','--npoints',nargs='*',type=int,default=[200,200,200],help='number of points for the interpolation')
	parser.add_argument('-s','--scale',type=float,default=1,help='scale for the interpolation')
	parser.add_argument('-c','--contours',type=int,default=10,help='number of contours on final plot')

	args=parser.parse_args()

	inputfile=args.input
	outputname=args.output
	periodic=args.periodicity
	interpolate=args.npoints
	scale=args.scale
	ncontours=args.contours

	if interpolate == [0]:
		interpolate = 0

	periodic = [False if x == 'False' else x for x in periodic]	

	if not (os.path.isfile(inputfile)):
		print('Input file is missing. Aborting! Please use -h to ask for help!')
		exit()

	print('Input file is :\n', inputfile, '\n and output file is : \n',outputname)
	print('CV1 periodicity is set to \n',periodic[0],' \n')
	if(len(periodic)>1):
		print('CV2 periodicity is set to (if it exists) \n',periodic[1],' \n')
	if(len(periodic)>2):
		print('CV3 periodicity is set to (if it exists) \n',periodic[2],' \n')
	print('Free energy will be interpolated with \n',interpolate,'\n points \n')
	print('Free energy will be plotted with \n',ncontours,'\n contours')
	f = open(outputname+'_integrated.txt','w')	

	vfield = []

	with open(inputfile,"r") as infile:
		for line in reversed(infile.readlines()):
			if "The columns" in line:
				headerinfo=re.findall(r'\d+',line)
				break
			elif not (line.strip()==""):
				vfield.append(list(map(float,line.split())))
	

	vfield = np.flipud(np.asarray(vfield))
	shapeinfo = vfield.shape
	headerfound = 0
	try:
		headerinfo = [int(i) for i in headerinfo]
		headerfound = 1
	except NameError:
		print("Header not found. Proceeding with analysis, but the file read in is likely not created by ABF directly. Proceed with caution.")
	
	if shapeinfo[1] == 2:
		print("Two columns detected in input, which translates to a 1-dimensional grid of "+str(shapeinfo[0])+ " points.")
		
		if interpolate == 0:
			A=np.zeros((vfield.shape[0],vfield.shape[0]))
			dx = vfield[2,0]-vfield[1,0]
			X = vfield[:,0]
			b = vfield[:,1]
			
		else:
			A=np.zeros((interpolate[0],interpolate[0]))
			dx = (vfield[-1,0]-vfield[0,0])/interpolate[0]
			X = np.linspace(vfield[0,0],vfield[-1,0],interpolate[0])
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

#		plt.show()
		plt.savefig(outputname+'_integrated.png')

	elif shapeinfo[1] == 4:	 
		if(len(periodic)<2):
			print('Please enter periodicity for each dimension.')
			exit()
		
		unique,counts=np.unique(vfield[:,0],return_counts=True)
		nx = len(unique)
		unique,counts=np.unique(vfield[:,1],return_counts=True)
		ny = len(unique)
	
		boundx = [vfield[0,0],vfield[-1,0]]
		boundy = [vfield[0,1],vfield[-1,1]]

		plot1 = plt.figure()
		plt.quiver(vfield[:,0],vfield[:,1],vfield[:,2],vfield[:,3],width = (0.0002*(vfield[nx*ny-1,0]-vfield[0,0])), headwidth=3)
		plt.axis([boundx[0],boundx[1],boundy[0],boundy[1]])
		plt.title('Gradient Field of the Free Energy')
		plot1.savefig(outputname+'_GradientField.png')

		if interpolate != 0:

			grid_x, grid_y = np.mgrid[boundx[0]:boundx[1]:interpolate[0]*1j,boundy[0]:boundy[1]:interpolate[1]*1j]
		
			dx = (boundx[1]-boundx[0])/interpolate[0]
			dy = (boundy[1]-boundy[0])/interpolate[1]
	
			nx = interpolate[0]
			ny = interpolate[1]

			fx=interp.griddata(vfield[:,0:2], vfield[:,2],(grid_x,grid_y), method='cubic')
			fy=interp.griddata(vfield[:,0:2], vfield[:,3],(grid_x,grid_y), method='cubic')

			fx=fx.T
			fy=fy.T
	
		else:			

			fy = vfield[:,3]
			fx = vfield[:,2]		

			dx = vfield[1,0] - vfield[0,0]
			dy = vfield[nx,1] - vfield[0,1]

			grid_x = vfield[:,0].reshape((nx,ny),order='F')
			grid_y = vfield[:,1].reshape((nx,ny),order='F')		
			

		fhat = intgrad2(fx,fy,nx,ny,dx,dy,1,periodic[0],periodic[1])
		fhat = fhat*scale

		zmin = min(fhat)
		zmax = max(fhat)
		fhat = fhat.reshape((nx,ny),order='F')
		
		X = grid_x
		Y = grid_y

		for k in range(0,ny):
			for n in range(0,nx): 
				f.write("{0:4f} {1:4f} {2:4f} \n".format(X[n,k],Y[n,k],-fhat[n,k]))


		plot3 = plt.figure()	
		plt.pcolormesh(X,Y,-fhat)
		plt.colorbar()
		plt.axis([min(vfield[:,0]),max(vfield[:,0]),min(vfield[:,1]),max(vfield[:,1])])
		plot3.savefig(outputname+'_EnergySurface.png')
	
		plot4=plt.figure()
		CS=plt.contour(X,Y,-fhat,antialiased=True,levels=np.linspace(np.amin(-fhat),np.amax(-fhat),30))

		plot4.savefig(outputname+'_contourmap.png')
		
		plot5 = plt.figure()
		black_contours=np.linspace(0,np.amax(-fhat-np.amin(-fhat)),ncontours+1)
		gray_contours=black_contours+np.amax(-fhat-np.amin(-fhat))/(2*(ncontours+1))
		plt.pcolormesh(X,Y,-fhat-np.amin(-fhat))
		plt.colorbar(label='kJ/mol')
		plt.axis([min(vfield[:,0]),max(vfield[:,0]),min(vfield[:,1]),max(vfield[:,1])])
		plt.contour(X,Y,(-fhat-np.amin(-fhat)),black_contours,colors='black',linewidths=0.75)
		plt.contour(X,Y,(-fhat-np.amin(-fhat)),gray_contours,colors='grey',linewidths=0.75)
		plt.tight_layout()
		plot5.savefig(outputname+'_merged.png')

	elif shapeinfo[1] == 6:
		if(len(periodic)<3):
			print('Please enter periodicity for each dimension.')
			exit()

		unique,counts=np.unique(vfield[:,0],return_counts=True)
		nx = len(unique)
		unique,counts=np.unique(vfield[:,1],return_counts=True)
		ny = len(unique)
		unique,counts=np.unique(vfield[:,2],return_counts=True)
		nz = len(unique)

		boundx = [vfield[0,0],vfield[-1,0]]
		boundy = [vfield[0,1],vfield[-1,1]]
		boundz = [vfield[0,2],vfield[-1,2]]
		
		

		if interpolate != 0:

			grid_x, grid_y, grid_z = np.mgrid[0:nx:interpolate[0]*1j,0:ny:interpolate[1]*1j,0:nz:interpolate[2]*1j]

			X, Y, Z = np.mgrid[boundx[0]:boundx[1]:interpolate[0]*1j,boundy[0]:boundy[1]:interpolate[1]*1j,boundz[0]:boundz[1]:interpolate[2]*1j]

			dx = X[1,0,0]-X[0,0,0]
			dy = Y[0,1,0]-Y[0,0,0]
			dz = Z[0,0,1]-Z[0,0,0]

			x = np.linspace(boundx[0],boundx[1],nx)
			y = np.linspace(boundy[0],boundy[1],ny)
			z = np.linspace(boundz[0],boundz[1],nz)

			datax = vfield[:,3].reshape((nx,ny,nz))
			datay = vfield[:,4].reshape((nx,ny,nz))
			dataz = vfield[:,5].reshape((nx,ny,nz))
	
			nx = interpolate[0]
			ny = interpolate[1]
			nz = interpolate[2]
		
			fx= interp2.map_coordinates(datax,np.vstack((np.ravel(grid_x),np.ravel(grid_y),np.ravel(grid_z))),order=3,mode='nearest')
			fy= interp2.map_coordinates(datay,np.vstack((np.ravel(grid_x),np.ravel(grid_y),np.ravel(grid_z))),order=3,mode='nearest')
			fz= interp2.map_coordinates(dataz,np.vstack((np.ravel(grid_x),np.ravel(grid_y),np.ravel(grid_z))),order=3,mode='nearest')

		else:			
		

			fx = vfield[:,3]
			fy = vfield[:,4]
			fz = vfield[:,5]		

			dx = vfield[1,0] - vfield[0,0]
			dy = vfield[nx,1] - vfield[0,1]
			dz = vfield[ny*nx,2] - vfield[0,2]

			grid_x = vfield[:,0].reshape((nx,ny,nz),order='F')
			grid_y = vfield[:,1].reshape((nx,ny,nz),order='F')
			grid_z = vfield[:,2].reshape((nx,ny,nz),order='F')
		
		fhat = intgrad3double(fx,fy,fz,nx,ny,nz,dx,dy,dz,1,periodic[0],periodic[1],periodic[2])
		fhat = fhat*scale

		zmin = min(fhat)
		zmax = max(fhat)
		fhat = fhat.reshape((nx,ny,nz),order='F')
		
		X = grid_x
		Y = grid_y
		Z = grid_z


		for l in range(0,nz):
			for k in range(0,ny):
				for j in range(0,nx): 
					f.write("{0:4f} {1:4f} {2:4f} {3:4f} \n".format(X[j,k,l],Y[j,k,l],Z[j,k,l],-fhat[j,k,l]))
	
	else:
		print("Input file does not contain the corrent number of dimensions. This code only supports 1D, 2D and 3D currently.")
		f.close()
		exit()
			
	f.close()

if __name__ == "__main__":
   main()
	 




