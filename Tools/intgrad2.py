import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl

def intgrad2(fx,fy,nx,ny,dx,dy,intconst,per1,per2):
	
	rhs = np.ravel((fx,fy))
	
	Af=np.zeros((4*nx*ny,3))
	
	n=0
	#Equations in x
	for i in range(0,ny):
		#Leading edge
		Af[2*nx*i][0] = 2*nx*i/2
		if(per2):
			Af[2*nx*i][1] = nx*i+(nx-1)
		else:
			Af[2*nx*i][1] = nx*i
		Af[2*nx*i][2] = -0.5/dx
	
		Af[2*nx*i+1][0] = 2*nx*i/2
		Af[2*nx*i+1][1] = nx*i+1
		Af[2*nx*i+1][2] = 0.5/dx

		#Loop over inner space
		for j in range(1,nx-1):
			Af[2*nx*i+2*j][0] = int((2*nx*i+2*j)/2)
			Af[2*nx*i+2*j][1] = nx*i+j
			Af[2*nx*i+2*j][2] = -1/dx
	
			Af[2*nx*i+2*j+1][0] = int((2*nx*i+2*j)/2)
			Af[2*nx*i+2*j+1][1] = nx*i+j+1
			Af[2*nx*i+2*j+1][2] = 1/dx

		#Trailing edge
		Af[2*nx*(i+1)-2][0] = int((2*nx*(i+1)-2)/2)
		Af[2*nx*(i+1)-2][1] = nx*i+(nx-2)
		Af[2*nx*(i+1)-2][2] = -0.5/dx
	
		Af[2*nx*(i+1)-1][0] = int((2*nx*(i+1)-2)/2)
		if(per2):
			Af[2*nx*(i+1)-1][1] = nx*i
		else:
			Af[2*nx*(i+1)-1][1] = nx*i+(nx-1)
		Af[2*nx*(i+1)-1][2] = 0.5/dx
	
	
	n=2*nx*ny
	
	#Equations in y
	#Leading edge
	for j in range(0,nx):

		Af[2*j+n][0] = 2*j/2 + n/2
		
		if(per1):
			Af[2*j+n][1] = (ny-1)*nx+j
		else:
			Af[2*j+n][1] = j
		Af[2*j+n][2] = -0.5/dy
	
		Af[2*j+n+1][0] = 2*j/2 + n/2
		Af[2*j+n+1][1] = j+nx
		Af[2*j+n+1][2] = 0.5/dy
	
	#Loop over inner space
	for i in range(1,ny-1):
		for j in range(0,nx):
			
			Af[2*nx*i+2*j+n][0] = int((2*nx*i+2*j+n)/2)
			Af[2*nx*i+2*j+n][1] = j+(i)*nx
			Af[2*nx*i+2*j+n][2] = -1/dy
	
			Af[2*nx*i+2*j+n+1][0] = int((2*nx*i+2*j+n)/2)
			Af[2*nx*i+2*j+n+1][1] = j+(i+1)*nx
			Af[2*nx*i+2*j+n+1][2] = 1/dy
			a=2*nx*i+2*j+n+1
	n=n+2*(ny-1)*nx
	
	#Trailing edge
	for j in range(0,nx):
		Af[2*j+n][0] = int((2*j+n)/2)
		Af[2*j+n][1] = (ny-2)*nx+j
		Af[2*j+n][2] = -0.5/dy
	
		Af[2*j+n+1][0] = int((2*j+n)/2)
		if(per1):
			Af[2*j+n+1][1] = j
		else:
			Af[2*j+n+1][1] = (ny-1)*nx+j
		Af[2*j+n+1][2] = 0.5/dy


	#Boundary conditions
	Af[0][2]=1
	Af[1][:]=0
	rhs[0] = intconst

	#Solve
	A=sps.csc_matrix((Af[:,2],(Af[:,0],Af[:,1])),shape=(2*ny*nx,ny*nx))
	fhat=spsl.lsmr(A,rhs)
	fhat=fhat[0]
	
	return fhat
