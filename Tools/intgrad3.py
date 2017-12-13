import numpy as np
import scipy.sparse as sps
import scipy.sparse.linalg as spsl

def intgrad3bw(fx,fy,fz,nx,ny,nz,dx,dy,dz,intconst,per1,per2,per3):		
	
	rhs = np.ravel((fx,fy,fz))
	A=np.zeros((6*nx*ny*nz,3))

	#A[][0] is the eqn. number
	#A[][1] is the point(s) from which the numerical derivative is calculated
	#A[][2] is -0.5 +0.5 or similar depending on eqn.

	#So, each A[][0] should have two corresponding entries in A[] (one for plus one for minus)

	#Equations in x
	n=0
	for i in range(0,nz):
		for j in range(0,ny):
			#Leading edge
			A[2*(nx*j+i*nx*ny)][0] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)][1] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)][2] = -1.0/dx
		
			A[2*(nx*j+i*nx*ny)+1][0] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)+1][1] = (nx*j+i*nx*ny)+1
			A[2*(nx*j+i*nx*ny)+1][2] = 1.0/dx		
			
			#Loop over inner space 
			for k in range(1,nx-1):
				A[2*(nx*j+i*nx*ny)+2*k][0] = (nx*j+i*nx*ny)+k
				A[2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny)+k-1
				#A[2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dx
				A[2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dx
				
			
				A[2*(nx*j+i*nx*ny)+2*k+1][0] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny)+k+1
				A[2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dx
				A[2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dx
			
			#Trailing edge
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][0] = (nx*j+i*nx*ny)+(nx-1) 
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][1] = (nx*j+i*nx*ny)+(nx-1)-1
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][2] = -1.0/dx
		
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][0] = (nx*j+i*nx*ny)+(nx-1) 
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][1] = (nx*j+i*nx*ny)+(nx-1)
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][2] = 1.0/dx
		
			
	n=2*nx*ny*nz
	#Equations in y
	for i in range(0,nz):
		#Leading edge (j=0)
		for k in range(0,nx):
			A[n+2*(i*nx*ny)+2*k][0] = n/2 + (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k][1] = (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k][2] = -1.0/dy
		
			A[n+2*(i*nx*ny)+2*k+1][0] = n/2 + (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k+1][1] = (i*nx*ny) + k + nx
			A[n+2*(i*nx*ny)+2*k+1][2] = 1.0/dy
	
		#Inner space
		for j in range(1,ny-1):
			for k in range(0,nx):
				A[n+2*(nx*j+i*nx*ny)+2*k][0] = n/2 + (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k - nx
				#A[n+2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dy				
				A[n+2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dy
		
				A[n+2*(nx*j+i*nx*ny)+2*k+1][0] = n/2 + (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k + nx
				A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dy
				A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dy
	
		#Trailing edge (j=ny-1)
		for k in range(0,nx):
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][0] = n/2 + (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][1] = (nx*(ny-1)+i*nx*ny)+k-nx
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][2] = -1.0/dy
		
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][0] = n/2 + (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][1] = (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][2] = 1.0/dy

	n=4*nx*ny*nz
	#Equations in z
	#Leading edge (i=0)
	for j in range(0,ny):
		for k in range(0,nx):
			A[n+2*(nx*j)+2*k][0] = n/2 + (nx*j) + k
			A[n+2*(nx*j)+2*k][1] = (nx*j) + k
			A[n+2*(nx*j)+2*k][2] = -1.0/dz
		
			A[n+2*(nx*j)+2*k+1][0] = n/2 + (nx*j) + k
			A[n+2*(nx*j)+2*k+1][1] = (nx*j) + k + nx*ny
			A[n+2*(nx*j)+2*k+1][2] = 1.0/dz

	#Inner space
	for i in range(1,nz-1):
		for j in range(0,ny):
			for k in range(0,nx):
				A[n+2*(nx*j+i*nx*ny)+2*k][0] = n/2 + (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k - nx*ny
				#A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dz				
				A[n+2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dz
			
				A[n+2*(nx*j+i*nx*ny)+2*k+1][0] = n/2 + (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k +nx*ny
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dz
				A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dz

	#Trailing edge (i=nz-1)
	for j in range(0,ny):
		for k in range(0,nx):			
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][0] = n/2 + (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][1] = (nx*j+(nz-1)*nx*ny) + k - nx*ny
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][2] = -1.0/dz
			
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][0] = n/2 + (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][1] = (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][2] = 1.0/dz
		


	#Solve
	A[0][:] = 0
	A[1][:] = 0
	A[0][0]=0
	A[0][1]=0
	A[0][2]=1
	rhs[0] = intconst	

	AA=sps.csc_matrix((A[:,2],(A[:,0],A[:,1])),shape=(3*nx*ny*nz,nx*ny*nz))
	fhat=spsl.lsmr(AA,rhs)
	fhat=fhat[0]

	return fhat

def intgrad3fw(fx,fy,fz,nx,ny,nz,dx,dy,dz,intconst,per1,per2,per3):		
	
	rhs = np.ravel((fx,fy,fz))
	A=np.zeros((6*nx*ny*nz,3))

	#A[][0] is the eqn. number
	#A[][1] is the point(s) from which the numerical derivative is calculated
	#A[][2] is -0.5 +0.5 or similar depending on eqn.

	#So, each A[][0] should have two corresponding entries in A[] (one for plus one for minus)

	#Equations in x
	n=0
	for i in range(0,nz):
		for j in range(0,ny):
			#Leading edge
			A[2*(nx*j+i*nx*ny)][0] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)][1] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)][2] = -1.0/dx
		
			A[2*(nx*j+i*nx*ny)+1][0] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)+1][1] = (nx*j+i*nx*ny)+1
			A[2*(nx*j+i*nx*ny)+1][2] = 1.0/dx		
			
			#Loop over inner space 
			for k in range(1,nx-1):
				A[2*(nx*j+i*nx*ny)+2*k][0] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny)+k-1
				A[2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dx
				A[2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dx
				
			
				A[2*(nx*j+i*nx*ny)+2*k+1][0] = (nx*j+i*nx*ny)+k
				A[2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny)+k+1
				#A[2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dx
				A[2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dx
			
			#Trailing edge
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][0] = (nx*j+i*nx*ny)+(nx-1) 
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][1] = (nx*j+i*nx*ny)+(nx-1)-1
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][2] = -1.0/dx
		
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][0] = (nx*j+i*nx*ny)+(nx-1) 
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][1] = (nx*j+i*nx*ny)+(nx-1)
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][2] = 1.0/dx
		
			
	n=2*nx*ny*nz
	#Equations in y
	for i in range(0,nz):
		#Leading edge (j=0)
		for k in range(0,nx):
			A[n+2*(i*nx*ny)+2*k][0] = n/2 + (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k][1] = (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k][2] = -1.0/dy
		
			A[n+2*(i*nx*ny)+2*k+1][0] = n/2 + (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k+1][1] = (i*nx*ny) + k + nx
			A[n+2*(i*nx*ny)+2*k+1][2] = 1.0/dy
	
		#Inner space
		for j in range(1,ny-1):
			for k in range(0,nx):
				A[n+2*(nx*j+i*nx*ny)+2*k][0] = n/2 + (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k - nx
				#A[n+2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dy				
				A[n+2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dy
		
				A[n+2*(nx*j+i*nx*ny)+2*k+1][0] = n/2 + (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k + nx
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dy
				A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dy
	
		#Trailing edge (j=ny-1)
		for k in range(0,nx):

			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][0] = n/2 + (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][1] = (nx*(ny-1)+i*nx*ny)+k-nx
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][2] = -1.0/dy
		
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][0] = n/2 + (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][1] = (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][2] = 1.0/dy

	n=4*nx*ny*nz
	#Equations in z
	#Leading edge (i=0)
	for j in range(0,ny):
		for k in range(0,nx):
			A[n+2*(nx*j)+2*k][0] = n/2 + (nx*j) + k
			A[n+2*(nx*j)+2*k][1] = (nx*j) + k
			A[n+2*(nx*j)+2*k][2] = -1.0/dz
		
			A[n+2*(nx*j)+2*k+1][0] = n/2 + (nx*j) + k
			A[n+2*(nx*j)+2*k+1][1] = (nx*j) + k + nx*ny
			A[n+2*(nx*j)+2*k+1][2] = 1.0/dz

	#Inner space
	for i in range(1,nz-1):
		for j in range(0,ny):
			for k in range(0,nx):
				A[n+2*(nx*j+i*nx*ny)+2*k][0] = n/2 + (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k - nx*ny
				A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dz				
				A[n+2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dz
			
				A[n+2*(nx*j+i*nx*ny)+2*k+1][0] = n/2 + (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k +nx*ny
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dz
				A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dz

	#Trailing edge (i=nz-1)
	for j in range(0,ny):
		for k in range(0,nx):			
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][0] = n/2 + (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][1] = (nx*j+(nz-1)*nx*ny) + k - nx*ny
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][2] = -1.0/dz
			
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][0] = n/2 + (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][1] = (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][2] = 1.0/dz
		


	#Solve
	A[0][:] = 0
	A[1][:] = 0
	A[0][0]=0
	A[0][1]=0
	A[0][2]=1
	rhs[0] = intconst	

	AA=sps.csc_matrix((A[:,2],(A[:,0],A[:,1])),shape=(3*nx*ny*nz,nx*ny*nz))
	fhat=spsl.lsmr(AA,rhs)
	fhat=fhat[0]

	return fhat

def intgrad3double(fx,fy,fz,nx,ny,nz,dx,dy,dz,intconst,per1,per2,per3):		
	
	rhs = np.ravel((fx,fy,fz))
	A=np.zeros((6*nx*ny*nz,3))

	#A[][0] is the eqn. number
	#A[][1] is the point(s) from which the numerical derivative is calculated
	#A[][2] is -0.5 +0.5 or similar depending on eqn.

	#So, each A[][0] should have two corresponding entries in A[] (one for plus one for minus)

	#Equations in x
	n=0
	for i in range(0,nz):
		for j in range(0,ny):
			#Leading edge
			A[2*(nx*j+i*nx*ny)][0] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)][1] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)][2] = -1.0/dx
		
			A[2*(nx*j+i*nx*ny)+1][0] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)+1][1] = (nx*j+i*nx*ny)+1
			A[2*(nx*j+i*nx*ny)+1][2] = 1.0/dx		
			
			#Loop over inner space 
			for k in range(1,nx-1):
				A[2*(nx*j+i*nx*ny)+2*k][0] = (nx*j+i*nx*ny)+k
				A[2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny)+k-1
				#A[2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dx
				A[2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dx
				
			
				A[2*(nx*j+i*nx*ny)+2*k+1][0] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny)+k+1
				A[2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dx
				A[2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dx
			
			#Trailing edge
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][0] = (nx*j+i*nx*ny)+(nx-1) 
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][1] = (nx*j+i*nx*ny)+(nx-1)-1
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][2] = -1.0/dx
		
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][0] = (nx*j+i*nx*ny)+(nx-1) 
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][1] = (nx*j+i*nx*ny)+(nx-1)
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][2] = 1.0/dx
		
			
	n=2*nx*ny*nz
	#Equations in y
	for i in range(0,nz):
		#Leading edge (j=0)
		for k in range(0,nx):
			A[n+2*(i*nx*ny)+2*k][0] = n/2 + (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k][1] = (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k][2] = -1.0/dy
		
			A[n+2*(i*nx*ny)+2*k+1][0] = n/2 + (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k+1][1] = (i*nx*ny) + k + nx
			A[n+2*(i*nx*ny)+2*k+1][2] = 1.0/dy
	
		#Inner space
		for j in range(1,ny-1):
			for k in range(0,nx):
				A[n+2*(nx*j+i*nx*ny)+2*k][0] = n/2 + (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k - nx
				#A[n+2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dy				

				A[n+2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dy
		
				A[n+2*(nx*j+i*nx*ny)+2*k+1][0] = n/2 + (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k + nx
				A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dy
				A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dy
	
		#Trailing edge (j=ny-1)
		for k in range(0,nx):
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][0] = n/2 + (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][1] = (nx*(ny-1)+i*nx*ny)+k-nx
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][2] = -1.0/dy
		
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][0] = n/2 + (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][1] = (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][2] = 1.0/dy

	n=4*nx*ny*nz
	#Equations in z
	#Leading edge (i=0)
	for j in range(0,ny):
		for k in range(0,nx):
			A[n+2*(nx*j)+2*k][0] = n/2 + (nx*j) + k
			A[n+2*(nx*j)+2*k][1] = (nx*j) + k
			A[n+2*(nx*j)+2*k][2] = -1.0/dz
		
			A[n+2*(nx*j)+2*k+1][0] = n/2 + (nx*j) + k
			A[n+2*(nx*j)+2*k+1][1] = (nx*j) + k + nx*ny
			A[n+2*(nx*j)+2*k+1][2] = 1.0/dz

	#Inner space
	for i in range(1,nz-1):
		for j in range(0,ny):
			for k in range(0,nx):
				A[n+2*(nx*j+i*nx*ny)+2*k][0] = n/2 + (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k - nx*ny
				#A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dz				
				A[n+2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dz
			
				A[n+2*(nx*j+i*nx*ny)+2*k+1][0] = n/2 + (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k +nx*ny
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dz
				A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dz

	#Trailing edge (i=nz-1)
	for j in range(0,ny):
		for k in range(0,nx):			
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][0] = n/2 + (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][1] = (nx*j+(nz-1)*nx*ny) + k - nx*ny
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][2] = -1.0/dz
			
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][0] = n/2 + (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][1] = (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][2] = 1.0/dz
		


	#Solve
	A[0][:] = 0
	A[1][:] = 0
	A[0][0]=0
	A[0][1]=0
	A[0][2]=1
	rhs[0] = intconst	

	AA=sps.csc_matrix((A[:,2],(A[:,0],A[:,1])),shape=(3*nx*ny*nz,nx*ny*nz))
	fhat=spsl.lsmr(AA,rhs)
	fhat=fhat[0]

	rhs = np.ravel((fx,fy,fz))
	A=np.zeros((6*nx*ny*nz,3))

	#A[][0] is the eqn. number
	#A[][1] is the point(s) from which the numerical derivative is calculated
	#A[][2] is -0.5 +0.5 or similar depending on eqn.

	#So, each A[][0] should have two corresponding entries in A[] (one for plus one for minus)

	#Equations in x
	n=0
	for i in range(0,nz):
		for j in range(0,ny):
			#Leading edge
			A[2*(nx*j+i*nx*ny)][0] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)][1] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)][2] = -1.0/dx
		
			A[2*(nx*j+i*nx*ny)+1][0] = (nx*j+i*nx*ny)
			A[2*(nx*j+i*nx*ny)+1][1] = (nx*j+i*nx*ny)+1
			A[2*(nx*j+i*nx*ny)+1][2] = 1.0/dx		
			
			#Loop over inner space 
			for k in range(1,nx-1):
				A[2*(nx*j+i*nx*ny)+2*k][0] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny)+k-1
				A[2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dx
				A[2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dx
				
			
				A[2*(nx*j+i*nx*ny)+2*k+1][0] = (nx*j+i*nx*ny)+k
				A[2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny)+k+1
				#A[2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny)+k
				#A[2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dx
				A[2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dx
			
			#Trailing edge
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][0] = (nx*j+i*nx*ny)+(nx-1) 
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][1] = (nx*j+i*nx*ny)+(nx-1)-1
			A[2*(nx*j+i*nx*ny)+2*(nx-1)][2] = -1.0/dx
		
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][0] = (nx*j+i*nx*ny)+(nx-1) 
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][1] = (nx*j+i*nx*ny)+(nx-1)
			A[2*(nx*j+i*nx*ny)+2*(nx-1)+1][2] = 1.0/dx
		
			
	n=2*nx*ny*nz
	#Equations in y
	for i in range(0,nz):
		#Leading edge (j=0)
		for k in range(0,nx):
			A[n+2*(i*nx*ny)+2*k][0] = n/2 + (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k][1] = (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k][2] = -1.0/dy
		
			A[n+2*(i*nx*ny)+2*k+1][0] = n/2 + (i*nx*ny) + k
			A[n+2*(i*nx*ny)+2*k+1][1] = (i*nx*ny) + k + nx
			A[n+2*(i*nx*ny)+2*k+1][2] = 1.0/dy
	
		#Inner space
		for j in range(1,ny-1):
			for k in range(0,nx):
				A[n+2*(nx*j+i*nx*ny)+2*k][0] = n/2 + (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k - nx
				#A[n+2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dy				
				A[n+2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dy
		
				A[n+2*(nx*j+i*nx*ny)+2*k+1][0] = n/2 + (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k + nx
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dy
				A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dy
	
		#Trailing edge (j=ny-1)
		for k in range(0,nx):
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][0] = n/2 + (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][1] = (nx*(ny-1)+i*nx*ny)+k-nx
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k][2] = -1.0/dy
		
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][0] = n/2 + (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][1] = (nx*(ny-1)+i*nx*ny)+k
			A[n+2*(nx*(ny-1)+i*nx*ny)+2*k+1][2] = 1.0/dy

	n=4*nx*ny*nz
	#Equations in z
	#Leading edge (i=0)
	for j in range(0,ny):
		for k in range(0,nx):
			A[n+2*(nx*j)+2*k][0] = n/2 + (nx*j) + k
			A[n+2*(nx*j)+2*k][1] = (nx*j) + k
			A[n+2*(nx*j)+2*k][2] = -1.0/dz
		
			A[n+2*(nx*j)+2*k+1][0] = n/2 + (nx*j) + k
			A[n+2*(nx*j)+2*k+1][1] = (nx*j) + k + nx*ny
			A[n+2*(nx*j)+2*k+1][2] = 1.0/dz

	#Inner space
	for i in range(1,nz-1):
		for j in range(0,ny):
			for k in range(0,nx):
				A[n+2*(nx*j+i*nx*ny)+2*k][0] = n/2 + (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k - nx*ny
				A[n+2*(nx*j+i*nx*ny)+2*k][1] = (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k][2] = -0.5/dz				
				A[n+2*(nx*j+i*nx*ny)+2*k][2] = -1.0/dz
			
				A[n+2*(nx*j+i*nx*ny)+2*k+1][0] = n/2 + (nx*j+i*nx*ny) + k
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k
				A[n+2*(nx*j+i*nx*ny)+2*k+1][1] = (nx*j+i*nx*ny) + k +nx*ny
				#A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 0.5/dz
				A[n+2*(nx*j+i*nx*ny)+2*k+1][2] = 1.0/dz

	#Trailing edge (i=nz-1)
	for j in range(0,ny):
		for k in range(0,nx):			
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][0] = n/2 + (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][1] = (nx*j+(nz-1)*nx*ny) + k - nx*ny
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k][2] = -1.0/dz
			
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][0] = n/2 + (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][1] = (nx*j+(nz-1)*nx*ny) + k
			A[n+2*(nx*j+(nz-1)*nx*ny)+2*k+1][2] = 1.0/dz
		


	#Solve
	A[0][:] = 0
	A[1][:] = 0
	A[0][0]=0
	A[0][1]=0
	A[0][2]=1
	rhs[0] = intconst	

	AA=sps.csc_matrix((A[:,2],(A[:,0],A[:,1])),shape=(3*nx*ny*nz,nx*ny*nz))
	fhat2=spsl.lsmr(AA,rhs)
	fhat2=fhat2[0]

	return (fhat+fhat2)/2.0 
	
