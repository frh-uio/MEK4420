from numpy import *
from matplotlib.pyplot import *
import warnings
warnings.filterwarnings('ignore')

N_list=[100, 200, 400, 1000] # Number of segments

def square(a, N):
	"""
	Function for calculating the potential and added 
	mass of coefficients for a square.
	"""
	print('For a square with side length of 2a = %.2f' % (2*a))
	
	m11_ = []
	m22_ = []
	m66_ = []

	for N in N_list:
		print 'With %d segments we have:' % N
		# Matrix for storing values achieved from the lhs of the equation
		A = zeros((N,N))  

		# Array(vector) on the rhs os the equation 
		B11 = zeros(N)                                          
		B22 = zeros(N)											
		B66 = zeros(N)

		x = zeros(N+1)
		y = zeros(N+1)

		S1 = -a # Negative half of a side
		S2 = a  # Positive half of the side

		N = N/4 * 4 	# Total number of segments for all the four sides										 
		N1 = N/4 
				# Number of segments for one side										 
		x = zeros(N+1)
		y = zeros(N+1)

		#making points in all four sides
		for i in range(N1+1):
			inc = a*(1-(cos(pi/N1*i))) # from a to 2a

			# Bottom side (Constant y=-a while x goes from -a to a)
			x[i] = S1 + inc
			y[i] = -a

			# Right side (Constant x=a while y goes from -a to a)
			x[i+N1] = a
			y[i+N1] = S1 + inc

			# Top side (Constant y=a while x goes from a to -a)
			x[i+2*N1] = -(S1 + inc)
			y[i+2*N1] = a

			# Left side (Constant x=-a while y goes from a to -a)
			x[i+3*N1] = -a
			y[i+3*N1] = -(S1 + inc)

		xbar = (x[1:] + x[:-1])/2.0
		ybar = (y[1:] + y[:-1])/2.0

		# Length of each segment dS = sqrt(d(x0,y0)^2 - d(x1,y1)^2)
		# (x,y) position minus the next (x,y) position
		dS = linalg.norm(array([x[1:],y[1:]])-array([x[:-1],y[:-1]]), axis=0) 

		# Normal vector components of the segments
		n1 = -(y[1:] - y[:-1])/dS 	#-dy/dS					       	 		 
		n2 = (x[1:] - x[:-1])/dS 	#dx/dS							               
		n6 = (xbar*n2 - ybar*n1) 

		for i in range(N):
			# array transpose to get [x,y] position of a point on the circumference
			# Distance from midpoint in segment i to the starting/ending point of the
			# next/current segment
			r1 = linalg.norm(array([x[:-1],y[:-1]]).T - array([xbar[i], ybar[i]]), axis=1)

			r2 = linalg.norm(array([x[1:],y[1:]]).T - array([xbar[i], ybar[i]]), axis=1)

			theta = -arccos((dS**2 - r2**2 - r1**2)/(-2*r2*r1)) 
			theta[isnan(theta)] = 0	

			#Calculates the right-hand side of the integral (24)
			h11 = (log(r1)+log(r2))*0.5*dS             
			h22 = (log(r1)+log(r2))*0.5*dS		     
			h66 = (log(r1)+log(r2))*0.5*dS

			#Adds the angles to the matrix A	    
			A[i] = theta # N matrices that are NxN
			fill_diagonal(A,-pi) #replace diagonal entries with -pi

			#Adds rhs to the B-arrays					    
			B11[i] = sum(n1*h11) 										 
			B22[i] = sum(n2*h22) 								 
			B66[i] = sum(n6*h66) 

		# Calculates phi for the three directions
		# Solve the linear matrix equation A*phi=B
		phi11 = linalg.solve(A,B11)							 
		phi22 = linalg.solve(A,B22)								 
		phi66 = linalg.solve(A,B66)

		# Calculates the added mass coefficients
		m11 = sum(phi11*n1*dS)
		m11_.append(m11)

		m22 = sum(phi22*n2*dS)
		m22_.append(m22)

		m66 = sum(phi66*n6*dS)
		m66_.append(m66)	

		Exact_m11 = 4.754*a**2
		Error_m11 = Exact_m11-m11
		print('The error in the added mass coefficient m11 is %.5f' % Error_m11)

		Exact_m22 = 4.754*a**2
		Error_m22 = Exact_m22-m22
		print('The error in the added mass coefficient m22 is %.5f' % Error_m22)

		Exact_m66 = 0.725*a**4
		Error_m66 = Exact_m66-m66
		print('The error in the added mass coefficient m66 is %.5f' % Error_m66)
		print

	return m11_, Exact_m11, m22_, Exact_m22, m66_, Exact_m66



