from numpy import *
from matplotlib.pyplot import *
import warnings
warnings.filterwarnings('ignore')
import exceptions	

N_list = [100, 200, 400] # Number of segments								 

#N_list = [100, 200, 300] # Number of segments	

def Ellipse_Circle(a,b,N):
	"""
	Function for calculating the potential and added mass coefficients 
	for an ellipse where a is the half major axis and b is the half 
	minor axis. When a=b we have a circle.
	"""
	print

	if a==b:
		print('For a circle with radius R0 = %.2f' % a)
		print

	elif a>b:
		print('For an ellipse with half major axis a0 = %.2f' '\n'
				'and half minor axis b0 = %.2f' % (a, b))
		print

	else:
		raise exceptions.AssertionError('The major axis a must be larger than the minor axis b')
		# As required in the text for the assignment

	# Values for the potential and the added mass coefficients. achieved from 
	# different choice of N that needs to be stored in lists for later use
	phi11_ = []
	Exact_phi_ = []
	m11_ = []
	m22_ = []
	m66_ = []

	for N in N_list:
		print 'With %d segments we have:' % N

		# Matrix for storing values achieved from the lhs of the equation
		A = zeros((N,N))  

		# Array(vector) for storing values for the rhs of the equation 
		B11 = zeros(N)                                          
		B22 = zeros(N)											
		B66 = zeros(N)	

		dtheta = linspace(0, 2*pi, N+1) # division into N segments
		# Evaluation points
		x = a*cos(dtheta) # the x position of the start/end point of a segment
		y = b*sin(dtheta) # the y position of the start/end point of a segment

		# Collocation points
		# Centroid of each segment S(i)
		xbar = (x[1:] + x[:-1])/2.0
		ybar = (y[1:] + y[:-1])/2.0

		# Length of each segment dS = sqrt(d(x0,y0)^2 - d(x1,y1)^2)
		# (x,y) position - next (x,y) position
		dS = linalg.norm(array([x[1:],y[1:]])-array([x[:-1],y[:-1]]), axis=0) 

		# Normal vector components of the segments
		n1 = -(y[1:] - y[:-1])/dS 	# -dy/dS				       	 		 
		
		n2 = (x[1:] - x[:-1])/dS	# dx/dS	

		n6 = (xbar*n2 - ybar*n1) 	# x*n2 - y*n1

		for i in range(N):
			# Array transpose to get [x,y] position of a point on the circumference
			r1 = linalg.norm(array([x[:-1],y[:-1]]).T - array([xbar[i], ybar[i]]), axis=1)

			r2 = linalg.norm(array([x[1:],y[1:]]).T - array([xbar[i], ybar[i]]), axis=1)

			# Opening angle of segment S(i)
			theta = -arccos((dS**2 - r2**2 - r1**2)/(-2*r2*r1)) 
			theta[isnan(theta)] = 0	

			#Calculates the right-hand side of the integral eq.(24)
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
		phi11_.append(phi11)	

		phi22 = linalg.solve(A,B22)								 
		phi66 = linalg.solve(A,B66)

		# Calculates the added mass coefficients
		m11 = sum(phi11*n1*dS)
		m11_.append(m11)
		
		m22 = sum(phi22*n2*dS)
		m22_.append(m22)

		m66 = sum(phi66*n6*dS)
		m66_.append(m66)	

		if a>b:
			Exact_m11 = pi*b**2
			Error_m11 = Exact_m11-m11
			plot(N, Error_m11,'r*', markersize=6, label='m_{11}')
			hold(True)
			print('The error in the added mass coefficient m11 is %.5f' % Error_m11)

			Exact_m22 = pi*a**2
			Error_m22 = Exact_m22-m22
			plot(N, Error_m22, 'k^', markersize=6, label='m_{22}' )
			hold(True)
			print('The error in the added mass coefficient m22 is %.5f' % Error_m22)

			Exact_m66 = (1.0/8.0)*pi*((a**2) - (b**2))**2
			Error_m66 = Exact_m66-m66
			plot(N, Error_m66,'gs', markersize=6, label='m_{66}')
			hold(True)
			print('The error in the added mass coefficient m66 is %.5f' % Error_m66)
			print

		if a==b:
			Exact_phi = -(a**2*xbar)/(xbar**2 + ybar**2)
			Exact_phi_.append(Exact_phi)
			Error_phi = abs(Exact_phi-phi11).max()
			print('The maximum error between the exact and the numerical potential is %.5f' % Error_phi)
			
			Exact_m11 = pi*a**2
			Error_m11 = Exact_m11-m11
			print('The maximum error between the exact and the numerical added mass coefficient m11 is %.5f' % Error_m11)
			print
		if N == 20:
			legend(loc='best', numpoints = 1)
	title('Error as a function of the resolution')
	xlabel('Number of segments N')
	ylabel('Error in the computed values compared to the exact ones')
	xlim((0, 1600))
	ylim((-10, 600))
	show()

	return m11_, m22_, m66_, Exact_m11, Exact_m22, Exact_m66 	# in case of an ellipse
	#return phi11_, Exact_phi_, m11_, Exact_m11 		# in case of a circle



