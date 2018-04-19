#Circle

from Ellipse_function import *

colors = ['r', 'c', 'k', 'g', 'b']


phi11_, Exact_phi_, m11_, Exact_m11 = Ellipse_Circle(1, 1, N_list)

# Comparing analytical and numerical results for the potential
for i in range(len(N_list)):
	plot(phi11_[i], '-', color=colors[i], label='N=%d' % N_list[i])
	plot(Exact_phi_[i], '-', color=colors[i+1], label='Analytical solution')
	title('Analytical vs Numerical solution for the potential ' '$\phi_1$')
	xlabel('$Number$ $of$ $segments$ $N$', fontsize=18)
	ylabel('$\phi_1$', fontsize=24)
	legend(loc='upper right', numpoints = 1)
	savefig('Analytical_vs_numerical.png')
	show()


for i in range(len(N_list)):
	plot(1.0/N_list[i], m11_[i]/Exact_m11,'*', color=colors[i], markersize=16, linewidth=10, label='N=%d' %N_list[i])
	title('Computed added mass ' '$m_{11}$' ' divided by the analytical ' '$m_{11}$')
	hold(True)
	xlim((0, 0.011))
	ylabel('$m_{11}$' '/' ' Exact ' '$m_{11}$', fontsize=16)
	xlabel('1/N', fontsize=20)
	legend(loc='upper right', numpoints = 1)
show()


"""
Output:

For a circle with radius R0 = 1.00

With 100 segments we have:
The maximum error between the exact and the numerical potential is 0.01452
The maximum error between the exact and the numerical added mass coefficient m11 is 0.04460

With 200 segments we have:
The maximum error between the exact and the numerical potential is 0.00710
The maximum error between the exact and the numerical added mass coefficient m11 is 0.02204

With 400 segments we have:
The maximum error between the exact and the numerical potential is 0.00351
The maximum error between the exact and the numerical added mass coefficient m11 is 0.01095

With 1000 segments we have:
The maximum error between the exact and the numerical potential is 0.00139
The maximum error between the exact and the numerical added mass coefficient m11 is 0.00437
"""