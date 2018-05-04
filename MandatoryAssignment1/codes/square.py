#Square

from square_function import *

colors = ['r', 'c', 'k', 'g', 'b']

m11_, Exact_m11, m22_, Exact_m22, m66_, Exact_m66 = square(1, N_list)

#checking for convergence
for i in range(len(N_list)):
	plot(1.0/N_list[i], m11_[i]/Exact_m11,'*', color=colors[i], markersize=14, label='$m_{11} = m_{22}$')
	hold(True)
	plot(1.0/N_list[i], m66_[i]/Exact_m66,'^', color=colors[i], markersize=14, label='$m_{66}$')

	if i==0:
		legend(loc='best', numpoints = 1)

title('Convergence for the computed added mass coefficients ' '$m_{11}$' ', $m_{22}$' ' and ' '$m_{66}$')
hold(True)
xlim((0, 0.011))
ylim((0.90 , 1.005))
ylabel('Numerical added mass coefficients/Exact values', fontsize=14)
xlabel('  1/N , (same colors show the same number of segments)', fontsize=14)
show()

"""
#output

With 100 segments we have:
The error in the added mass coefficient m11 is 0.11908
The error in the added mass coefficient m22 is 0.11908
The error in the added mass coefficient m66 is 0.06435

With 200 segments we have:
The error in the added mass coefficient m11 is 0.06120
The error in the added mass coefficient m22 is 0.06120
The error in the added mass coefficient m66 is 0.03207

With 400 segments we have:
The error in the added mass coefficient m11 is 0.03116
The error in the added mass coefficient m22 is 0.03116
The error in the added mass coefficient m66 is 0.01620

With 1000 segments we have:
The error in the added mass coefficient m11 is 0.01272
The error in the added mass coefficient m22 is 0.01272
The error in the added mass coefficient m66 is 0.00673
"""