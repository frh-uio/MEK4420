#Ellipse

from Ellipse_function import *


colors = ['r', 'g', 'k', 'c']

#m11_, m22_, m66_, Exact_m11, Exact_m22, Exact_m66 = Ellipse_Circle(2, 1, N_list)

m11_, m22_, m66_, Exact_m11, Exact_m22, Exact_m66 = Ellipse_Circle(10, 1, N_list)


#checking for convergence
for i in range(len(N_list)):
	plot(1.0/N_list[i], m11_[i]/Exact_m11,'*', color=colors[i], markersize=14, label='$m_{11}$')
	hold(True)
	plot(1.0/N_list[i], m22_[i]/Exact_m22,'s', color=colors[i], markersize=14, label='$m_{22}$')
	hold(True)
	plot(1.0/N_list[i], m66_[i]/Exact_m66,'^', color=colors[i], markersize=14, label='$m_{66}$')

	if i==0:
		legend(loc='best', numpoints = 1)

title('Convergence for the computed added mass coefficients ' '$m_{11}$' ', $m_{22}$' ' and ' '$m_{66}$')
hold(True)
xlim((0, 0.011))
ylim((0.91 , 1.005))
ylabel('Numerical added mass coefficients/Exact values', fontsize=14)
xlabel('  1/N , (same colors show the same number of segments)', fontsize=14)
show()


"""
For an ellipse with half major axis a0 = 2.00
and half minor axis b0 = 1.00

With 100 segments we have:
The error in the added mass coefficient m11 is 0.03396
The error in the added mass coefficient m22 is 0.26348
The error in the added mass coefficient m66 is 0.11808

With 200 segments we have:
The error in the added mass coefficient m11 is 0.01666
The error in the added mass coefficient m22 is 0.13118
The error in the added mass coefficient m66 is 0.05710

With 400 segments we have:
The error in the added mass coefficient m11 is 0.00825
The error in the added mass coefficient m22 is 0.06546
The error in the added mass coefficient m66 is 0.02806

With 1000 segments we have:
The error in the added mass coefficient m11 is 0.00328
The error in the added mass coefficient m22 is 0.02615
The error in the added mass coefficient m66 is 0.01110




For an ellipse with half major axis a0 = 10.00
and half minor axis b0 = 1.00

With 100 segments we have:
The error in the added mass coefficient m11 is 0.02544
The error in the added mass coefficient m22 is 23.62821
The error in the added mass coefficient m66 is 328.78271

With 200 segments we have:
The error in the added mass coefficient m11 is 0.01235
The error in the added mass coefficient m22 is 11.89113
The error in the added mass coefficient m66 is 162.94352

With 400 segments we have:
The error in the added mass coefficient m11 is 0.00608
The error in the added mass coefficient m22 is 5.96638
The error in the added mass coefficient m66 is 81.09127

With 1000 segments we have:
The error in the added mass coefficient m11 is 0.00241
The error in the added mass coefficient m22 is 2.39177
The error in the added mass coefficient m66 is 32.34333

"""