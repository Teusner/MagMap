import numpy as np
import vibes
from pyibex import *
from tubex_lib import *
import time
from tank import Tank


if __name__ == "__main__":
	
	# Time
	t0 = 0
	tf = 10
	h = 1/20
	T = np.arange(t0, tf, h)

	# Tubes
	phi = IntervalVector([Interval(-np.pi/2, np.pi/2), Interval(-1, 1)])
	n = 10
	A, B = IntervalVector(n), IntervalVector(n)

	A[0] = Interval(-1, 1)
	print(A)

	# Contractor
	ctc_f = CtcFunction(Function("phi", "dphi", "dphi+sin(phi)"))

	for t in T:
		lb, ub, r = phi[0].lb(), phi[0].ub(), phi[0].rad()
		for k in range(n-1):
			A[k] = Interval(lb + k * r / n, lb + (k+1) * r / n)
		for i, p in enumerate(A):
			cn = ContractorNetwork()
			dp = Interval(-1, 1)
			cn.add(ctc_f, [p, dp])
			cn.add(ctc.deriv, [p, dp])
			cn.contract()
			B[i] = p
		print(B)

	# Vehicle
	# vehicle = Tank(X=np.array([[0.], [0.], [0.], [np.pi/4]]), h=h)

	# Loop
	# for t in T:
	# 	U = np.array([[1], [1]])
	# 	vehicle.step(U)
	# 	vehicle.draw()
		
		# Cn adding data
		# cn.add_data(v, t, Interval((U[0, 0] + U[1, 0])/2))
		# cn.add_data(omega, t, Interval(U[1, 0] - U[0, 0]))

		#time.sleep(0.05)

	cn.contract()
	print(phi)
		
	# Graphics
	# beginDrawing()
	# fig = VIBesFigTube("Tube")
	# fig.set_properties(100, 100, 600, 300)
	# fig.add_tube(phi[0], "phi", "#376D7C[lightGray]")
	# fig.show(True)
	# endDrawing()
