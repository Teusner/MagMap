import numpy as np
from vibes import *
from pyibex import *
from pyibex.thickset import *
from tubex_lib import *

if __name__ == "__main__":
	# Loading data
	Su = np.load("./data/control.npy")
	t, u = Su[:, 0].T, Su[:, 1:]

	# Getting time in data
	t0, tf, h = t[0], t[-1], t[1]-t[0]


	#### FIRST METHOD : using CtcEval
	# Creating tubes
	tdomain = Interval(t0, tf)
	U = TubeVector(tdomain, h, IntervalVector(2))
	dU = TubeVector(tdomain, h, IntervalVector(2))
	dx3 = Tube(tdomain, h, Interval())
	x3 = Tube(tdomain, h, Interval())
	x3.set(Interval(-0.5, 0.5), 0)

	trajectory_u = TrajectoryVector(dict(zip(t, u.tolist())))
	trajectory_u.truncate_tdomain(tdomain)
	U = TubeVector(trajectory_u, h)

	print(U.tdomain)

	# Contractor Network
	cn = ContractorNetwork()

	# Deriv constrain between U and dU
	ctc_deriv = CtcDeriv()
	cn.add(ctc_deriv, [U, dU])
	cn.add(ctc_deriv, [x3, dx3])
	
	ctc_f3 = CtcFunction(Function("v2", "u1", "u2", "v2+u1-u2"))
	cn.add(ctc_f3, [dx3, U[0], U[1]])


	# Contract the 
	cn.contract(True)
	
	beginDrawing()

	# Showing U1
	fig_dist = VIBesFigTube("Control using CtcEval")
	fig_dist.set_properties(100, 430, 600, 300)
	fig_dist.add_tube(U[0], "U0")
	fig_dist.add_tube(U[1], "U1")
	fig_dist.add_tube(x3, "x3")
	fig_dist.show()

	endDrawing()