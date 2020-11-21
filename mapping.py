import numpy as np
from vibes import *
from pyibex import *
from tubex_lib import *

from tank import Tank


class Mapping:
	def __init__(self, t0, tf, h):
		self.h = h

		# Trajectory Tubes
		self.tdomain = Interval(t0, tf+h)
		self.x = TubeVector(self.tdomain, self.h, IntervalVector(2))
		self.v = TubeVector(self.tdomain, self.h, IntervalVector(2, Interval(0)))
		self.v.inflate(10)
		self.a = TubeVector(self.tdomain, self.h, IntervalVector(2, Interval(0)))
		self.a.inflate(10)

		# Contractor Network
		self.cn = ContractorNetwork()
		ctc_deriv1, ctc_deriv2 = CtcDeriv(), CtcDeriv()
		self.cn.add(ctc_deriv1, [self.x, self.v])
		self.cn.add(ctc_deriv2, [self.v, self.a])
		self.cn.contract()

	def add_position(self, t, position, accuracy):
		# Creating a contractor
		ctc_eval = CtcEval()
		ctc_eval.enable_time_propag(False)

		n = len(t)
		for i in range(n):
			p = IntervalVector([Interval(position[i, 0]), Interval(position[i, 1])])
			p.inflate(accuracy)
			ctc_eval.contract(Interval(t[i], t[i]+h), p, self.x, self.v)
		ctc.deriv.contract(self.x, self.v)
		
	def add_velocity(self, t, velocity, accuracy):
		# Creating a contractor
		ctc_eval = CtcEval()
		ctc_eval.enable_time_propag(False)

		n = len(t)
		for i in range(n):
			V = IntervalVector([Interval(velocity[i, 0]), Interval(velocity[i, 1])])
			V.inflate(accuracy)
			ctc_eval.contract(Interval(t[i], t[i]+h), V, self.v, self.a)
		ctc.deriv.contract(self.v, self.a)
		ctc.deriv.contract(self.x, self.v)


if __name__ == "__main__":
	# Time
	t0, tf, h = 0, 10.1, 1/20
	T = np.arange(t0, tf, h)

	m = Mapping(t0, tf, h)
	L = 3
	tank = Tank(L=L, h=h, show=False)

	# Lu = [U]
	Lu = []

	# Lg = [p, theta, v]
	Lg = []

	for t in T:
		# Command generating
		if t % 2 > 1:
			U = np.array([[1], [0.7]])
		else:
			U = np.array([[0.7], [1]])

		tank.step(U)

		# Data storage
		Lu.append(U.flatten())
		Lg.append(tank.g(U).flatten())
	
	Lu, Lg = np.asarray(Lu), np.asarray(Lg)
	
	# Adding positions in m
	accuracy_p = 0.5
	gnss = np.vstack((np.unique(Lg[:,0]), np.unique(Lg[:,1]))).T
	t_gnss = T[np.nonzero(np.r_[1, np.diff(Lg[:,0])[:-1]])]
	m.add_position(t_gnss, gnss, accuracy_p)

	# Adding velocities in m
	accuracy_v = 0.03
	m.add_velocity(T, Lg[:, [3,4]], accuracy_v)

	# Contracting
	m.cn.contract()

	beginDrawing()
	fig_map = VIBesFigMap("Saturne")
	fig_map.set_properties(100, 100, 600, 300)
	fig_map.smooth_tube_drawing(True)
	fig_map.add_tube(m.x, "x*", 0, 1)
	fig_map.axis_limits(-2.5,2.5,-0.1,0.1, True)
	fig_map.show()
	endDrawing()
