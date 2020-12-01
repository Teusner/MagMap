import numpy as np
from vibes import *
from pyibex import *
from pyibex.thickset import *
from tubex_lib import *

from tank import Tank
from differential_integration import *
import time as time


class Mapping:
	def __init__(self, t0, tf, h, L):
		# Variables
		self.t0 = t0
		self.tf = tf
		self.h = h
		self.L = L
		self.measurement_range = 5.0

		# Trajectory Tubes
		self.tdomain = Interval(self.t0, self.tf+h)
		self.x = TubeVector(self.tdomain, self.h, IntervalVector(2))
		self.v = TubeVector(self.tdomain, self.h, IntervalVector(2, Interval(0)))
		self.v.inflate(10)
		self.a = TubeVector(self.tdomain, self.h, IntervalVector(2, Interval(0)))
		self.a.inflate(10)

		# Tube enclosing the magnetometer
		self.x_m = TubeVector(self.tdomain, self.h, IntervalVector(2))
		self.phi = Tube(self.tdomain, self.h, Interval(0))
		self.phi.inflate(np.pi/2)

		# Contractor Network
		self.cn = ContractorNetwork()
		ctc_deriv = CtcDeriv()
		self.cn.add(ctc_deriv, [self.x, self.v])
		self.cn.add(ctc_deriv, [self.v, self.a])
		
		# Magnetometer contractors
		ctc_polar = CtcPolar()
		self.ctc_f = CtcFunction(Function("x1", "x2", "x3", "x1-x2+x3"))

		# Temporary variable
		self.d = TubeVector(self.tdomain, self.h, IntervalVector(2))

		# Constrains
		self.cn.add(ctc_polar, [self.d[0], self.d[1], Interval(self.L), self.phi])
		self.cn.add(self.ctc_f, [self.x_m[0], self.x[0], self.d[0]])
		self.cn.add(self.ctc_f, [self.x_m[1], self.x[1], self.d[1]])

	def add_position(self, t, position, accuracy):
		# Creating a contractor
		ctc_eval = CtcEval()
		ctc_eval.enable_time_propag(False)

		n = len(t)
		for i in range(n):
			p = IntervalVector([Interval(position[i, 0]), Interval(position[i, 1])])
			p.inflate(accuracy)
			ctc_eval.contract(Interval(t[i]), p, self.x, self.v)
		
	def add_velocity(self, t, velocity, accuracy):
		# Creating a contractor
		ctc_eval = CtcEval()
		ctc_eval.enable_time_propag(False)

		n = len(t)
		for i in range(n):
			V = IntervalVector([Interval(velocity[i, 0]), Interval(velocity[i, 1])])
			V.inflate(accuracy)
			ctc_eval.contract(Interval(t[i]), V, self.v, self.a)
	
	def add_control(self, t, U):
		# Contracting phi
		alpha = 0.01		
		n = len(U)
		for i in range(n):
			v, theta = Interval((U[i, 0]+U[i, 1])/2), Interval(U[i, 1] - U[i, 0])
			I, out = integrate(self.phi(i), alpha, v, theta, self.h)
			lb, ub = I.lb(), I.ub()
			self.phi.set(Interval(lb, ub), i+1)

	def process_coverage(self):
		# The map
		self.m_map = IntervalVector(2, Interval(-20, 20))
		y = Interval(-oo, 0)
		f_dist = Function("x1", "x2", "p1", "p2", "p3", "(x1-p1)^2+(x2-p2)^2-p3^2")

		# Initialisation
		p = IntervalVector(3, Interval())
		p[0] = self.x_m[0](0)
		p[1] = self.x_m[1](0)
		p[2] = Interval(self.measurement_range)
		S_dist = ThickSep_from_function(f_dist, p, y)

		n = int((self.tf-self.t0)/self.h) + 1
		for i in range(1, n):
			p = IntervalVector(3, Interval())
			p[0] = self.x_m[0](i)
			p[1] = self.x_m[1](i)
			p[2] = Interval(self.measurement_range)
			S_dist = S_dist | ThickSep_from_function(f_dist, p, y)
		return ThickPaving(self.m_map, ThickTest_from_ThickSep(S_dist), 0.1)


if __name__ == "__main__":
	# Time
	t0, tf, h = 0, 10.15, 1/20
	T = np.arange(t0, tf, h)

	L = 3
	m = Mapping(t0, tf, h, L)
	tank = Tank(L=L, h=h, show=False)

	# Lu = [U], Lg = [p, theta, v]
	Lu, Lg = [], []

	for t in T:
		# Command generating
		U = np.array([[1.0], [1.0]])

		tank.step(U)

		# Data storage
		Lu.append(U.flatten())
		Lg.append(tank.g(U).flatten())
	
	Lu, Lg = np.asarray(Lu), np.asarray(Lg)
	
	# Benchmark
	t_bench = time.time()

	# Adding positions in m
	accuracy_p = 0.5
	gnss = np.vstack((np.unique(Lg[:,0]), np.unique(Lg[:,1]))).T
	t_gnss = T[np.nonzero(np.r_[1, np.diff(Lg[:,0])[:-1]])]
	m.add_position(t_gnss, gnss, accuracy_p)

	# Adding velocities in m
	accuracy_v = 0.03
	m.add_velocity(T, Lg[:, [3,4]], accuracy_v)

	# Adding control
	m.add_control(T, Lu)

	print("Constrains adding time {:.4} s".format(time.time()-t_bench))

	# Contracting
	m.cn.contract(True)

	beginDrawing()
	fig_map = VIBesFigMap("Saturne")
	fig_map.set_properties(100, 100, 600, 300)
	fig_map.smooth_tube_drawing(True)
	fig_map.add_tube(m.x, "x*", 0, 1)
	fig_map.axis_limits(-2.5,2.5,-0.1,0.1, True)
	fig_map.show()

	fig_map2 = VIBesFigMap("Magnetometer")
	fig_map2.set_properties(700, 100, 600, 300)
	fig_map2.smooth_tube_drawing(True)
	fig_map2.add_tube(m.x_m, "x*", 0, 1)
	fig_map2.axis_limits(-2.5,2.5,-0.1,0.1, True)
	fig_map2.show()

	fig_dist = VIBesFigTube("Phi angle")
	fig_dist.set_properties(100, 430, 600, 300)
	fig_dist.add_tube(m.phi, "y*")
	fig_dist.show()

	endDrawing()

	# Coverage computation
	vibes.beginDrawing()
	P = m.process_coverage()
	vibes.setFigureSize(m.m_map.max_diam(), m.m_map.max_diam())
	vibes.newFigure("Mapping Coverage")
	vibes.setFigureProperties({"x": 700, "y": 430, "width": 600, "height": 600})
	P.visit(ToVibes(figureName="Mapping Coverage", color_map=My_CMap))
	vibes.endDrawing()