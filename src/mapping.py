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
		# Physical variables
		self.L = L
		self.measurement_range = 5.0

		# Temporal variables
		self.t0 = t0
		self.tf = tf
		self.h = h
		self.tdomain = Interval(self.t0, self.tf)
		
		# Contractor Network
		self.cn = ContractorNetwork()

		# Creating tubes
		self.create_tubes()
		#self.create_tubes_magnetometer()

	def create_tubes(self):
		# State tubes
		self.x = TubeVector(self.tdomain, self.h, IntervalVector(3))
		self.v = TubeVector(self.tdomain, self.h, IntervalVector(3, Interval(-100, 100)))
		self.cn.add(ctc.deriv, [self.x, self.v])

		# Control tube
		self.u = TubeVector(self.tdomain, self.h, IntervalVector(2))
	
	def create_tubes_magnetometer(self):
		# Tube enclosing the magnetometer
		self.x_m = TubeVector(self.tdomain, self.h, IntervalVector(2))
		self.phi = Tube(self.tdomain, self.h, Interval(-np.pi/2, np.pi/2))

		# Magnetometer contractors
		ctc_polar = CtcPolar()
		self.ctc_f = CtcFunction(Function("x1", "x2", "x3", "x1-x2+x3"))

		# Temporary variable
		self.d = TubeVector(self.tdomain, self.h, IntervalVector(2))
		self.theta = Tube(self.tdomain, self.h, Interval())

		# Constrains
		self.cn.add(self.ctc_f, [self.theta, self.phi, self.x[2]])
		self.cn.add(ctc_polar, [self.d[0], self.d[1], Interval(self.L), self.theta])
		self.cn.add(self.ctc_f, [self.x_m[0], self.x[0], self.d[0]])
		self.cn.add(self.ctc_f, [self.x_m[1], self.x[1], self.d[1]])

	def add_f(self):
		ctc_f1 = CtcFunction(Function("v0", "vx", "vy", "x2", "v0-sqrt(vx^2+vy^2)*cos(x2)"))
		ctc_f2 = CtcFunction(Function("v1", "vx", "vy", "x2", "v1-sqrt(vx^2+vy^2)*sin(x2)"))
		ctc_f3 = CtcFunction(Function("v2", "u1", "u2", "v2+u1-u2"))

		self.cn.add(ctc_f1, [self.v[0], self.v[0], self.v[1], self.x[2]])
		self.cn.add(ctc_f2, [self.v[1], self.v[0], self.v[1], self.x[2]])
		self.cn.add(ctc_f3, [self.v[2], self.u[0], self.u[1]])

	def add_position(self, t, position, accuracy):
		# Creating a contractor
		ctc_eval = CtcEval()
		ctc_eval.enable_time_propag(False)

		n = len(t)
		for i in range(n):
			p = IntervalVector([Interval(position[i, 0]), Interval(position[i, 1])])
			p.inflate(accuracy)
			self.cn.add(ctc_eval, [Interval(t[i]), p[0], self.x[0], self.v[0]])
			self.cn.add(ctc_eval, [Interval(t[i]), p[1], self.x[1], self.v[1]])

	def add_heading(self, t, heading, accuracy):
		trajectory_h = Trajectory(dict(zip(t, heading)))
		trajectory_h.truncate_tdomain(self.tdomain)
		self.h_t = Tube(trajectory_h, self.h)
		self.h_t.inflate(accuracy)

		ctc_equal = CtcFunction(Function("a", "b", "a-b"))

		# Storage of the velocity in the tube
		self.cn.add(ctc_equal, [self.x[2], self.h_t])
			
	def add_velocity(self, t, velocity, accuracy):
		trajectory_v = TrajectoryVector(dict(zip(t, velocity.tolist())))
		trajectory_v.truncate_tdomain(self.tdomain)

		self.v_t = TubeVector(trajectory_v, self.h)
		self.v_t.inflate(accuracy)

		# ctc_equal = CtcFunction(Function("a", "b", "a-b"))

		# Storage of the velocity in the tube
		# self.cn.add(ctc_equal, [self.v[0], self.v_t[0]])
		# self.cn.add(ctc_equal, [self.v[1], self.v_t[1]])

		self.v[0] &= trajectory_v[0]
		self.v[1] &= trajectory_v[1]

		# Inflating these tube with accuracy
		self.v[0].inflate(accuracy)
		self.v[1].inflate(accuracy)

	def add_control(self, t, U):
		trajectory_u = TrajectoryVector(dict(zip(t, U.tolist())))
		trajectory_u.truncate_tdomain(self.tdomain)
		self.u &= trajectory_u

	def process_coverage(self):
		# The map
		self.m_map = IntervalVector(2, Interval(-35, 35))
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

	# Loading control and measurements data
	Su = np.load("./data/control.npy")
	Sg = np.load("./data/gnss.npy")
	St = np.load("./data/angle.npy")
	Sv = np.load("./data/velocity.npy")

	# Getting time parameters
	T = Su[:, 0].T
	t0, tf, h = T[0], T[-1], T[1]-T[0]

	L = 3
	m = Mapping(t0, tf, h, L)
	tank = Tank(L=L, h=h, show=False)
	
	# Benchmark
	t_bench = time.time()

	# Adding evolution function
	m.add_f()

	# Adding positions in m
	accuracy_p = 0.5
	m.add_position(Sg[:, 0].T, Sg[:, 1:], accuracy_p)

	# Adding velocities in m
	accuracy_v = 0.03
	m.add_velocity(Sv[:, 0].T, Sv[:, 1:], 65*accuracy_v)

	# Adding heading in m
	accuracy_t = 0.01
	m.add_heading(St[:, 0].T, St[:, 1], 26*accuracy_t)

	# Adding control
	m.add_control(Su[:, 0].T, Su[:, 1:])

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

	# fig_map2 = VIBesFigMap("Magnetometer")
	# fig_map2.set_properties(700, 100, 600, 300)
	# fig_map2.smooth_tube_drawing(True)
	# fig_map2.add_tube(m.x_m, "x*", 0, 1)
	# fig_map2.axis_limits(-2.5,2.5,-0.1,0.1, True)
	# fig_map2.show()

	fig_dist = VIBesFigTube("Control")
	fig_dist.set_properties(100, 430, 600, 300)
	fig_dist.add_tube(m.u[0], "u0", "blue")
	fig_dist.add_tube(m.u[1], "u1", "red")
	fig_dist.show()

	endDrawing()

	# Coverage computation
	# vibes.beginDrawing()
	# P = m.process_coverage()
	# vibes.setFigureSize(m.m_map.max_diam(), m.m_map.max_diam())
	# vibes.newFigure("Mapping Coverage")
	# vibes.setFigureProperties({"x": 700, "y": 430, "width": 600, "height": 600})
	# P.visit(ToVibes(figureName="Mapping Coverage", color_map=My_CMap))
	# vibes.endDrawing()