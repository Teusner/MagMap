import numpy as np
from vibes import *
from pyibex import *
from pyibex.thickset import *
from tubex_lib import *

if __name__ == "__main__":

##### LOADING DATA
	# Loading data (heading)
	Sh = np.load("./data/angle.npy")
	t_h, heading = Sh[:, 0].T, Sh[:, 1]
	accuracy_heading = 0.1

	# Loading data (control)
	Su = np.load("./data/control.npy")
	t_u, control = Su[:, 0].T, Su[:, 1:]
	accuracy_control = 0.5

	# Loading data (velocity)
	Sv = np.load("./data/velocity.npy")
	t_v, velocity = Sv[:, 0].T, Sv[:, 1:]
	accuracy_velocity = 3

	# Loading data (gnss)
	Sg = np.load("./data/gnss.npy")
	t_g, gnss = Sg[:, 0].T, Sg[:, 1:]
	accuracy_gnss = 1

##### CREATING TUBES
	# Getting time in data
	t0, tf, h = t_h[0], t_h[-1], t_h[1]-t_h[0]

	# Creating tubes x and v
	tdomain = Interval(t0, tf)
	x = TubeVector(tdomain, h, IntervalVector(4))
	v = TubeVector(tdomain, h, IntervalVector(4))
	u = TubeVector(tdomain, h, IntervalVector(2))

##### UPDATING TUBES WITH DATA
	# Heading trajectory
	trajectory_h = Trajectory(dict(zip(t_h, heading)))
	trajectory_h.truncate_tdomain(tdomain)

	# Updating heading x[2] with the Trajectory
	x[2] &= trajectory_h
	x[2].inflate(accuracy_heading)

	# Control trajectory
	trajectory_u = TrajectoryVector(dict(zip(t_u, control.tolist())))
	trajectory_u.truncate_tdomain(tdomain)

	# Updating control u with the Trajectory
	u &= trajectory_u
	u.inflate(accuracy_control)

	# Velocity trajectory
	trajectory_v = TrajectoryVector(dict(zip(t_v, velocity.tolist())))
	trajectory_v.truncate_tdomain(tdomain)

	# Updating velocity v with the Trajectory
	v[0] &= trajectory_v[0]
	v[1] &= trajectory_v[1]
	v[0].inflate(accuracy_velocity)
	v[1].inflate(accuracy_velocity)

	# Setting the initial condition on phi
	x[3].set(Interval(-np.pi/4, np.pi/4), 0.)

##### Adding contrains
	# Contractor Network
	cn = ContractorNetwork()

	# Deriv constrain between x and v
	ctc_deriv = CtcDeriv()
	cn.add(ctc_deriv, [x, v])

	# f constrain
	ctc_f1 = CtcFunction(Function("v0", "vx", "vy", "x2", "v0-sqrt(vx^2+vy^2)*cos(x2)"))
	ctc_f2 = CtcFunction(Function("v1", "vx", "vy", "x2", "v1-sqrt(vx^2+vy^2)*sin(x2)"))
	ctc_f3 = CtcFunction(Function("v2", "u1", "u2", "v2+u1-u2"))
	ctc_f4 = Function("x", "-sin(x)")

	cn.add(ctc_f1, [v[0], v[0], v[1], x[2]])
	cn.add(ctc_f2, [v[1], v[0], v[1], x[2]])
	cn.add(ctc_f3, [v[2], u[0], u[1]])
	ctc_lohner = CtcLohner(ctc_f4)
	cn.add(ctc_lohner, [x[3]])

	# GNSS constrains
	ctc_eval = CtcEval()
	ctc_eval.enable_time_propag(False)
	for i in range(len(t_g)):
		p = IntervalVector([Interval(gnss[i, 0]), Interval(gnss[i, 1])])
		p.inflate(accuracy_gnss)
		cn.add(ctc_eval, [Interval(t_g[i]), p[0], x[0], v[0]])
		cn.add(ctc_eval, [Interval(t_g[i]), p[1], x[1], v[1]])

	# Contract the 
	cn.contract(True)
	
	beginDrawing()

	# Showing trajectory
	beginDrawing()
	fig_map = VIBesFigMap("Saturne Trajectory")
	fig_map.set_properties(100, 100, 600, 300)
	fig_map.smooth_tube_drawing(True)
	fig_map.add_tube(x, "x*", 0, 1)
	fig_map.axis_limits(-2.5,2.5,-0.1,0.1, True)
	fig_map.show()

	# Showing U1
	fig_dist = VIBesFigTube("Heading and Control Tubes")
	fig_dist.set_properties(100, 430, 600, 300)
	fig_dist.add_tube(x[2], "Heading", "#2c3e50[#2c3e50]")
	fig_dist.add_tube(u[0], "U1", "#2980b9[#2980b9]")
	fig_dist.add_tube(u[1], "U2", "#c0392b[#c0392b]")
	fig_dist.add_tube(x[3], "x3", "red")
	fig_dist.show()

	endDrawing()