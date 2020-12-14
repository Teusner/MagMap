import numpy as np
from vibes import *
import time
from itertools import cycle
from pyibex import *
from tubex_lib import *


"""
TODO :
 - docstring updating
 - supress some class variable which could be replaced as local variables
 - fix the loose rope detector with the 0.5 as threshold which is not clean
 - Adding Interval analysis to estimate the state of the robot
 - Showing the uncertainty on phi
 - Computing the trajectory of the magnetometer using Interval analysis
 - Adding the figure as class variable
 - Adding some figures : Mapping Coverage, Sled trajectory, Tank trajectory
 - Storage of the trajectories of the robot and the sled in a Trajectory object :
	y = TrajectoryVector(2)
	for t in np.arange(0., 10., 0.6):
		y.set([np.cos(t),np.sin(t)], t)
 - Better import of vibes
 - Better import of pyibex and of tubex_lib
 - Set the figure windows depending on the domain of the trajectories and refreshing it
 - Improvements :
	draw other vehicles with trajectories
	add groups in VIBesFigMap
	add some drawing functions like draw_vehicle or draw_triangle
"""


class Tank():
	def __init__(self, X=np.array([[0.], [0.], [0.], [0.]], dtype=np.float), L=3, h=1/20, show=True):
		"""
		Constructor of the vehicle.

		Parameters
		----------
		X : numpy.ndarray
			State vector of the vehicle [[x], [y], [theta]], x and y are in meters, theta is in radians
		U : numpy.ndarray
			Control vector of the vehicle [[u1], [u2]], u1 is for the left wheels and u2 for the right wheels
		M : numpy.ndarray
			Magnetometer position [[xm], [ym]], where xm is the abscissa and ym is the ordinate
		L : float
			The length of the rope
		IX : IntervalVector(2)
			Box enclosing the vehicle position in 2D in meters.
		"""
		# Integration step
		self.h = h

		# Measurement iterator
		self.gnss_measurement = cycle([k==0 for k in range(20)])

		# State representation
		self.X = X
		self.Y = self.g(np.array([[0], [0]]))
		self.M = X[:2] - L/2 * np.array([[np.cos(X[2, 0]+X[3, 0])], [np.sin(X[2, 0]+X[3, 0])]])

		# Rope variables
		self.L = L
		self.loose_rope = True

		# Path of the vehicle
		self.path = X[:2]
		self.anchor_point = self.X[:2] - 0.3*np.array([[np.cos(self.X[2, 0])], [np.sin(self.X[2, 0])]])

		# Trajectory
		self.t = 0
		self.robot_trajectory = TrajectoryVector(2)
		self.magnetometer_trajectory = TrajectoryVector(2)
		self.robot_trajectory.set(self.X[:2].flatten().tolist(), self.t)
		self.magnetometer_trajectory.set(self.M[:2].flatten().tolist(), self.t)
		
		# Initialize VIBes
		if show:
			self.__init_vibes__()

	def f(self, U):
		"""
		Dynamic function of the vehicle.

		Parameters
		----------
		X : numpy.ndarray
			State of the robot
		U : numpy.ndarray
			Control vector
		dx : float
			Derivative of x
		dy : float
			Derivative of y
		dtheta : float
			Derivative of theta
		dphi : float
			Derivative of phi

		Returns
		-------
		numpy.ndarray
			The derivative of the state vector with the current control vector.
		"""
		v = (U[0, 0] + U[1, 0])/2
		dx = v * np.cos(self.X[2, 0])
		dy = v * np.sin(self.X[2, 0])
		dtheta = U[1, 0] - U[0, 0]
		dphi = - v * np.sin(self.X[3, 0]) / self.L - dtheta
		return np.array([[dx], [dy], [dtheta], [dphi]])

	def g(self, U):
		"""
		Compute the observation vector which is the measurement of the vehicle's state.

		Parameters
		----------
		X : numpy.ndarray
			State of the robot
		U : numpy.ndarray
			Control vector
		p : numpy.ndarray
			Measured position of the robot p = [[xm], [ym]]
		theta : numpy.ndarray
			Measured angle of the robot
		v : numpy.ndarray
			Measured linear speed of the robot v = [[vx], [vy]]

		Returns
		-------
		numpy.ndarray
			The observation vector.
		"""

		# Uncertainties
		self.sigma_theta = 0.01
		self.sigma_gnss = 0.5
		self.sigma_v = 0.03

		if next(self.gnss_measurement):
			self.p = self.X[:2] + self.sigma_gnss * (np.random.random_sample((2, 1)) - 0.5)
		theta = self.X[2] + self.sigma_theta * (np.random.random() - 0.5)
		v = (U[0, 0] + U[1, 0])/2 * np.array([[np.cos(self.X[2, 0])], [np.sin(self.X[2, 0])]]) + self.sigma_v * (np.random.random_sample((2, 1)) - 0.5)
		return np.vstack((self.p, theta, v))

	def step(self, U):
		# Processing the new state
		self.t += self.h
		self.X += self.h * self.f(U)
		self.Y = self.g(U)

		# Processing the new anchor point position
		self.anchor_point = self.X[:2] - 0.3*np.array([[np.cos(self.X[2, 0])], [np.sin(self.X[2, 0])]])

		# Updating the path of the vehicle
		self.path = np.append(self.path, self.X[:2], axis=1)
		self.robot_trajectory.set(self.X[:2].flatten().tolist(), self.t)
		self.magnetometer_trajectory.set(self.M[:2].flatten().tolist(), self.t)

		# Rope processing
		rope = self.anchor_point - self.M
		v = (U[0, 0] + U[1, 0])/2

		# Loose rope testing
		if np.sign(v) * np.array([[np.cos(self.X[2, 0])], [np.sin(self.X[2, 0])]]).T @ rope > 0 and np.linalg.norm(rope) >= self.L - 0.5:
			self.loose_rope = False
			# Processing the new Magnetometer position
			R = np.array([[np.cos(self.X[2, 0]), - np.sin(self.X[2, 0])], [np.sin(self.X[2, 0]), np.cos(self.X[2, 0])]])
			self.M = self.X[:2] - self.L * R @ np.array([[np.cos(self.X[3, 0])], [np.sin(self.X[3, 0])]])
		else:
			self.loose_rope = True
			# Updating phi
			self.X[3, 0] = np.arctan2(rope[1, 0], rope[0, 0]) - self.X[2, 0]

	def __init_vibes__(self):
		"""
		Initialize VIBes software in order to represent the vehicle and the magnetometer.
		This create a new figure and 3 groups to draw the vehicle, his path and the
		magnetometer.
		"""
		# Initialize VIBes
		vibes.beginDrawing()
		vibes.newFigure("MagMap")
		vibes.axisLimits(-10, 10, -10, 10)
		vibes.axisEqual()
		vibes.axisLabels("x axis", "y axis", figure="MagMap")
		vibes.setFigureProperties({"x": 200, "y": 200, "width": 600, "height": 600})

		# Creating drawing groups
		vibes.newGroup("vehiclePath")
		vibes.newGroup("vehicle")
		vibes.newGroup("magnetometer")

		# beginDrawing()
		# self.figure = VIBesFigMap("MagMap_Trajectory")
		# self.figure.set_properties(100, 100, 600, 600)
		# self.figure.add_trajectory(self.robot_trajectory, "*x", 0, 1, "#2980b9")
		# self.figure.add_trajectory(self.magnetometer_trajectory, "*xm", 0, 1, "#8e44ad")
		# self.figure.show(0.01)


	def draw(self):
		"""
		Draw the scene with the vehicle, the magnetometer and the rope.
		The rope is drawn in red when it's loose and in green otherwise.
		The drawings are separated on three layers in order to be able 
		to erase some traces during the simulation.
		"""
		# Draw the vehicle
		vibes.clearGroup("vehicle")
		vibes.drawVehicle(self.X[0, 0], self.X[1, 0], self.X[2, 0]*180/np.pi, 1.0, color="black", figure="MagMap", group="vehicle")
		vibes.drawCircle(self.anchor_point[0, 0], self.anchor_point[1, 0], 0.08, color="black[grey]", group="vehicle")

		vibes.drawPoint(self.Y[0, 0], self.Y[1, 0], radius=3, color="#AAAAAA[#77116677]", group="estimatedState")
		
		# Draw the path
		vibes.drawLine(np.transpose(self.path[-2:]).tolist(), color="lightGray", group="vehiclePath")  

		# Drawing the magnetometer
		rope_color = "darkRed" if self.loose_rope else "darkGreen"
		vibes.clearGroup("magnetometer")
		vibes.drawLine([self.anchor_point[:, 0].tolist(), self.M[:, 0].tolist()], color=rope_color, group="magnetometer")
		vibes.drawCircle(self.M[0, 0], self.M[1, 0], 0.2, color="black[darkCyan]", group="magnetometer")

		# self.figure.draw_polygon([self.X[0, 0], self.X[1, 0], self.X[2, 0]*180/np.pi])
		# self.figure.draw_circle(self.anchor_point[0, 0], self.anchor_point[1, 0], 0.08, color="black[grey]")
		# self.figure.show(0.01)


if __name__=="__main__":
	t = Tank()
	h = 1/20
	for i in range (500):
		if i % 100 > 50:
			U = np.array([[1], [0.7]])
		else:
			U = np.array([[0.7], [1]])
		t.step(U)
		t.draw()
		time.sleep(0.05)
