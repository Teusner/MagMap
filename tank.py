import numpy as np
from vibes import vibes
from kalman import kalman
from scipy.linalg import sqrtm
import time

"""
TODO :
 - docstring updating
 - supress some class variable which could be replaced as local variables
 - fix the loose rope detector with the 0.5 as threshold which is not clean
 - tune the kalman filter
 - supress set_cmd(U) and give U for each step
 - fix the kalmann filter with B and U
 - Showing the uncertainty on phi
 - Computing the trajectory of the magnetometer using Interval analysis
"""


class Tank():
	def __init__(self, X=np.array([[0.], [0.], [0.], [0.]], dtype=np.float), width=0.6, wheel_radius=0.5, alpha=1, L=3, h=1/20):
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
		width : float
			Width of the robot between the two wheels in meters
		wheel_radius : float
			Wheel's radius in meters
		alpha : float
			Coefficient between the motor's control and the wheel's spinning rate
		"""
		self.X = X
		self.M = X[:2] - L/2 * np.array([[np.cos(X[2, 0])], [np.sin(X[2, 0])]])
		self.L = L
		self.width = width
		self.wheel_radius = wheel_radius
		self.alpha = alpha
		self.path = X[:2]
		self.anchor_point = self.X[:2] - 0.3*np.array([[np.cos(self.X[2, 0])], [np.sin(self.X[2, 0])]])
		self.loose_rope = True

		# Integration step
		self.h = h

		# State estimation
		self.Xhat = np.zeros(X.shape)
		self.Γx = 1000 * np.eye(4)

		self.__init_vibes__()
		

	def f(self, U):
		"""
		Dynamic function of the vehicle.

		Parameters
		----------
		dx : float
			Derivative of x
		dy : float
			Derivative of y
		dtheta : float
			Derivative of theta

		Returns
		-------
		numpy.ndarray
			The derivative of the state vector with the current control vector.
		"""
		v = (U[0, 0] + U[1, 0])/2
		dx = self.wheel_radius * self.alpha * v * np.cos(self.X[2, 0])
		dy = self.wheel_radius * self.alpha * v * np.sin(self.X[2, 0])
		dtheta = self.wheel_radius*self.alpha/self.width*(U[1, 0] - U[0, 0])
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

		Returns
		-------
		numpy.ndarray
			The observation vector.
		"""
		rn = np.array([[1], [1], [0.05]]) * np.random.randn(3, 1)
		return self.X[:3] + rn

	def jacobian(self, U):
		"""
		Compute the Jacobian matrix of the system df(X, U)/dX.

		Parameters
		----------
		X : numpy.ndarray
			State of the robot
		U : numpy.ndarray
			Control vector

		Returns
		-------
		numpy.ndarray
			The Jacobian matrix of the system.
		"""
		v = (U[0, 0] + U[1, 0])/2
		J = np.zeros((4, 4))
		J[0, 2] = - self.wheel_radius * self.alpha * v * np.sin(self.X[2, 0])
		J[1, 2] = self.wheel_radius * self.alpha * v * np.cos(self.X[2, 0])
		J[2, 3] = - v * np.cos(self.X[3, 0]) / self.L
		return J

	def state_estimation(self, U):
		"""
		Vehicle state estimation using a kalman filter.
		"""
		A = np.eye(4) + h * self.jacobian(U)
		B = np.zeros((4, 2))
		C = np.eye(4)[:3, :]
		Y = self.g(U)
		U = self.h * B @ U
		Γα = self.h * np.diag((1, 1, 0.05, 0.1))
		Γβ = np.diag((1, 1, 0.1))
		self.Xhat, self.Γx = kalman(self.Xhat, self.Γx, U, Y, Γα, Γβ, A, C)

	def step(self, U):
		# Processing the new state
		self.X += self.h * self.f(U)

		# State estimation
		self.state_estimation(U)

		# Processing the new anchor point position
		self.anchor_point = self.X[:2] - 0.3*np.array([[np.cos(self.X[2, 0])], [np.sin(self.X[2, 0])]])

		# Updating the path of the vehicle
		self.path = np.append(self.path, self.X[:2], axis=1)

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
		vibes.newGroup("estimatedState")
		vibes.newGroup("vehiclePath")
		vibes.newGroup("vehicle")
		vibes.newGroup("magnetometer")

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

		# Draw the estimated state
		η = 0.9
		if np.linalg.norm(self.Γx[:2, :2]) == 0:
			self.Γx = self.Γx + 0.01 * np.eye(4)
		
		A = sqrtm(- 2 * np.log(1 - η) * self.Γx[:2, :2])
		_, v = np.linalg.eig(A)
		f1, f2 = A @ v[:, 0], A @ v[:, 1]
		a, b, α = 2 * np.linalg.norm(f1), 2 * np.linalg.norm(f2), np.arctan2(v[1, 0], v[0, 0]) * 180 / 3.14

		vibes.clearGroup("estimatedState")
		vibes.drawEllipse(self.X[0, 0], self.X[1, 0], a, b, α, color="#c7ecee[#dff9fb]", group="estimatedState")
		vibes.drawPoint(self.Xhat[0, 0], self.Xhat[1, 0], radius=1, color="#eb4d4b") #eb4d4b #95afc0

		# Draw the estimated rope angle
		p = self.anchor_point - 0.5 * np.array([[np.cos(self.Xhat[2, 0] + self.Xhat[3, 0])], [np.sin(self.Xhat[2, 0] + self.Xhat[3, 0])]])
		vibes.drawPoint(p[0, 0], p[1, 0], radius=3, color="red[red]", group="estimatedState")
		#vibes.drawPie(self.anchor_point.tolist(), [0.5, 0.52], [np.pi - self.Γx[3, 3]/2, np.pi + self.Γx[3, 3]/2], color="purple", use_radian=True, group="estimatedState")

		# Draw the path
		vibes.drawLine(np.transpose(self.path[-2:]).tolist(), color="lightGray", group="vehiclePath")  

		# Drawing the magnetometer
		rope_color = "darkRed" if self.loose_rope else "darkGreen"
		vibes.clearGroup("magnetometer")
		vibes.drawLine([self.anchor_point[:, 0].tolist(), self.M[:, 0].tolist()], color=rope_color, group="magnetometer")
		vibes.drawCircle(self.M[0, 0], self.M[1, 0], 0.2, color="black[darkCyan]", group="magnetometer")


if __name__=="__main__":
	t = Tank(alpha=3)
	h = 1/20
	for i in range (1000):
		if i % 70 > 35:
			U = np.array([[1], [0.7]])
		else:
			U = np.array([[0.7], [1]])
		t.step(U)
		t.draw()
		time.sleep(0.05)
