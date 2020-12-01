import sys
from pathlib import Path
sys.path.append(str(Path().resolve()))

from tank import Tank
import numpy as np
import matplotlib.pyplot as plt
import time


class Controller:
	def __init__(self, t0=0, tf=128, h=0.1):
		self.t0, self.tf, self.h = 0, 128, 0.1
		self.T = np.arange(self.t0, self.tf, self.h)

		self.t = Tank()

		self.x = 20*np.sin(2*self.T/20)
		self.y = 20*np.sin(3*self.T/20)
		self.dx = 2*np.cos(2*self.T/20)
		self.dy = 3*np.cos(3*self.T/20)
		self.ddx = -0.05*4*np.sin(2*self.T/20)
		self.ddy = -0.05*9*np.sin(3*2*self.T/20)

		self.Su = np.empty(shape=(0, 2))
		self.Sg = np.empty(shape=(0, 5))

	def show_trajectory(self):
		plt.plot(self.x, self.y)
		plt.quiver(self.x, self.y, self.dx, self.dy, color="green")
		plt.show()

	def process_data(self):
		z = 1.
		for xi, yi, dxi, dyi, ddxi, ddyi in zip(self.x, self.y, self.dx, self.dy, self.ddx, self.ddy):
			xt, yt, theta = self.t.X[0, 0], self.t.X[1, 0], self.t.X[2, 0]
			A = np.array([[np.cos(theta), -z*np.sin(theta)], [np.sin(theta), z*np.cos(theta)]])
			v = np.array([[xi-xt + 2*(dxi - z*np.cos(theta)) + ddxi], [yi - yt + 2*(dyi - z*np.sin(theta)) + ddyi]])
			c = np.linalg.solve(A, v)

			z += self.h*c[0, 0]
			
			A = np.array([[0.5, 0.5], [-1, 1]])
			B = np.array([[z], c[1]])
			U = np.linalg.solve(A, B)
			self.Su = np.vstack((self.Su, U.flatten()))
			self.t.step(U)
			self.Sg = np.vstack((self.Sg, self.t.g(U).flatten()))
			self.t.draw()
		return self.Su, self.Sg


if __name__ == "__main__":
	c = Controller()
	Su, Sg = c.process_data()
	np.savetxt("control.txt", Su)
	np.savetxt("observations.txt", Sg)