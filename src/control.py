from tank import Tank
import numpy as np
import matplotlib.pyplot as plt
import time


class Controller:
	def __init__(self, t0=0, tf=128, h=0.1):
		self.t0, self.tf, self.h = t0, tf, h
		self.T = np.arange(self.t0, self.tf, self.h)

		self.t = Tank()

		self.x = 20*np.sin(2*self.T/20)
		self.y = 20*np.sin(3*self.T/20)
		self.dx = 2*np.cos(2*self.T/20)
		self.dy = 3*np.cos(3*self.T/20)
		self.ddx = -0.05*4*np.sin(2*self.T/20)
		self.ddy = -0.05*9*np.sin(3*2*self.T/20)

		self.Su = []
		self.Sg = []
		self.Sv = []
		self.St = []

	def show_trajectory(self):
		plt.plot(self.x, self.y)
		plt.quiver(self.x, self.y, self.dx, self.dy, color="green")
		plt.show()

	def process_data(self):
		z = 1.
		for ti, xi, yi, dxi, dyi, ddxi, ddyi in zip(self.T, self.x, self.y, self.dx, self.dy, self.ddx, self.ddy):
			xt, yt, theta = self.t.X[0, 0], self.t.X[1, 0], self.t.X[2, 0]
			A = np.array([[np.cos(theta), -z*np.sin(theta)], [np.sin(theta), z*np.cos(theta)]])
			v = np.array([[xi-xt + 2*(dxi - z*np.cos(theta)) + ddxi], [yi - yt + 2*(dyi - z*np.sin(theta)) + ddyi]])
			c = np.linalg.solve(A, v)

			z += self.h*c[0, 0]
			
			A = np.array([[0.5, 0.5], [-1, 1]])
			B = np.array([[z], c[1]])
			U = np.linalg.solve(A, B)
			self.Su.append(np.hstack((np.array([ti]), U.flatten())))
			self.t.step(U)

			g = self.t.g(U)

			self.St.append(np.hstack((np.array([ti]), g[2, 0])))
			self.Sv.append(np.hstack((np.array([ti]), g[3:, 0])))

			if len(self.Sg) == 0:
				self.Sg.append(np.hstack((np.array([ti]),g[:2, 0])))
			else:
				if np.linalg.norm(self.Sg[-1][1:] - g[:2, 0]) > 0.1:
					self.Sg.append(np.hstack((np.array([ti]), g[:2, 0])))
			self.t.draw()
		return np.asarray(self.Su), np.asarray(self.Sg), np.asarray(self.St), np.asarray(self.Sv)


if __name__ == "__main__":
	c = Controller(tf=10)
	Su, Sg, St, Sv = c.process_data()
	
	# Saving data
	np.save("./data/control.npy", Su)
	np.save("./data/gnss.npy", Sg)
	np.save("./data/angle.npy", St)
	np.save("./data/velocity.npy", Sv)
