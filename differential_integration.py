import numpy as np
from pyibex import *
from tubex_lib import *
from vibes import vibes


def sub_interval(phi, n):
	if n <= 1:
		return [phi]
	A = []
	lb, d = phi.lb(), phi.diam()
	for k in range(n):
		A.append(Interval(lb + k * d / n, lb + (k+1) * d / n))
	return A

def union(I):
	A = Interval.EMPTY_SET 
	for i in I:
		A = A | i
	return A

def integrate_box(phi, v, dtheta, h):
	dphi = Tube(phi.tdomain(), h, Interval())
	ctc_f = CtcFunction(Function("phi", "dphi", "v", "dtheta", "dphi+v*sin(phi)-dtheta"))
	
	cn = ContractorNetwork()
	cn.add(ctc_f, [phi(0), dphi, v, dtheta])
	cn.add(ctc.deriv, [phi, dphi])
	cn.contract()
	
	return phi(1)

def integrate(phi, alpha, v, theta, h):
	out = []
	tdomain = Interval(0, h)
	n = int(np.sqrt(phi.diam()//alpha + 1))
	splited = sub_interval(phi, n)

	for interval in splited:
		phi = Tube(tdomain, h/2, Interval())
		phi.set(interval, Interval(0, h/2))
		out.append(integrate_box(phi, v, theta, h/2))
	
	return union(out), out


if __name__ == "__main__":
	h = 1/20
	phi = Interval(-np.pi/2, np.pi/2)

	alpha = 0.01
	deltat = 300

	vibes.beginDrawing()
	vibes.newFigure("MagMap")
	vibes.setFigureProperties({"x": 200, "y": 200, "width": 600, "height": 600})
	vibes.axisLimits(-h, deltat*h, -np.pi/2, np.pi/2)

	for k in range(deltat):
		v = Interval(1.0)
		theta = Interval(0.4*np.sin(k/20))
		I, out = integrate(phi, alpha, v, theta, h)
		lb, ub = I.lb(), I.ub()
		vibes.drawBox(k*h, (k+1)*h, lb, ub, "darkCyan[darkCyan]")
		phi = Interval(lb, ub)
	