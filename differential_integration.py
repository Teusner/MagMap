import numpy as np
from pyibex import *
from tubex_lib import *
from vibes import vibes


def sub_interval(phi, n):
	if n <= 1:
		return phi
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

def integrate_step(phi, v, dtheta):
	dphi = Tube(phi.tdomain(), phi.tdomain().ub(), Interval())
	ctc_f = CtcFunction(Function("phi", "dphi", "v", "dtheta", "dphi+v*sin(phi)-dtheta"))
	
	cn = ContractorNetwork()
	cn.add(ctc_f, [phi, dphi, v, dtheta])
	cn.add(ctc.deriv, [phi, dphi])
	cn.contract()

	return phi

def integrate(phi, n, h):
	out = []
	tdomain = Interval(0, h)
	splited = sub_interval(phi, n)

	v, dtheta = Interval(1.0), Interval(0.)

	for interval in splited:
		out.append(integrate_step(Tube(tdomain, h, interval), v, dtheta))
	
	return union(out), out


if __name__ == "__main__":
	n = 50
	h = 1/20
	phi = Interval(-np.pi/2, np.pi/2)

	I, out = integrate(phi, n, h)

	vibes.beginDrawing()
	vibes.newFigure("MagMap")
	vibes.setFigureProperties({"x": 200, "y": 200, "width": 600, "height": 600})
	vibes.axisLimits(-np.pi/4, 3*np.pi/4, -np.pi/2, np.pi/2)
	for interval in sub_interval(phi, n):
		vibes.drawBox(0, h, interval[0], interval[1], "black[darkCyan]")
	for interval in out:
		vibes.drawBox(h, 2*h, interval.codomain()[0], interval.codomain()[1], "black[darkCyan]")

	vibes.drawBox(3*h, 4*h, -np.pi/2, np.pi/2, "black[darkCyan]")
	vibes.drawBox(4*h, 5*h, I.codomain().lb(), I.codomain().ub(), "black[darkCyan]")