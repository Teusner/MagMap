import numpy as np
from pyibex import *
from tubex_lib import *

if __name__ =="__main__":
    h = 0.1
    tdomain = Interval(0, 5)
    x = TubeVector(tdomain, h, 1)
    x.set(IntervalVector(1, Interval(-np.pi/4, np.pi/4)), 0.)

    f = Function("x", "-sin(x)")

    ctc_lohner = CtcLohner(f)
    ctc_lohner.contract(x)

    beginDrawing()
    fig2 = VIBesFigTube("Lohner integration")
    fig2.set_properties(100, 550, 800, 400)
    fig2.add_tube(x[0], "X")
    fig2.show()
    endDrawing()
