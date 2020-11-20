from pyibex import *
from pyibex.thickset import *
import numpy as np


class ThickSep:
    def __init__(self, Ssub, Ssup): 
        self.Ssub, self.Ssup = Ssub, Ssup

    def separate(self, X):
       xin_sub, xout_sub = X.copy(), X.copy()
       xin_sup, xout_sup = X.copy(), X.copy()
       self.Ssub.separate(xin_sub, xout_sub)
       self.Ssup.separate(xin_sup, xout_sup)
       return (xin_sup, xin_sub | xout_sup, xout_sub)

    def __and__(self, thickSep_y): 
        return ThickSep(self.Ssub & thickSep_y.Ssub, self.Ssup & thickSep_y.Ssup)

    def __invert__(self): 
        return ThickSep(~self.Ssup, ~self.Ssub)

    def __or__(self, thickSep_y): 
        return ThickSep(self.Ssub | thickSep_y.Ssub, self.Ssup | thickSep_y.Ssup)

def ThickSep_from_function(f, p, y):
    S = SepFwdBwd(f, y)
    Ssub = ~SepProj(~S, p)
    Ssup = SepProj(S, p)
    return ThickSep(Ssup, Ssub)

def ThickTest_from_ThickSep(S):
    def test(X):
       Xin, Xu, Xout = S.separate(X)
       if Xin.is_empty(): return IN
       elif Xout.is_empty(): return OUT
       elif Xu.is_empty(): return MAYBE
       return UNK
    return test

if __name__ == "__main__":

    # The map
    m_map = IntervalVector(2, Interval(-10, 10))

    # All the initial feasible values
    y = Interval(-oo, 0)

    # Contracting y
    f_dist = Function("x1", "x2", "p1", "p2", "p3", "(x1-p1)^2+(x2-p2)^2-p3^2")

    traj_x = np.array([0]) #np.arange(-15, 15, 0.2)
    traj_y = 10 * np.exp(-(traj_x+5)/10) * np.sin(traj_x/4)

    # Initialisation because I am not able to setup S_dist an empty separator
    p = IntervalVector(3, Interval())
    p[0] = Interval(traj_x[0])
    p[1] = Interval(traj_y[0])
    p[2] = Interval(5)
    p[0].inflate(2)
    p[1].inflate(2)
    S_dist = ThickSep_from_function(f_dist, p, y)

    for k in range(1, len(traj_x)):
        # Parameters of the function circle function [[x0], [y0], [r0]]
        p = IntervalVector(3, Interval())
        p[0] = Interval(traj_x[k])
        p[0].inflate(2)
        p[1] = Interval(traj_y[k])
        p[1].inflate(2)
        p[2] = Interval(5)
        S_dist = S_dist | ThickSep_from_function(f_dist, p, y)

    # Paving the thickset
    P = ThickPaving(m_map, ThickTest_from_ThickSep(S_dist), 0.5, display=True)

    # vibes.beginDrawing()
    # v = ToVibes(P)
    # vibes.endDrawing()