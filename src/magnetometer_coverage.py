from pyibex import *
from pyibex.thickset import *
from pyibex.geometry import *
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

My_CMap = {
  str(IN) : "#c0392b[#c0392b]",
  str(OUT): "#27ae60[#27ae60]",
  str(MAYBE): "#e67e22[#e67e22]",
  str(MAYBE_IN): "#d35400[#e67e22]",
  str(MAYBE_OUT) :"#d35400[#e67e22]",
  str(UNK) :  "#f1c40f[#f1c40f]",
  str(EMPTY) : "[g]"
}

if __name__ == "__main__":

    # The map
    m_map = IntervalVector(2, Interval(-5, 5))

    # All the initial feasible values
    y = Interval(-1, 0)

    f_dist = Function("x1", "x2", "p1", "p2", "(x1-p1)^2+(x2-p2)^2")
    f_angle = Function("x1", "x2", "p1", "p2", "atan((x2-p2)/(x1-p1))")

    L = Interval(3)
    L.inflate(0.1)
    phi = Interval(-np.pi/4, np.pi/4) + Interval(np.pi)
    robot_box = IntervalVector([[-1, 1], [-1, 1]])

    # Polar separator

    # S=SepPolarXY(L, phi)
    # S2 = ThickTest_from_ThickSep(S_circle)
    
    # Paving the thickset
    S = SepPolarXY(L, phi)
    f = Function("x1", "x2", "p1", "p2", "(x1-p1)^2+(x2-p2)^2")
    g = Function("x1", "x2", "p1", "p2", "atan((x2-p2)/(x1-p1))")
    S2 = SepTransform(S,f_dist,f_dist)
    #P = ThickPaving(m_map, ThickTest_from_ThickSep(S_polar), 0.05)
    
    vibes.beginDrawing()
    vibes.setFigureSize(m_map.max_diam(), m_map.max_diam())
    vibes.newFigure("ThickSet")
    vibes.setFigureProperties({'x':0,'y':0, 'width':500, 'height':500})
    pySIVIA(m_map, S, 0.1)
    #P.visit(ToVibes(figureName="ThickSet", color_map=My_CMap))
    vibes.endDrawing()
