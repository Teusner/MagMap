import numpy as np
from vibes import vibes
from pyibex import *
from tubex_lib import *


class Mapping:
    def __init__(self, t0, tf, h):
        self.h = h

        # Trajectory Tubes
        self.tdomain = Interval(t0, tf)
        self.x = TubeVector(self.tdomain, self.h, IntervalVector(2))
        self.v = TubeVector(self.tdomain, self.h, IntervalVector(2))

        # Contractor Network
        self.cn = ContractorNetwork()
        ctc_deriv = CtcDeriv()
        self.cn.add(ctc_deriv, [self.x, self.v])

    def add_position(self, t, position):
        n = len(t)
        for k in range(n):
            self.cn.add(ctc_eval, [t[i], position[i], self.x, self.v])
    
    def add_velocity(self, t, velocity):
        n = len(t)
        for k in range(n):
            self.v.set(velocity[i], Interval(t[i], t[i] + self.h))


if __name__ == "__main__":
    m = Mapping(0, 10, 1/20)