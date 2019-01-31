#!/usr/bin/python3

import numpy as np
from ..util import *

class basic():
    def __init__(self, cw=None, w=None):
        if cw:
            self.x = self.make_x_scale(cw)
        else:
            self.x = None

        if w:
            self.w = w
        else:
            self.w = None


    def instrument_function(self, x, xc, transition=None, params = None, w = None):
        if not w and not self.w:
            w = 0.001
        elif self.w:
            w = self.w
        y = gauss_function(x,xc,w)
        return y

    def make_x_scale(self, cw):
        x = np.linspace(cw-50,cw+50,5000)
        return x


    # alias
    instr = instrument_function
    instrumental_function = instrument_function



