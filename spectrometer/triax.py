#!/usr/bin/python3

import numpy as np
from ..util import *
from .base_spectrometer import *

class triax(base_spectrometer):
    def __init__(self, cw, order=1, grating_constant=50):
        self.cw = cw
        self.order = order
        self.grating_constant = grating_constant
        self.x = None
        # self.fine_x = self.make_x_scale(cw,order,64*1024)


    def instrument_function(self, x, xc, transition=None, params = None):
        if params:
            w,mu = params
        elif transition:
           w,mu = self.get_instrumental_params(transition)
        y = psd_voigt(x,xc,w,mu)
        return y


    def get_instrumental_params(self, transition):
        "Database for instrumental functions for lines."
        if self.order == 1 and self.grating_constant == 50:
            return (2.48216407e+00, 3.63312977e-01)


    # alias
    instr = instrument_function
    instrumental_function = instrument_function
    instrumental_profile = instrument_function