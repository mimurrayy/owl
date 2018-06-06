#!/usr/bin/python3

import numpy as np
from ..util import *

class base_spectrometer():
    def __init__(self):
        self.x = None

    def instrument_function(self, x, xc, transition=None, params = None):
        y = np.zeros(x)
        y[int(len(x)/2)] = 1
        return y

    # alias
    instr = instrument_function
    instrumental_function = instrument_function



