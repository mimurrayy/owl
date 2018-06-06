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