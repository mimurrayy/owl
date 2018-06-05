#!/usr/bin/python

import numpy as np
from ..util import *
from .gigosos_loader import *

class stark():
    def __init__(self, transition, ion_perp = None):
        self.transition = transition
        self.fast = False
        self.perp = ion_perp

    def get_profile(self,x, ne, Te=None, perp=None):
        middle_wl = x[int(len(x)/2)]
        if self.transition.element == "H":
            if round(self.transition.wl,0) == 656:
                if self.fast or not Te:
                    # Gigosos et al., 2003: Computer simulated Balmer-alpha,
                    # -beta and -gamma Stark line profiles for non-equilibrium
                    # plasmas diagnostics
                    w = ((ne/1e23)**(0.67965)) * 1.098
                    return lorentz_function(x,middle_wl,w)

                else:
                    this_loader = gigosos_loader(self.transition)
                    gigosos_x,y = this_loader.load(ne, Te, perp)
                    gigosos_x = gigosos_x + middle_wl # to nm
                    y = interpol(gigosos_x,y,x)
                    return y


            if round(self.transition.wl, 0) == 486:
                if self.fast or not Te:
                    # Gigosos et al., 2003: Computer simulated Balmer-alpha,
                    # -beta and -gamma Stark line profiles for non-equilibrium
                    # plasmas diagnostics
                    w = ((ne/1e23)**(0.68116)) * 4.8
                    return lorentz_function(x,xc,w)
                else:
                    this_loader = gigosos_loader(self.transition)
                    gigosos_x,y = this_loader.load(ne, Te, perp)
                    y = interpol(gigosos_x,y,x)
                    return y


    