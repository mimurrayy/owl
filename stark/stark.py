#!/usr/bin/python

import numpy as np
from ..util import *
from .gigosos_loader import *
from .gigosos_he_loader import *
from .griem import *

class stark():
    def __init__(self, transition, ion_pert = None):
        self.transition = transition
        self.fast = False
        self.pert = ion_pert


    def get_profile(self,x, ne, Te=None, pert=None):
        middle_wl = x[int(len(x)/2)]
        ele = self.transition.element

        ################ Hydrogen ########
        if ele == "H":
            if round(self.transition.wl,0) == 656:
                if self.fast or not Te:
                    # Gigosos et al., 2003: Computer simulated Balmer-alpha,
                    # -beta and -gamma Stark line profiles for non-equilibrium
                    # plasmas diagnostics
                    w = ((ne/1e23)**(0.67965)) * 1.098
                    return lorentz_function(x,middle_wl,w)
                else:
                    this_loader = gigosos_loader(self.transition)
                    gigosos_x,y = this_loader.load(ne, Te, pert)
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
                    gigosos_x,y = this_loader.load(ne, Te, pert)
                    gigosos_x = gigosos_x + middle_wl # to nm
                    y = interpol(gigosos_x,y,x)
                    return y

        ################ Helium ##########
        if ele == "He":
            if round(self.transition.wl,0) == 447:
                this_loader = gigosos_he_loader(self.transition)
                gigosos_x,y = this_loader.load(ne, Te, pert)
                gigosos_x = gigosos_x + middle_wl # to nm
                y = interpol(gigosos_x,y,x)
                return y

            if round(self.transition.wl, 0) == 492:
                this_loader = gigosos_he_loader(self.transition)
                gigosos_x,y = this_loader.load(ne, Te, pert)
                gigosos_x = gigosos_x + middle_wl # to nm
                y = interpol(gigosos_x,y,x)
                return y

        ################ Others: Griem ##########
        if ele == "O" or "Ar":
            this_griem = griem(self.transition)
            w,d = this_griem.get_width_shift(ne, Te)
            # shift is ignored for now
            return lorentz_function(x,middle_wl,w)


    def get_width(self, ne, Te=None, pert=None):
        ele = self.transition.element
        if ele == "O" or ele == "Ar":
            this_griem = griem(self.transition)
            w,d = this_griem.get_width_shift(ne, Te)
            return w

    def get_shift(self, ne, Te=None, pert=None):
        ele = self.transition.element
        if ele == "O" or ele == "Ar" :
            this_griem = griem(self.transition)
            w,d = this_griem.get_width_shift(ne, Te)
            return d


