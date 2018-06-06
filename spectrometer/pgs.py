#!/usr/bin/python3

import numpy as np
from ..util import *
from .base_spectrometer import *

class pgs(base_spectrometer):
    def __init__(self, cw, order=1, grating_constant=1302.26, cam_shift=49):
        self.cw = cw
        self.order = order
        self.grating_constant = grating_constant
        self.cam_shift = cam_shift
        self.x = self.make_x_scale(cw,order,1024)
        self.fine_x = self.make_x_scale(cw,order,64*1024)     

    def make_x_scale(self, cw=None, order=None, size=1024):
        """Calculates the wavelength scale. Can use bigger sizes to
        create a better resolution for the plots of the results.
        Does not influence calculations though. Use sizes of 2^n"""
        if cw == None:
            cw = self.cw
        if order == None:
            order = self.order
        grating_constant = self.grating_constant        
        cam_shift = self.cam_shift

        pixel_array = np.arange(size)
        angle = np.arcsin(grating_constant/2*(cw*10**-6)*order)
        delta_wl = ((np.cos(angle) *
            1/(grating_constant)*2*(10**6))*(0.026*1024/size*1/order*
            (pixel_array-(size/2-1)+cam_shift*size/1024)/4150)) # formula to get wavelenght scale
        if order == 1:
            return (cw + delta_wl)
        else:
            delta_wl = ((np.cos(angle) *
            1/(grating_constant)*2*(10**6))*(0.026*1024/size*1/order*
            (pixel_array - (size/2)- (cam_shift)*size/1024)/4150)) # formula to get wavelenght scale
            return (cw + delta_wl)

    def get_instrumental_params(self, transition):
        "Database for instrumental functions for lines."
        if transition.element == "Ti" and round(transition.wl,0) == 453.0:
            return (0.00653, 0.00173, 0.7795, 1.0, 1.0295, 0.005, 0.0115)
        if transition.element == "Ti" and transition.charge == 0 and round(transition.wl,0) == 399:
            return (0.00653, 0.00173, 0.7795, 1.0, 1.0295, 0.005, 0.0115)
        if transition.element == "Ar" and transition.charge == 1 and round(transition.wl,0) == 473:
            return (0.00416, 0.00322, 1.2025, 0.0, 2.5765, 0.0029, 0.0139)
        else:
            return (0.00416, 0.00322, 1.2025, 0.0, 2.5765, 0.0029, 0.0139)

    def instrument_function(self, x, xc, transition=None, params = None):
        """ Derived with a HeNe-Laser for 1,2 Order and a HollowCathode Ti Line for 3 Order
        and a Ne II Line for the 4th Order.
        Fit done with asymmetric linear pseudo Voigt,
        for no physical reasons. Just seemed to fit very well."""
       # H = 1
        if params:
            w1,w2,mu1,mu2,Hg,dx,w3 = params
        elif transition:
            w1,w2,mu1,mu2,Hg,dx,w3 = self.get_instrumental_params(transition)
        lx = x[x<=xc]
        left_part = (0.5 * w1 / (mu1/np.pi + (1-mu1) * np.sqrt(np.log(2)/np.pi))
            * ( mu1 * (2/np.pi) * (w1 / (4*(lx-xc)**2 + w1**2))
            + (1 - mu1) * (np.sqrt(4*np.log(2)) / (np.sqrt(np.pi) * w1))
            * np.exp(-(4*np.log(2)/w1**2)*(lx-xc)**2) ))
        rx = x[x>xc]
        right_part = (0.5 * w2 / (mu2/np.pi + (1-mu2) * np.sqrt(np.log(2)/np.pi))
            * ( mu2 * (2/np.pi) * (w2 / (4*(rx-xc)**2 + w2**2))
            + (1 - mu2) * (np.sqrt(4*np.log(2)) / (np.sqrt(np.pi) * w2))
            * np.exp(-(4*np.log(2)/w2**2)*(rx-xc)**2) ))
        y = np.append(left_part, right_part)
        if Hg:
            y = y + Hg*gauss_function(x,xc+dx,w3)*1e-3
        return y

    # alias
    instr = instrument_function
    instrumental_function = instrument_function



