#!/usr/bin/python3

import numpy as np
from ..util import *
from .basic import *

class pgs(basic):
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

    def get_instrumental_params(self, transition = None, xc = None):
        "Database for instrumental functions for lines."
        if transition:
            wl = transition.wl
        if xc:
            wl = xc
        else:
            wl = 0

        if self.order == 3 and round(wl,0) == 453.0:
            return (0.00653, 0.00173, 0.7795, 1.0, 1.0295, 0.005, 0.0115)
        if self.order == 3 and round(wl,0) == 399:
            return (0.00653, 0.00173, 0.7795, 1.0, 1.0295, 0.005, 0.0115)
        if self.order == 3 and round(wl,0) == 395:
            return (0.00810, 0.00349, 0.7094, 0.9947, 0.4724, 0.0074, 0.0082)
        if self.order == 3 and round(wl,0) == 473:
            return (0.00476, 0.00324, 1.1660, 0.6125, 0.9649, 0.0063, 0.0098)
        if round(wl,0) == 633 and self.order == 2:
            return (0.00506, 0.00598, 0.7877, 0.5498, 0.6312, -5e-06, 0.0442)
        elif self.order == 3:
            return (0.00653, 0.00173, 0.7795, 1.0, 1.0295, 0.005, 0.0115)
        elif self.order == 1:
            return (0.01478, 0.01478, 0.5845, 0.5845, 0, 0, 0)
        else:
            return (0.00506, 0.00598, 0.7877, 0.5498, 0.6312, -5e-06, 0.0442)

    def instrument_function(self, x, xc, transition=None, params = None):
        """ Derived with a HeNe-Laser for 1,2 Order and a HollowCathode Ti Line for 3 Order
        and a Ne II Line for the 4th Order.
        Fit done with asymmetric linear pseudo Voigt,
        for no physical reasons. Just seemed to fit very well."""
       # H = 1
        if params:
            w1,w2,mu1,mu2,Hg,dx,w3 = params
        elif transition:
            w1,w2,mu1,mu2,Hg,dx,w3 = self.get_instrumental_params(
                transition=transition)
        elif xc:
            w1,w2,mu1,mu2,Hg,dx,w3 = self.get_instrumental_params(xc=xc)

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

    def load_spec(self, f):
        num_lines = sum(1 for line in open(f))
        columns = np.genfromtxt(f,skip_footer=num_lines-1027).T#[2][::-1]
        x = self.x
        data = []
        data.append(x)
        for y in columns[1:]:
            if self.order == 1:
                data.append(y)
            else:
                data.append(y[::-1])
        return np.array(data)

    # alias
    instr = instrument_function
    instrumental_function = instrument_function
    instrumental_profile = instrument_function


