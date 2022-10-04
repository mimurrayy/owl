#!/usr/bin/python3

import numpy as np
from ..util import *
from .basic import *
from scipy.interpolate import interp1d
import os  

class newport(basic):
    def __init__(self, cw, order=1, grating_constant=1200, size=1024):

        self.cw = cw
        self.order = order
        self.grating_constant = grating_constant
        self.x = self.make_x_scale(cw,order,size)
        self.fine_x = self.make_x_scale(cw,order,64*1024)


    def make_x_scale(self, cw=None, order=1, size=1024):
        """Calculates the wavelength scale. Can use bigger sizes to
        create a better resolution for the plots of the results.
        Only tested for the 1200er grating in first order."""
        if cw == None:
            cw = self.cw
        if order == None:
            order = self.order

        if np.round(self.grating_constant) == 1200:
            n = 1197e3 # effective grating constant
        L = 260e-3 # focal length
        D = -23.66 * (np.pi/180) # (fixed) angle between the grating and the two mirrors
        w = 24.56e-3 # width cam chip
        shift = -0.8
        xc0 = cw*1e-9 + shift*1e-9

        alpha = np.linspace(-np.pi/2, np.pi/2, 1000)
        beta = alpha + D
        wl = 2/(order*n) * np.sin(D/2 + alpha) * np.cos(D/2)  

        beta_wl = interp1d(wl, beta)

        z = np.linspace(-w/2, w/2, size)
        dz = z[1] - z[0]

        phi = -np.arctan(z/L)
        dwl = 1/(order*n*L) * (np.cos(beta_wl(xc0) - phi)*np.cos(phi)**2) * dz

        x = np.cumsum(dwl) 
        x = x - 0.5*(x[int(size/2)] + x[int(size/2-1)])
        x = x + xc0
        return x*1e9 


    def instrument_function(self, x, xc, transition=None, params = None):
        if params:
            w,mu = params
        else:
           w,mu = self.get_instrumental_params(transition)
        y = psd_voigt(x,xc,w,mu)
        return y


    def get_instrumental_params(self, transition):
        "Database for instrumental functions for lines."
        if self.order == 1 and self.grating_constant == 1200:
            return (0.15872448825374236, 0.40376419128985563)
        
        if self.order == 2 and self.grating_constant == 1200:
            return (0.08889650144805517, 0.5249123107116069)


    def load_spec(self, path):
        if os.path.isdir(path):
            data = []
            for f in sorted(os.listdir(path)):
                y = np.genfromtxt(path + "/" + f, delimiter=",")[:-1][::-1]
                data.append(y)
            return np.array(data).T
                
        elif os.path.isfile(path):  
            return np.genfromtxt(path, delimiter=",")[:-1][::-1]
        
        else:
            print("Not a directory or file path.")
            


    # alias
    instr = instrument_function
    instrumental_function = instrument_function
    instrumental_profile = instrument_function
    load = load_spec