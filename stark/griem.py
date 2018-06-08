#!/usr/bin/python3

import numpy as np
from scipy import constants as const
from ..util import *
import os

class griem():
    def __init__(self, transition):
        self.transition = transition

    def get_width_shift(self, ne, Te):
        if self.transition.element == "O":
            A,B,we,de = self.params_O777(Te)

        w = self.width(ne,Te,A,we)
        d = self.shift(ne,Te,A,we,de)

        return w,d


    def width(self,ne,Te,A,we):
        """ Te in K, ne in m^-3
        we, de in Angstrom @ 1e22/m^3 density
        Returns w in nm
        """

        # we need Te in K
        T = Te * const.eV / const.k
        # we need ne in cm^-3
        n = ne * 1e-6

        w = 2 * we * n * 1e-16 * (1 +
            1.75 * 1e-4 * n**(1/4) * A * (1 - 0.068*n**(1/6) * T**(-1/2) ) )

        # w is Angstrom
        w = w/10
        return w


    def shift(self,ne,Te,A,we,de):
        """ Te in K, ne in m^-3
        we, de in Angstrom @ 1e22/m^3 density """

        # we need Te in K
        T = Te * const.eV / const.k
        # we need ne in cm^-3
        n = ne * 1e-6

        d = n*1e-16 * (de + 2e-4 * n**(1/4) * A * we * (1-0.068 * n**(1/6) * T**(-1/2)))

        # d is Angstrom
        d = d/10
        return d


    def params_O777(self,Te):
        """ Parameter from:
         H.R. Griem: Spectral Line Broadening by Plasmas
         returns A,B and we,de in nm
         All values at ne = 10^16 /cm^3
         """

        # we need Te in K
        T = Te * const.eV / const.k

        B = 0.00023718 * T # recommended by Griem
        A = 0.27661461 * T**(-0.33635335) # Just fits okay.
        we = 4.36767799e-04 * T**(0.465538767) # Fitted; in Angstrom

        a0,a1,a2,a3,a4 = [1.39308003e-02,1.88771484e-07,-2.17304732e-11,
        6.29115391e-16,-6.33972112e-21] # very much overdefined...
        de = a0 + a1*T + a2*T**2 + a3*T**3 + a4*T**4

        return A,B,we,de