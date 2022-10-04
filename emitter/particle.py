#!/usr/bin/python

import os
from scipy import constants as const
from ..util import *

class particle():
    def __init__(self, element, charge=0):
        self.T = None
        self.element,self.charge = parse_spectroscopic_name(element)

        self.name = self.element

        mass = {'H'  :   1.00794,
                'He' :   4.00260,
                'C'  :  12.011,
                'N'  :  14.0067,
                'O'  :  15.999,
                'H2O':  18.015,
                'Al' :  26.981539,
                'Si' :  28.0855,
                'Cl' :  35.45,
                'Ar' :  39.948,
                'Ti' :  47.867,
                'V'  :  50.9415,
                'Cr' :  51.9961,
                'Fe' :  55.845,
                'Cu' :  63.546,
                'Nb' :  92.90638,
                'Mo' :  95.94,
                'Ta' : 180.94788,
                'W'  : 183.84,
                }

        E_ion = {'H'  : 13.59844,
                 'He' : 24.58741,
                 'C'  : 11.26030,
                 'N'  : 14.53414,
                 'O'  : 13.61806,
                 'Al' :  5.98577,
                 'Si' :  8.15169,
                 'Cl' : 12.96764, 
                 'Ar' : 15.7596,
                 'Ti' :  6.8281,
                 'V'  :  6.7462,
                 'Cr' :  6.7665,
                 'Fe' :  7.9024,
                 'Cu' :  7.72638,
                 'Nb' :  6.75885,
                 'Mo' :  7.09243,
                 'Ta' :  7.5496,
                 'W'  :  7.8640,
                 }

        self.m = mass[self.element]
        try:
            self.Ei = E_ion[self.element]
        except:
            print("Ionization energy unknown.")

    def spectroscopic_name(self):
        return get_spectroscopic_name(self.element, self.charge)

