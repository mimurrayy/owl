#!/usr/bin/python

import os
from scipy import constants as const
from ..util import *

class particle():
    def __init__(self, element, charge=0):
        self.T = None
        if len(element) > 2:
            self.element,self.charge = parse_spectroscopic_name(element)
        else:
            self.element = element.title()
            self.charge = charge  

        mass = {'H'  :   1.00794,
                'C'  :  12.011,
                'N'  :  14.0067,
                'O'  :  15.999,
                'Al' :  26.981539,
                'Ar' :  39.948,
                'Ti' :  47.867,
                'V'  :  50.9415,
                'Cr' :  51.9961,
                'Cu' :  63.546,
                'Nb' :  92.90638,
                'Mo' :  95.94,
                'Ta' : 180.94788,
                'W'  : 183.84,
                }


        E_ion = {'H'  : 13.59844,
                 'C'  : 11.26030,
                 'N'  : 14.53414,
                 'O'  : 13.61806,
                 'Al' :  5.98577,
                 'Ar' : 15.7596,
                 'Ti' :  6.8281,
                 'V'  :  6.7462,
                 'Cr' :  6.7665,
                 'Cu' :  7.72638,
                 'Nb' :  6.75885,
                 'Mo' :  7.09243,
                 'Ta' :  7.5496,
                 'W'  :  7.8640,
                 }


        self.m = mass[self.element]
        self.Ei = E_ion[self.element]

    def spectroscopic_name(self):
        return get_spectroscopic_name(self.element, self.charge)

