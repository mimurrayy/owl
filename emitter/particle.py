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

        mass = {'H'  : 1.00794,
                'Ar' :  39.948,
                'Ti' :  47.867,
                }


        E_ion = {'H'  : 13.59844,
                 'Ar' :  15.7596,
                 'Ti' :   6.8281,}


        self.m = mass[self.element]
        self.Ei = E_ion[self.element]
