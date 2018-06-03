#!/usr/bin/python

import os
from scipy import constants as const
from ..util import *

class particle():
    def __init__(self, element, charge=0):
        if len(element) > 2:
            self.element,self.charge = parse_spectroscopic_name(element)
        else:
            self.element = element.title()
            self.charge = charge  
        mass = {'Ti' :  47.867,
                'Ar' :  39.948 }
        self.m = mass[self.element]
        E_ion = {'Ti' :   6.8281,
                 'Ar' :  15.7596 }
        self.Ei = E_ion[self.element]

