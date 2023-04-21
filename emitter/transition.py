#!/usr/bin/python

import os
import numpy as np
from ..util import parse_spectroscopic_name, get_spectroscopic_name
from .level import level
import mendeleev 
from astroquery.nist import Nist
import astropy.units as u

class transition():
    def __init__(self, emitter_name, wavelength, debug=False):
        self.name, self.charge = parse_spectroscopic_name(emitter_name)
        self.spec_name = emitter_name
        self.emitter = self.particle = mendeleev.element(self.name)
        self.emitter.charge = self.charge
        self.emitter.m = self.emitter.mass
        self.emitter.Ei = self.emitter.ionenergies[1]
        self.charge = self.emitter.charge
        
        # Load data tables from NIST (+- 1 nm around requested wl)
        nist_lines = Nist.query((wavelength-1)*u.nm, (wavelength+1)*u.nm, 
                                linename=emitter_name, wavelength_type='vac+air')
        
        # select closest line
        line_idx = np.argmin((nist_lines['Observed']-wavelength)**2)
        row = nist_lines[line_idx]

        self.upperE = float(row['Ei           Ek'].split('-')[1]) # Energy in eV
        self.lowerE = float(row['Ei           Ek'].split('-')[0])
        # self.upperJ = row['Upper level'].split('|')[-1] # Angular momentum
        # self.lowerJ = row['Lower level'].split('|')[-1]
        
        upper_conf_end = row['Upper level'].split('|')[0].split('.')
        if "(" in upper_conf_end[-1]:
            upper_conf_end.pop()
        lower_conf_end = row['Lower level'].split('|')[0].split('.')
        if "(" in lower_conf_end[-1]:
            lower_conf_end.pop()
        try:
            self.upperl = self.l_name_to_num(upper_conf_end[-1][1])
            self.lowerl = self.l_name_to_num(lower_conf_end[-1][1]) 
        except:
            self.upperl = None
            self.lowerl = None
        
        self.Aik    = row['Aki']
        self.upperg = float(row['gi   gk'].split('-')[1])
        self.lowerg = float(row['gi   gk'].split('-')[0])
        self.wl = row['Observed']
        

        self.upper_level, self.lower_level = self.levels()
        self.upper, self.lower = self.upper_level, self.lower_level

    def l_name_to_num(self, name):
        chars = ["s","p","d","f","g","h","i","j"]
        num = chars.index(name)
        return num

    def levels(self):
        upper_level = level(self.spec_name, self.upperE)
        lower_level = level(self.spec_name, self.lowerE)
        self.upperJ = upper_level.J
        self.upperG = upper_level.G
        self.lowerJ = lower_level.J
        self.lowerG = lower_level.G
        return upper_level, lower_level

    