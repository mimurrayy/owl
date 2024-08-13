#!/usr/bin/python

import numpy as np
import re
from ..util import parse_spectroscopic_name, get_spectroscopic_name
from .level import level

class transition():
    def __init__(self, emitter_name, wavelength, wl_type="Observed", debug=False):
        """Return the closest transition for the specified emitter to the
        specified wavelength.
        emitter_name in spectroscopic notation, e.g. O I or Ar II
        wl_type can be "Observed" or "Ritz" or "either"
        
        Returns a transition object with entries: upperE, lowerE, upperl, 
        lowerl, upperg, lowerg, Aik, wl and the upper 
        and lower levels of the transition: upper_level/lower_level."""
        
        self.name, self.charge = parse_spectroscopic_name(emitter_name)
        self.spec_name = get_spectroscopic_name(self.name, self.charge)
        from mendeleev import element
        self.emitter = self.particle = element(self.name)
        self.emitter.charge = self.charge
        self.emitter.m = self.emitter.mass
        self.emitter.Ei = self.emitter.ionenergies[1]
        self.charge = self.emitter.charge
        # clean up the wl column from NIST levels (e.g. remove '[') below
        non_decimal = re.compile(r'[^\d.]+')
                
        # Load data tables from NIST (+- 1 nm around requested wl)
        from astroquery.nist import Nist
        import astropy.units as u
        nist_lines = Nist.query((wavelength-1)*u.nm, (wavelength+1)*u.nm, 
                                linename=self.spec_name, wavelength_type='vac+air')
        
        # select closest line
        if wl_type == "Observed":
            line_idx = np.argmin((nist_lines['Observed']-wavelength)**2)
            self.wl = nist_lines['Observed'][line_idx]

        if wl_type == "Ritz":   
            wls = []
            for this_wl in nist_lines['Ritz']:
                try:
                    wls.append(float(non_decimal.sub('', str(this_wl))))
                except:
                    wls.append(0.0)
            line_idx = np.argmin((wls-wavelength)**2)
            self.wl = wls[line_idx]

        if wl_type == "either":
            wls = []
            for this_wl in nist_lines['Ritz']:
                try:
                    wls.append(float(non_decimal.sub('', str(this_wl))))
                except:
                    wls.append(0.0)
            this_wl = np.max((np.array(nist_lines['Observed'], dtype=float), wls), axis=0)
            line_idx = np.argmin((this_wl-wavelength)**2)
            self.wl = this_wl[line_idx]


        row = nist_lines[line_idx]
        
        self.upperE = row['Ei           Ek'].split('-')[1]
        self.upperE = float(non_decimal.sub('', str(self.upperE)))
        self.lowerE = row['Ei           Ek'].split('-')[0]
        self.lowerE = float(non_decimal.sub('', str(self.lowerE)))
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

    