#!/usr/bin/python

import os
import re
from ..util import parse_spectroscopic_name
import numpy as np

class level():
    def __init__(self, emitter_name, energy, debug=False):
        self.name, self.charge = parse_spectroscopic_name(emitter_name)
        self.spec_name = emitter_name
        from mendeleev import element
        self.emitter = self.particle = element(self.name)
        self.emitter.charge = self.charge
        self.emitter.m = self.emitter.mass
        self.emitter.Ei = self.emitter.ionenergies[1]
        self.charge = self.emitter.charge
        
        from ..nist_levels import NistLevels
        nist_levels = NistLevels.query(linename=emitter_name)

        # clean up the energy column from NIST levels (e.g. remove '[')
        non_decimal = re.compile(r'[^\d.]+')
        level_energies = []

        for entry in nist_levels['Level (eV)']:
            this_entry = non_decimal.sub('', str(entry))
            try: 
                this_entry = float(this_entry)
            except:
                this_entry = np.nan
            level_energies.append(this_entry)
        level_energies = np.array(level_energies, dtype=float)

        # select upper and lower energy levels
        level_idx = np.argmin((np.nan_to_num(level_energies)-energy)**2)
        row = nist_levels[level_idx]
        
        self.E = level_energies[level_idx] # Energy in eV
        try:
            self.J = float(eval(row['J'])) # Angular momentum, can be 5/2 or so...
        except:
            self.J = np.nan
        
        try:
            conf_end = row['Configuration'].split('.')
            if "(" in conf_end[-1]:
                conf_end.pop()
            self.l = int(self.l_name_to_num(conf_end[-1][1]))
        except:
            self.l = np.nan
        
        self.g = float(nist_levels[level_idx]['g']) # statistical weight
        self.conf = row['Configuration'] # configuration and term string
        
        if 'Landé' in nist_levels.keys():
            self.G = float(nist_levels[level_idx]['Landé']) #Lande g factor
        else:
            self.G = None
        


    def l_name_to_num(self, name):
        chars = ["s","p","d","f","g","h","i","j"]
        num = chars.index(name)
        return num
    
    
    def get_transitions(self, kind="all", debug=False):
        """Return the transition objects that belong to the energy level
        Type is a string and may be "from", "to" or "all", to select transitions
        from the level, to the level or both."""
        from .transition import transition # needed to avoid circular import

        # clean up the wl column from NIST levels (e.g. remove '[') below
        non_decimal = re.compile(r'[^\d.]+')

        # Load data tables from NIST (+- 1 nm around requested wl)
        from astroquery.nist import Nist # import here for startup perfromance
        import astropy.units as u
        nist_lines = Nist.query(1*u.nm, 99999*u.nm, 
                                linename=self.spec_name, wavelength_type='vac+air')

        transitions_to = []
        transitions_from = []
        transitions_all = []

        for entry in nist_lines:
            if not np.ma.is_masked(entry['Ei           Ek']) \
                                    and not (np.ma.is_masked(entry['Observed']) and np.ma.is_masked(entry['Ritz']) ):
                
                if '-' in entry['Ei           Ek']:
                    upperE = float(entry['Ei           Ek'].split('-')[1]) # Energy in eV
                    lowerE = float(entry['Ei           Ek'].split('-')[0])
                    if upperE == self.E:
                        try:
                            wl = float(entry['Observed'])
                        except:
                            try:
                                wl = float(non_decimal.sub('', str(entry['Ritz'])))
                            except:
                                break
                        if not np.isnan(wl):
                            this_transition = transition(self.spec_name, wl, wl_type="either")
                            transitions_from.append(this_transition)                    
                            transitions_all.append(this_transition)
                    
                    if lowerE == self.E:
                        try:
                            wl = float(entry['Observed'])
                            if np.isnan(wl):
                                wl = float(non_decimal.sub('', str(entry['Ritz'])))
                        except:
                            try:
                                wl = float(non_decimal.sub('', str(entry['Ritz'])))
                            except:
                                break
                        if not np.isnan(wl):
                            this_transition = transition(self.spec_name, wl, wl_type="either")
                            transitions_to.append(this_transition)                    
                            transitions_all.append(this_transition)

        if kind=="all":
            return transitions_all
        if kind=="from":
            return transitions_from        
        if kind=="to":
            return transitions_to
        
        return None


    def get_lifetime(self, debug=False):
        Aiks = []
        for transition in self.get_transitions(kind="from", debug=debug):
            if transition.Aik:
                Aiks.append(transition.Aik)
        sum_Aik = np.sum(Aiks)
        return 1/sum_Aik
