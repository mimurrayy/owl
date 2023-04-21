#!/usr/bin/python

import os
import re
from ..util import parse_spectroscopic_name
import numpy as np
import mendeleev 
from astroquery.nist import Nist
from astroquery.nist_levels import NistLevels
import astropy.units as u

class level():
    def __init__(self, emitter_name, energy, debug=False):
        self.name, self.charge = parse_spectroscopic_name(emitter_name)
        self.spec_name = emitter_name
        self.emitter = self.particle = mendeleev.element(self.name)
        self.emitter.charge = self.charge
        self.emitter.m = self.emitter.mass
        self.emitter.Ei = self.emitter.ionenergies[1]
        self.charge = self.emitter.charge
        
        nist_levels = NistLevels.query(linename=emitter_name)

        # clean up the energy column from NIST levels (e.g. remove '[')
        non_decimal = re.compile(r'[^\d.]+')
        level_energies = []
        for entry in nist_levels['Level (eV)']:
            level_energies.append(non_decimal.sub('', str(entry)))
        level_energies = np.array(level_energies, dtype=float)
        
        # select upper and lower energy levels
        level_idx = np.argmin((level_energies-energy)**2)
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
            self.G = nist_levels[level_idx]['Landé'] #Lande g factor
        else:
            self.G = None
            print("Warning: No Landé-g factors in the NIST asd database.")
        


    def l_name_to_num(self, name):
        chars = ["s","p","d","f","g","h","i","j"]
        num = chars.index(name)
        return num


    # def get_transitions(self, kind="all", debug=False):
    #     """Return the transition objects that belong to the energy level
    #     Type is a string and may be "from", "to" or "all", to select transitions
    #     from the level, to the level or both."""
    #     from .transition import transition # needed to avoid circular import

    #     folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "nist-db")
    #     lines_file = os.path.join(folder, (self.element.lower() + "-lines.txt"))
    #     spec_name = get_spectroscopic_name(self.element, self.charge)

    #     wl_col = 1
    #     E_col = 6

    #     transitions = []

    #     for line in open(lines_file, 'r').readlines():
    #         if "Unc." in line and int_col == 3: # some NIST files contain WL uncertainties.
    #             E_col = E_col + 2
    #         if line.startswith(spec_name):
    #             line = line.replace(" ", "")
    #             array = line.split("|")
    #             this_col = array[E_col]
    #             this_col = this_col.replace("[","")
    #             this_col = this_col.replace("]","")
    #             this_col = this_col.replace("(","")
    #             this_col = this_col.replace(")","")
    #             upperE = this_col.split("-")[1]
    #             lowerE = this_col.split("-")[0]

    #             if (str(self.E) in upperE) and (kind=="all" or kind=="from"):
    #                 try:
    #                     observed_wl = float(array[wl_col])
    #                 except:
    #                     observed_wl = float(array[wl_col+1]) # use calc WL : (
    #                 wl = round(observed_wl, 3)
    #                 this_transition = transition(self.particle, wl, debug=debug)
    #                 transitions.append(this_transition)

    #             if (str(self.E) in lowerE) and (kind=="all" or kind=="to"):
    #                 try:
    #                     observed_wl = float(array[wl_col])
    #                 except:
    #                     observed_wl = float(array[wl_col+1]) # use calc WL : (
    #                 wl = round(observed_wl, 3)
    #                 this_transition = transition(self.particle, wl, debug=debug)
    #                 transitions.append(this_transition)

    #     return transitions


    # def get_lifetime(self, debug=False):
    #     Aiks = []
    #     for transition in self.get_transitions(kind="from", debug=debug):
    #         if transition.Aik:
    #             Aiks.append(transition.Aik)
    #     sum_Aik = np.sum(Aiks)
    #     return 1/sum_Aik
