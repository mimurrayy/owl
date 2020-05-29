#!/usr/bin/python

import os
from ..util import *
from . import transition
import numpy as np

class level():
    def __init__(self, particle, energy, debug=False):
        self.particle = particle
        self.element = particle.element
        self.charge = particle.charge
        self.E = energy # Energy in eV
        self.J = None # Angular momentum
        self.G = None # Lande g factor
        self.l = None # l quantum number
        self.g = None # statistical weight
        self.conf = None # configuration and term string
        try:
            self.nist_info(debug=debug)
        except SyntaxError:
            print("Warning: EOF error. Should be harmless.")
        except Exception as e:
            print("No NIST Databse available.")
            print(e)

    def l_name_to_num(self, name):
        chars = ["s","p","d","f","g","h","i","j"]
        num = chars.index(name)
        return num

    def nist_info(self, debug=False):
        "Load NIST Tables and set J,G and conf"
        folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "nist-db")
        levels_file = os.path.join(folder, (self.element.lower() + str(self.charge+1) + "-levels.txt"))
        spec_name = get_spectroscopic_name(self.element, self.charge)

        for line in open(levels_file, 'r').readlines():
            line = line.replace(" ", "")
            array = line.split("|")
            if len(array) > 1 and len(array[0])>2:
                current_conf = array[0] + "," + array[1]
            if str(self.E) in line:
                self.J = float(eval(array[2]))
                self.conf = current_conf
                try:
                    self.G = float(array[4])
                except:
                    if debug:
                        print("No Lande-g in NIST DB")

    def get_transitions(self, kind="all", debug=False):
        """Return the transition objects that belong to the energy level
        Type is a string and may be "from", "to" or "all", to select transitions
        from the level, to the level or both."""
        folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "nist-db")
        lines_file = os.path.join(folder, (self.element.lower() + "-lines.txt"))
        spec_name = get_spectroscopic_name(self.element, self.charge)

        wl_col = 1
        E_col = 6

        transitions = []

        for line in open(lines_file, 'r').readlines():
            if "Unc." in line and int_col == 3: # some NIST files contain WL uncertainties.
                E_col = E_col + 2
            if line.startswith(spec_name):
                line = line.replace(" ", "")
                array = line.split("|")
                this_col = array[E_col]
                this_col = this_col.replace("[","")
                this_col = this_col.replace("]","")
                this_col = this_col.replace("(","")
                this_col = this_col.replace(")","")
                upperE = this_col.split("-")[1]
                lowerE = this_col.split("-")[0]

                if (str(self.E) in upperE) and (kind=="all" or kind=="from"):
                    try:
                        observed_wl = float(array[wl_col])
                    except:
                        observed_wl = float(array[wl_col+1]) # use calc WL : (
                    wl = round(observed_wl, 3)
                    this_transition = transition(self.particle, wl, debug=debug)
                    transitions.append(this_transition)

                if (str(self.E) in lowerE) and (kind=="all" or kind=="to"):
                    try:
                        observed_wl = float(array[wl_col])
                    except:
                        observed_wl = float(array[wl_col+1]) # use calc WL : (
                    wl = round(observed_wl, 3)
                    this_transition = transition(self.particle, wl, debug=debug)
                    transitions.append(this_transition)

        return transitions


    def get_lifetime(self, debug=False):
        Aiks = []
        for transition in self.get_transitions(kind="from", debug=debug):
            if transition.Aik:
                Aiks.append(transition.Aik)
        sum_Aik = np.sum(Aiks)
        return 1/sum_Aik
