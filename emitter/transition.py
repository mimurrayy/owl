#!/usr/bin/python

import os
from ..util import *

class transition():
    def __init__(self, particle, transition_wavelength):
        self.particle = particle
        self.element = particle.element
        self.charge = particle.charge
        self.wl = round(transition_wavelength, 3)
        self.upperE = None # Energy in eV
        self.lowerE = None
        self.upperJ = None # Angular momentum
        self.lowerJ = None
        self.upperG = None # Lande g factor
        self.lowerG = None
        self.upperl = None # l quantum number
        self.lowerl = None

        try:
            self.nist_info()
        except:
            print("No NIST Databse available.")

    def l_name_to_num(self, name):
        chars = ["s","p","d","f","g","h","i","j"]
        num = chars.index(name)
        return num

    def nist_info(self):
        "Load NIST Tables and set E,J and G"
        folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "nist-db")
        lines_file = os.path.join(folder, (self.element.lower() + "-lines.txt"))
        levels_file = os.path.join(folder, (self.element.lower() + str(self.charge+1) + "-levels.txt"))
        spec_name = get_spectroscopic_name(self.element, self.charge)
        name_col = 0
        wl_col = 1
        int_col = 3
        aik_col = 4
        E_col = 6
        lowconf_col = 7
        upconf_col = 10
        for line in open(lines_file, 'r').readlines():
            if "Unc." in line and int_col == 3: # some NIST files contain WL uncertainties.
                int_col = int_col + 2
                aik_col = aik_col + 2
                E_col = E_col + 2
                lowconf_col = lowconf_col + 2
                upconf_col = upconf_col + 2
            if line.startswith(spec_name):
                line = line.replace(" ", "")
                array = line.split("|")
                if len(array[wl_col]) > 3:
                    observed_wl = float(array[wl_col])
                    if self.wl == round(observed_wl, 3):
                        this_col = array[E_col]
                        this_col = this_col.replace("[","")
                        this_col = this_col.replace("]","")
                        this_col = this_col.replace("(","")
                        this_col = this_col.replace(")","")
                        self.upperE = float(this_col.split("-")[1])
                        self.lowerE = float(this_col.split("-")[0])
                        upper_conf = array[upconf_col]
                        lower_conf = array[lowconf_col]
                        self.upperl = self.l_name_to_num(upper_conf.split('.')[-1][1])
                        self.lowerl = self.l_name_to_num(lower_conf.split('.')[-1][1])

        for line in open(levels_file, 'r').readlines():
            if str(self.upperE) in line:
                line = line.replace(" ", "")
                array = line.split("|")
                self.upperJ = float(eval(array[2]))
                try:
                    self.upperG = float(array[4])
                except:
                    print("No Lande-g in NIST DB")
            if str(self.lowerE) in line:
                line = line.replace(" ", "")
                array = line.split("|")
                self.lowerJ = float(eval(array[2]))
                try:
                    self.lowerG = float(array[4])
                except:
                    print("No Lande-g in NIST DB")

def load_nist_lines(self, particle):
        emission_lines = []
        folder = os.path.join(os.path.dirname(os.path.realpath(__file__)),
            "emitter/nist-db")

        unwanted_chars = ["q","[","]","†","?"]

        lines_file = os.path.join(folder, (particle.element.lower() + "-lines.txt"))
        spec_name = particle.spectroscopic_name()
        name_col = 0
        wl_col = 1
        int_col = 3
        aik_col = 4
        E_col = 6
        for line in open(lines_file, 'r').readlines():
            if "Unc." in line: # some NIST files contain WL uncertainties.
                int_col = int_col + 2
                aik_col = aik_col + 2
                E_col = E_col + 2

            if line.startswith(spec_name):
                emission_line = {}
                emission_line['particle'] = spec_name
                observed_wl = rel_int = Aik = Ei = Ek = None
                line = line.replace(" ", "")
                array = line.split("|")
                if array[name_col] == spec_name.replace(" ",""):
                    if len(array[wl_col]) > 3:
                        emission_line['wl'] = float(array[wl_col])

                        substring = array[int_col].replace(" ","")
                        substring = substring.replace("q","") # important for Cr
                        if len(substring) > 0:
                            # Col sometimes contain strings, so...
                            try:
                                emission_line['rel_int'] = float(substring)
                            except:
                                emission_line['rel_int'] = -1
                        else:
                            emission_line['rel_int'] = -1

                        if len(array[aik_col]) > 1:
                            emission_line['Aik'] = float(array[4])

                        if len(array[E_col]) > 1:
                            ele = array[E_col]
                            for c in unwanted_chars:
                                ele = ele.replace(c,"")

                            emission_line['Ei'] = float(ele.split('-')[1])
                            emission_line['Ek'] = float(ele.split('-')[0])
                            # NIST does it the other way around
                        emission_lines.append(emission_line)

        return emission_lines
