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

        self.nist_info()


    def nist_info(self):
        "Load NIST Tables and set E,J and G"
        folder = os.path.join(os.path.dirname(os.path.realpath(__file__)), "nist-db")
        lines_file = os.path.join(folder, (self.element.lower() + "-lines.txt"))
        levels_file = os.path.join(folder, (self.element.lower() + str(self.charge+1) + "-levels.txt"))
        spec_name = get_spectroscopic_name(self.element, self.charge)
        for line in open(lines_file, 'r').readlines():
            if line.startswith(spec_name):
                line = line.replace(" ", "")
                array = line.split("|")
                if len(array[1]) > 3:
                    observed_wl = float(array[1])
                    if self.wl == round(observed_wl, 3):
                        self.upperE = float(array[6].split("-")[1])
                        self.lowerE = float(array[6].split("-")[0])

        for line in open(levels_file, 'r').readlines():
            if str(self.upperE) in line:
                line = line.replace(" ", "")
                array = line.split("|")
                self.upperJ = float(eval(array[2]))
                self.upperG = float(array[4])
            if str(self.lowerE) in line:
                line = line.replace(" ", "")
                array = line.split("|")
                self.lowerJ = float(eval(array[2]))
                self.lowerG = float(array[4])

    def get_instrumental_params(self):
        "Database for instrumental functions for lines."
        if self.element == "Ti" and (self.wl == 453.324 or self.wl == 453.478 or self.wl == 453.396):
            return (0.00653, 0.00173, 0.7795, 1.0, 1.0295, 0.005, 0.0115)
        if self.element == "Ti" and self.charge == 0 and (self.wl == 399.864 or self.wl == 398.976):
            return (0.00653, 0.00173, 0.7795, 1.0, 1.0295, 0.005, 0.0115)
        if self.element == "Ar" and self.charge == 1 and (self.wl == 472.687):
            return (0.00416, 0.00322, 1.2025, 0.0, 2.5765, 0.0029, 0.0139)





