#!/bin/python3
import numpy as np
from scipy import constants as const
from .util import *
import os
from .emitter import transition

class spectrum():
    """ Holds information about a complete emission spectrum of the given
    particles in the wavelength range defined by <wl_range> or the
    owl.spectrometer object.

    particles: list of owl.emitter.particle objects
    spectrometer: owl.spectrometer object
    wl_range: tuple or list of lower and upper wavelength of the spectrum"""

    def __init__(self, particles, spectrometer=None, wl_range=None):
        self.particles = particles
        self.spectrometer = spectrometer
        self.wl_range = wl_range
        self.linedata = None
        if spectrometer and not wl_range:
            if spectrometer.x is not None:
                self.wl_range = (spectrometer.x[0],spectrometer.x[-1])

    def load_nist_lines(self, particle, min_int=-1, min_Aik=-1):
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
                        substring = substring.replace("g","") # important for Al
                        if len(substring) > 0:
                            # Col sometimes contain strings, so...
                            try:
                                emission_line['rel_int'] = float(substring)
                            except:
                                emission_line['rel_int'] = -1
                        else:
                            emission_line['rel_int'] = -1

                        if len(array[aik_col]) > 1:
                            emission_line['Aik'] = float(array[aik_col])
                        else:
                            emission_line['Aik'] = -1

                        if len(array[E_col]) > 1:
                            ele = array[E_col]
                            for c in unwanted_chars:
                                ele = ele.replace(c,"")

                            emission_line['Ei'] = float(ele.split('-')[1])
                            emission_line['Ek'] = float(ele.split('-')[0])
                            # NIST does it the other way around
                            
                        if (emission_line['rel_int'] >= min_int)\
                                         and (emission_line['Aik'] >= min_Aik):
                            emission_lines.append(emission_line)
        return emission_lines


    def get_linedata(self, min_int=-1, min_Aik=-1):
        lines = []
        if isinstance(self.particles,(list,tuple)):
            for particle in self.particles:
                lines += self.load_nist_lines(particle, min_int, min_Aik)
        else:
            lines = self.load_nist_lines(self.particles, min_int, min_Aik)

        lines = sorted(lines, key=lambda k: k['wl'])
        if self.wl_range:
            lines = list(filter(lambda k: k['wl']>self.wl_range[0], lines))
            lines = list(filter(lambda k: k['wl']<self.wl_range[-1], lines))

        return lines


    def get_spectrum(self,x=None, width=0.02):
        """ Return simulated spectrum. Lines are Gaussian with
        the set width (FWHM) in nm  """
        try:
            if not x:
                x = self.spectrometer.x
        except:
            x = x
        lines = self.get_linedata()
        lines = list(filter(lambda k: k['wl']>self.wl_range[0], lines))
        lines = list(filter(lambda k: k['wl']<self.wl_range[-1], lines))
        spectrum = np.zeros(len(x))
        for line in lines:
            profile = line['rel_int']*gauss_function(x, line['wl'], width)
            spectrum = spectrum + profile

        return spectrum

    def get_transitions(self, debug=False):
        """ Returns transition objects in the wavelength range """
        transitions = []
        if isinstance(self.particles,(list,tuple)):
            for particle in self.particles:
                for line in self.load_nist_lines(particle):
                    if line['wl'] > self.wl_range[0] and line['wl'] < self.wl_range[-1]:
                        this_transition = transition(particle, line['wl'], debug=debug)
                        transitions.append(this_transition)
        else:
            for line in self.load_nist_lines(self.particles):
                if line['wl'] > self.wl_range[0] and line['wl'] < self.wl_range[-1]:
                    this_transition = transition(self.particles, line['wl'], debug=debug)
                    transitions.append(this_transition)
        return transitions
