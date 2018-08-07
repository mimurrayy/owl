#!/bin/python3
import numpy as np
from scipy import constants as const
from .util import *
import os

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

    def load_nist_lines(self,particle):
        emission_lines = []
        folder = os.path.join(os.path.dirname(os.path.realpath(__file__)),
            "emitter/nist-db")

        lines_file = os.path.join(folder, (particle.element.lower() + "-lines.txt"))
       # levels_file = os.path.join(folder, (self.element.lower() + str(self.charge+1) + "-levels.txt"))
        spec_name = particle.spectroscopic_name()

        for line in open(lines_file, 'r').readlines():
            if line.startswith(spec_name):
                emission_line = {}
                emission_line['particle'] = spec_name
                observed_wl = rel_int = Aik = Ei = Ek = None
                line = line.replace(" ", "")
                array = line.split("|")
                if array[0] == spec_name.replace(" ",""):
                    if len(array[1]) > 3:
                        emission_line['wl'] = float(array[1])

                        if len(array[3]) > 1:
                            # Col sometimes contain strings, so...
                            emission_line['rel_int'] = float(''.join(list(
                                filter(str.isdigit, array[3]))))
                        else:
                            emission_line['rel_int'] = -1

                        if len(array[4]) > 1:
                            emission_line['Aik'] = float(array[4])

                        if len(array[6]) > 1:
                            ele = array[6]
                            emission_line['Ei'] = float(ele.split('-')[1])
                            emission_line['Ek'] = float(ele.split('-')[0])
                            # NIST does it the other way around
                        emission_lines.append(emission_line)
        return emission_lines


    def get_linedata(self):
        lines = []
        if isinstance(self.particles,(list,tuple)):
            for particle in self.particles:
                lines += self.load_nist_lines(particle)
        else:
            lines = self.load_nist_lines(self.particles)

        lines = sorted(lines, key=lambda k: k['wl'])

        return lines


    def get_spectrum(self,x):
        """ stub """
        lines = self.get_linedata()
        lines = list(filter(lambda k: k['wl']>self.wl_range[0], lines))
        lines = list(filter(lambda k: k['wl']<self.wl_range[-1], lines))
        return None
