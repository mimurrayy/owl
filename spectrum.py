#!/bin/python3
import numpy as np
from scipy import constants as const
from .util import *
import os
from .emitter import transition
import mendeleev 
from astroquery.nist import Nist
from astroquery.nist_levels import NistLevels
import astropy.units as u
import re
import warnings


class spectrum():
    """ Holds information about a complete emission spectrum of the given
    particles in the wavelength range defined by <wl_range> or the
    owl.spectrometer object.

    particles: list of owl.emitter.particle objects
    spectrometer: owl.spectrometer object
    wl_range: tuple or list of lower and upper wavelength of the spectrum"""

    def __init__(self, emitter_name, spectrometer=None, wl_range=None):
        self.name, self.charge = parse_spectroscopic_name(emitter_name)
        self.spec_name = emitter_name
        self.emitter = self.particle = mendeleev.element(self.name)
        self.emitter.charge = self.charge
        self.emitter.m = self.emitter.mass
        self.emitter.Ei = self.emitter.ionenergies[1]
        self.charge = self.emitter.charge
        
        self.spectrometer = spectrometer
        self.wl_range = wl_range
        self.linedata = None
        if spectrometer and not wl_range:
            if spectrometer.x is not None:
                self.wl_range = (spectrometer.x[0],spectrometer.x[-1])


    def load_nist_lines(self, emitter_name, min_int=-1, min_Aik=-1):
        emission_lines = []
        folder = os.path.join(os.path.dirname(os.path.realpath(__file__)),
            "emitter/nist-db")

        unwanted_chars = ["q","[","]","†","?"]

        lines_file = os.path.join(folder, (self.name.lower() + "-lines.txt"))
        spec_name = emitter_name
        name_col = 0
        wl_col = 1
        int_col = 3
        aik_col = 4
        E_col = 6
        g_col = 13
        for line in open(lines_file, 'r').readlines():
            if "Unc." in line: # some NIST files contain WL uncertainties.
                int_col = int_col + 2
                aik_col = aik_col + 2
                E_col = E_col + 2
                g_col = g_col + 2

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
                        substring = substring.replace("r","") # important for Ag
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
                            
                        if len(array[g_col]) > 1:
                            g = array[g_col]
                            emission_line['gi'] = float(g.split('-')[1])
                            emission_line['gk'] = float(g.split('-')[0])

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


    def get_linedata(self):
        if self.wl_range is not None:
            wl_low = self.wl_range[0]
            wl_high = self.wl_range[-1]
        else:
            wl_low = 0
            wl_high = 9999
            
        nist_lines = Nist.query(wl_low*u.nm, wl_high*u.nm, 
                                linename=self.spec_name, wavelength_type='vac+air')

        return nist_lines


    def get_spectrum(self,x=None, width=0.02, mu=0, min_int=-1, min_Aik=-1):
        """ Return simulated spectrum. Lines are pseudo Voigt with
        the set width (FWHM) in nm and form parameter mu (0 = Gauss)  """
        try:
            if not x:
                x = self.spectrometer.x
        except:
            x = x
            
        nist_lines = self.get_linedata()
        
        spectrum = np.zeros(len(x))
        for line in  nist_lines:   
            rel_int = line['Rel.']
            wl =  line['Observed']
            if not np.ma.is_masked(rel_int) and not np.ma.is_masked(wl):
                rel_int = str(rel_int) # is byte type initally
                wl =  str(wl)
                rel_int = re.sub(r'[^\d.]+', '', rel_int)
                wl = re.sub(r'[^\d.]+', '', wl) # removes non-number chars
                try: # conversion to float can still fail...
                    rel_int = float(rel_int)
                    wl = float(wl)
                    if min_Aik > 0:
                        with warnings.catch_warnings(): # ignore masked element warning
                            warnings.simplefilter("ignore", category=UserWarning)
                            Aik =  float(line['Aki'])
                    else:
                        Aik = 1
                        
                    if rel_int > min_int and Aik > min_Aik:
                        profile = rel_int*psd_voigt_function(x, wl, width, mu)
                        spectrum = spectrum + profile
                except:
                    pass

        return spectrum


    def get_LTE_spectrum(self, x, Te, width=0.1, mu=0.5, norm=False, min_int=-1, min_Aik=-1):
        """ Return simulated spectrum with LTE line intnsities in units 
        proportional (!) to Photons/s (NOT W/cm²s)
        Lines are pseudo Voigt with the set width (FWHM) in nm and form 
        parameter mu (0 = Gauss).
        Set norm=True to normalize the intensity to 1 for easier fitting.
        """
        try:
            if not x:
                x = self.spectrometer.x
        except:
            x = x
        nist_lines = self.get_linedata()

        spectrum = np.zeros(len(x))
        for line in  nist_lines:   
            c1 = str(line['Observed'])
            # astroquery does not filter out headings in the middle of the table
            if 'Observed' in c1 or "Wavelength" in c1 or "nm" in c1:
                continue
            
            wl =  line['Observed']
            Aik = line['Aki']
            
            if not np.ma.is_masked(wl) and not np.ma.is_masked(Aik):
                
                wl =  str(wl)
                Aik = str(Aik)
                wl = re.sub(r'[^\d.]+', '', wl) # removes non-number chars
                Aik = re.sub(r'[^\d.]+', '', Aik) # removes non-number chars
                
                if min_int > 0:
                    with warnings.catch_warnings(): # ignore masked element warning
                        warnings.simplefilter("ignore", category=UserWarning)
                        rel_int =  float(line['Rel.'])
                else:
                    rel_int = 1

                try: # conversion to float can still fail...
                    wl = float(wl)
                    Aik = float(Aik)
                    gi = line['gi   gk'].split('-')[1]
                    Ei = line['Ei           Ek'].split('-')[1]
                    Ei = float(Ei)
                    gi = float(gi)

                    if rel_int > min_int and Aik > min_Aik:
                        intensity = Aik*gi*np.exp(-Ei/Te)
                        profile = intensity*psd_voigt_function(x, wl, width, mu)
                        spectrum = spectrum + profile
                except:
                    pass
        if norm == True:
            spectrum = spectrum/np.max(spectrum)
        return spectrum
        

    # def get_transitions(self, debug=False):
    #     """ Returns transition objects in the wavelength range """
    #     transitions = []
    #     if isinstance(self.particles,(list,tuple)):
    #         for particle in self.particles:
    #             for line in self.load_nist_lines(particle):
    #                 if line['wl'] > self.wl_range[0] and line['wl'] < self.wl_range[-1]:
    #                     this_transition = transition(particle, line['wl'], debug=debug)
    #                     transitions.append(this_transition)
    #     else:
    #         for line in self.load_nist_lines(self.particles):
    #             if line['wl'] > self.wl_range[0] and line['wl'] < self.wl_range[-1]:
    #                 this_transition = transition(self.particles, line['wl'], debug=debug)
    #                 transitions.append(this_transition)
    #     return transitions
