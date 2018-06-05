#!/bin/python3
import numpy as np
from scipy.signal import fftconvolve as convolve
from scipy import constants as const
from scipy import interpolate
#from . import emitter
from .util import *
from . import stark

class emission_line():
    def __init__(self, spectrometer, transition, wl,
    T = None, Eb = None, gamma = None, B = None, ne = None, Te = None, P = None):
        """ T in K, Eb,Te in eV, ne in m^-3"""
        self.wl = wl
        self.transition = transition
        particle = transition.particle
        self.m  = particle.m
        self.Ei = particle.Ei
        self.T  = T
        self.Eb = Eb
        self.gamma = gamma
        self.B  = B
        self.ne = ne
        self.P  = P
        self.spectrometer = spectrometer
        self.side = False # symmetric version of Thompson?
        self.last_profiles = []


    def get_profile(self, x, A,
            wl=None, T = None, Eb = None, gamma = None, B = None, P = None,
            ne = None, Te = None, ion = None):

        if wl == None and self.wl:
            wl = self.wl
        if T == None and self.T:
            T = self.T
        if Eb == None and self.Eb:
            Eb = self.Eb
        if gamma == None and self.gamma:
            gamma = self.gamma
        if B == None and self.B:
            B = self.B
        if ne == None and self.ne:
            ne = self.ne
        if P == None and self.P:
            P = self.P

        self.last_profiles = []
        resolution = abs((x[-1]-x[0])/(len(x)-1))
        orig_x = x
        x = self.spectrometer.fine_x
        middle_wl = x[int(len(x)/2)]

        # We let the instrumental profile determine center position.
        # ALL other components are shifted to the middle!
        instrumental_profile = self.spectrometer.instrument_function(x, wl, self.transition)
        instrumental_profile = instrumental_profile/resolution/1024
        components = []

        if T:
            if Eb and gamma:
                thompson = doppler_thompson(x, middle_wl, Eb, self.m, self.side)
                maxwell = doppler_maxwell(x, middle_wl, T, self.m)
                doppler = gamma*maxwell + (1-gamma)*thompson
                self.last_profiles.append(gamma*maxwell)
                self.last_profiles.append((1-gamma)*thompson)
                self.last_profiles.append(doppler)
                components.append(doppler)
            else:
                maxwell = doppler_maxwell(x, middle_wl, T, self.m)
                self.last_profiles.append(maxwell)
                components.append(maxwell)
        if Eb and not T:
            thompson = doppler_thompson(x, middle_wl, Eb, self.m, self.side)
            self.last_profiles.append(thompson)
            components.append(thompson)
        if B:
            t = self.transition
            zeeman_pattern = zeeman(x, middle_wl, B, t.upperJ, t.lowerJ, t.upperG, t.lowerG)/resolution
            self.last_profiles.append(zeeman_pattern)
            components.append(zeeman_pattern)

        if ne:
            this_stark = stark.stark(self.transition)
            stark_profile = this_stark.get_profile(x, ne, Te, ion)
            self.last_profiles.append(stark_profile)
            components.append(stark_profile)

        if P:
            print("TODO")

        profile = instrumental_profile
        for component in components:
            profile = convolve(profile, component, mode='same') * resolution
        y = A*profile
        self.last_profiles.append(y)
        lowres_y = interpol(x, y, orig_x)
        self.last_profiles.append(lowres_y)
        return lowres_y


    def get_instrumental_profile(self, x, wl=None):
        if wl == None and self.wl:
            wl = self.wl
        return self.spectrometer.instrument_function(x, wl, self.transition)