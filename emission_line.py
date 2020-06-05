#!/bin/python3
import numpy as np
from scipy.signal import fftconvolve as convolve
from scipy import constants as const
from scipy import interpolate
#from . import emitter
from .util import *
from . import stark
from . import vdW

class emission_line():
    def __init__(self, spectrometer, transition, wl,
    T = None, Eb = None, gamma = None, B = None, ne = None, Te = None,
    side=False, N = None, pert = None):
        """ T in K, Eb,Te in eV, ne in m^-3, N in m^-3"""
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
        self.N  = N
        self.pert = pert
        self.spectrometer = spectrometer
        self.side = side # symmetric version of Thompson?
        self.last_profiles = []


    def get_profile(self, x, A = 1,
            shift = False, wl=None, T = None, Eb = None, gamma = None, B = None,
            N = None, ne = None, Te = None, pert = None):

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
        if N == None and self.N:
            N = self.N
        if pert == None and self.pert:
            pert = self.pert

        self.last_profiles = []
        orig_x = x
        s = 0 # lineshift in nm
        try:
            x = self.spectrometer.fine_x
            if len(orig_x)<1024:
                x = x[x>orig_x[0]]
                x = x[x<orig_x[-1]]
        except:
            x = x

        resolution = abs((x[-1]-x[0])/(len(x)-1))

        middle_wl = x[int(len(x)/2)] - 0.5*abs((x[-1]-x[0])/(len(x)-1))
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
            stark_profile = this_stark.get_profile(x, ne, Te, pert)
            self.last_profiles.append(stark_profile)
            components.append(stark_profile)
        if N:
            this_vdW = vdW.vdW(self.transition, pert=pert)
            vdW_profile = this_vdW.get_profile(x, T, N)
            self.last_profiles.append(vdW_profile)
            components.append(vdW_profile)
            if shift:
                s = s + this_vdW.get_shift(x, N, T)

        # We let the instrumental profile determine center position.
        # ALL other components are shifted to the middle!
        instrumental_profile = self.spectrometer.instrument_function(x, wl+s, self.transition)
        instrumental_profile = instrumental_profile/resolution/1024
        profile = instrumental_profile

        for component in components:
            profile = convolve(profile, component, mode='same') * resolution
        y = A*profile
        self.last_profiles.append(y)
        lowres_y = interpol(x, y, orig_x)
        self.last_profiles.append(lowres_y)
        # normalize intensity to A
        y = A*lowres_y/np.sum(lowres_y)/resolution
        return y


    def get_instrumental_profile(self, x, wl=None):
        if wl == None and self.wl:
            wl = self.wl
        return self.spectrometer.instrument_function(x, wl, self.transition)
