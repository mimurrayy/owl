#!/bin/python3
import numpy as np
from scipy.signal import fftconvolve as convolve
from scipy import constants as const
from scipy import interpolate
from .util import zeeman, parse_spectroscopic_name, doppler_thompson,\
                        doppler_maxwell, gauss_function, psd_voigt
from . import stark
from . import vdW
import mendeleev 

class emission_line():
    def __init__(self, transition, wl, instr_func = None, w = None, mu = 0.0,
                 T = None, Eb = None, gamma = None, B = None, ne = None, 
                                Te = None, side = False, N = None, pert = None):
        
        """  Class to calculate the shape of an emission line based on the 
        physical situation. Providing information on densities and temperatures
        automatically switches on different broadening mechanisms.

        Supported broadening mechanisms: Doppler, Stark, Zeeman, van der Waals
        and instrumental broadening. Stark broadening is only supported for
        some transitions: H-alpha to gamma, He 447.1 nm, He 492.2 nm, O 777 nm,
        Ar 810.369 nm and Ar 738.398 nm.

        <transition>: An owl.emitter.transition object for the specific line.
        <wl>: Center wavelength in nm. Does not need to be the theoretical wl
              of the transition. Use to move simulation relative to measurement.
        <instr_func>: Instrumental function of the spectrometer. Must be a 
                      function that take parameters x and xc (wavelength axis 
                      and central position) and return the profile as a 
                      numpy array.
        <w>: Width of the instrumental function in nm. Only if instr_func=None.
        <mu>: Shape parameter of the instrumental fucntion. 0=Gauss, 1=Lorentz.
        Use together with w or use instr_func instead. Only if instr_func=None.
        <T>: Emitter temperature in Kelvin. Used for vdW and Doppler broadening.
        <Eb>: Surface binding energy in eV to calculate Doppler broadening of 
              sputtered particles.
        <gamma>: Thermalization degree of sputtered particles (used with Eb).
        <B>: Magnetic field strength for Zeeman splitting. Units of Tesla.
        <ne>: Electron density for Stark broadening. Units of m^-3. 
        <Te>: Electron temperature for stark broadening. Units of Kelvin.
              Has only minor impact on the results and can usually be guessed.
        <side>: Used for Zeeman and together with gamma. To be updated.
        <N>: Neutral density for van der Waals broadening. Units of m^-3.
             This is the density of species sourrounding the emitter.
        <pert>: Perturber species for vdW and Stark boradening. This is the ion
                or atom species surrounding the emitter particles. In case vdW 
                we look up the polarizability of the element. For Stark we only 
                need the mass. Specify like emitter species, e.g. 'Ar I'.  
        """

        self.wl = wl
        self.transition = transition
        self.instr = instr_func
        self.w = w
        self.mu = mu
        particle = transition.particle
        self.m  = particle.m
        self.Ei = particle.Ei
        self.T  = T
        self.Eb = Eb
        self.gamma = gamma
        self.B  = B
        self.ne = ne
        self.Te = Te
        self.N  = N
        self.pert = pert
        self.side = side # symmetric version of Thompson?
        self.profiles = dict.fromkeys(['Doppler', 'Zeeman', 'Stark', 'vdW', 'instrument'])

        if pert:          
            self.pert_name, self.pert_charge = parse_spectroscopic_name(pert)
            self.pert = mendeleev.element(self.pert_name)
            self.pert.charge = self.pert_charge
            self.pert.m = self.pert.mass
            self.pert.Ei = self.pert.ionenergies[1]
            if T:
                self.pert.T = T
            elif Te:
                self.pert.T = Te*const.eV/const.k



    def get_profile(self, x, A = 1, instr_func = None, w = None, mu = None,
            shift = False, wl=None, T = None, Eb = None, gamma = None, B = None,
            N = None, ne = None, Te = None, pert = None):

        if instr_func == None and self.instr:
            instr_func = self.instr
        if w == None and self.w:
            w = self.w
        if mu == None and self.mu:
            mu = self.mu
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
        if Te == None and self.Te:
            Te = self.Te
        if N == None and self.N:
            N = self.N
        if pert == None and self.pert:
            pert = self.pert
        
        elif pert:
            self.pert_name, self.pert_charge = parse_spectroscopic_name(pert)
            self.pert = mendeleev.element(self.pert_name)
            self.pert.charge = self.pert_charge
            self.pert.m = self.pert.mass
            self.pert.Ei = self.pert.ionenergies[1]
            if T:
                self.pert.T = T
            elif Te:
                self.pert.T = Te*const.eV/const.k
            pert = self.pert

        self.profiles = dict.fromkeys(['x', 'y', 'Doppler', 'Zeeman', 'Stark', 'vdW', 'instrument'])
        orig_x = np.copy(x)
        s = 0 # lineshift in nm
        fine_x = interpolate.make_interp_spline(np.linspace(0,1,len(x)), x)
        x = fine_x(np.linspace(0,1,len(x)*50)) # upsample x axis
        self.profiles['x'] = x

        resolution = abs((x[-1]-x[0])/(len(x)-1))
        middle_wl = x[int(len(x)/2)-1]
        components = []

        if T:
            if Eb and gamma:
                thompson = doppler_thompson(x, middle_wl, Eb, self.m, self.side)
                maxwell = doppler_maxwell(x, middle_wl, T, self.m)
                doppler = gamma*maxwell + (1-gamma)*thompson
                self.profiles['Doppler'] = doppler/np.max(doppler)
                components.append(doppler)
            else:
                maxwell = doppler_maxwell(x, middle_wl, T, self.m)
                self.profiles['Doppler'] = maxwell/np.max(maxwell)
                components.append(maxwell)
        if Eb and not T:
            thompson = doppler_thompson(x, middle_wl, Eb, self.m, self.side)
            self.profiles['Doppler'] = thompson/np.max(thompson)
            components.append(thompson)
        if B:
            t = self.transition
            zeeman_pattern = zeeman(x, middle_wl, B, t.upperJ, t.lowerJ, t.upperG, t.lowerG)/resolution
            self.profiles['Zeeman'] = zeeman_pattern/np.max(zeeman_pattern)
            components.append(zeeman_pattern)
        if ne:
            this_stark = stark.stark(self.transition)
            stark_profile = this_stark.get_profile(x, ne, Te, pert)
            self.profiles['Stark'] = stark_profile/np.max(stark_profile)
            components.append(stark_profile)
        if N:
            this_vdW = vdW.vdW(self.transition, pert=pert)
            vdW_profile = this_vdW.get_profile(x, T, N)
            self.profiles['vdW'] = vdW_profile/np.max(vdW_profile)
            components.append(vdW_profile)
            if shift:
                s = s + this_vdW.get_shift(x, N, T)

        # We let the instrumental profile determine center position.
        # ALL other components are shifted to the middle!
        if instr_func:
            instrumental_profile = instr_func(x, wl+s)
        elif w:   
            instrumental_profile = psd_voigt(x, wl+s, w, mu)
        else:
            instrumental_profile = gauss_function(x, wl+s, x[1]-x[0])
        
        self.profiles['instrument'] = instrumental_profile/np.max(instrumental_profile)

        profile = instrumental_profile
        for component in components:
            profile = convolve(profile, component, mode='same') * resolution

        y = A*profile
        self.profiles['y'] = y
        lowres_y = np.interp(orig_x, x, y)
        # normalize intensity to A
        y = A*lowres_y/np.sum(lowres_y)/resolution
        return y