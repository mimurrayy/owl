#!/usr/bin/python

import numpy as np
from ..util import *
from .. import emitter

class vdW():
    def __init__(self, transition, pert):
        self.transition = transition
        self.pert = pert

    def get_profile(self,x, T, n):
        middle_wl = x[int(len(x)/2)] - 0.5*abs((x[-1]-x[0])/(len(x)-1))
        if not T:
            print("Need to provide a gas temperture (in Kelvin).")

        if self.transition.particle.element == "H":
            return self.hydrogen_profile(x,T,n)

        xc = self.transition.wl
        m1 = self.transition.particle.m
        m2 = self.pert.m
        Eion = self.transition.particle.Ei
        Eu = self.transition.upperE
        El = self.transition.lowerE
        lu = self.transition.upperl
        ll = self.transition.lowerl
        a = self.alpha(self.pert)

        w = self.get_width(xc,m1,m2,T,n,Eion,Eu,El,lu,ll,a)
        y = lorentz_function(x,middle_wl,w)

        return y/max(y)


    def get_width(self,xc,m1,m2,T,n,Eion,Eu,El,lu,ll,alpha):
        """
        Energies in eV, alpha in m^2,
        lambda in nm, n in m^-3, T in K, m in amu.
        From: N. Konjevic & / Physics Reports 316 (1999) 339}401
        """
        mu = m1*m2/(m1+m2)

        ### upper
        EH = 13.59844
        nj = np.sqrt(EH/(Eion-Eu))
        Ru2 = 0.5*(nj**2) * (5 * (nj**2) + 1 - 3*lu * (lu + 1))

        ### lower
        nj = np.sqrt(EH/(Eion-El))
        Rl2 = 0.5*(nj**2) * (5 * (nj**2) + 1 - 3*ll * (ll + 1))

        R = np.sqrt(Ru2 - Rl2)
        alpha = alpha/1e6 # to cm^3
        n = n/1e6 # to cm^-3
        xc = xc/1e7 # to cm
        w = 8.18e-12 * (xc**2 ) * ((alpha * R**2)**(2/5)) * ((T/mu)**(3/10))*n
        return w*1e7 # to nm


    def get_shift(self,x, T, n):
        xc = self.transition.wl
        m1 = self.transition.particle.m
        m2 = self.pert.m
        Eion = self.transition.particle.Ei
        Eu = self.transition.upperE
        El = self.transition.lowerE
        lu = self.transition.upperl
        ll = self.transition.lowerl
        a = self.alpha(self.pert)

        w = self.get_width(xc,m1,m2,T,n,Eion,Eu,El,lu,ll,a)
        s = w*0.28 # citation needed
        return s

    def hydrogen_profile(self,x,T,n):
        """ Data from NIST: J. Phys. Chem. Ref. Data, Vol. 38, No. 3, 2009"""
        middle_wl = x[int(len(x)/2)] - 0.5*abs((x[-1]-x[0])/(len(x)-1))
        xc = self.transition.wl
        m1 = self.transition.particle.m
        m2 = self.pert.m
        Eion = self.transition.particle.Ei
        a = self.alpha(self.pert)

        if round(self.transition.wl, 0) == 656:
            H = self.transition.particle
            #[Aik*gi, wl, Eu, El, lu, ll]
            components = [
            [2.2448e-01*4, 656.2724, 12.087507, 10.19881, 1, 0],
            [2.2449e-01*2, 656.2771, 12.08749, 10.19881, 1, 0],

            [4.2097e-02*2, 656.2909, 12.087495, 10.19885, 0, 1],
            [2.1046e-02*2, 656.2752, 12.087495, 10.198806, 0, 1],

            [6.4651e-01*6, 656.2852, 12.08751, 10.19885, 0, 1],
            [5.3877e-01*4, 656.2701, 12.087507, 10.198806, 2, 1],
            [1.0775e-01*4, 656.2868, 12.087507, 10.19885, 2, 1]]

            y = np.zeros(len(x))
            for comp in components:
                Aik = comp[0]
                wl = comp[1]
                w = self.get_width(xc,m1,m2,T,n,Eion,*comp[2:],a)
                y = y + Aik*lorentz_function(x,wl-656.280+middle_wl,w)

        return y/max(y)

    def alpha(self,particle):
        """ Returns polarizability of partilce in m^3
        All Data from CRC Handbook of Chemistry and Physics """
        if particle.name == "H2O":
            return(1.45e-18)
        if particle.name == "H":
            return(0.666793e-18)
        if particle.name == "He":
            return(0.204956e-18)
        if particle.name == "O":
            return(0.802e-18)
        if particle.name == "Ar":
            return(1.6411e-18)
        if particle.name == "O2":
            return(1.5812e-18)
        if particle.name == "N2":
            return(1.7403e-18)