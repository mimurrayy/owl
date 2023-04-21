#!/usr/bin/python

import numpy as np
from ..util import *
from .. import emitter
import mendeleev 
from scipy import constants as const

class vdW():
    def __init__(self, transition, pert):
        self.transition = transition
        self.pert = pert


    def get_profile(self,x, T, n):
        middle_wl = x[int(len(x)/2)] - 0.5*abs((x[-1]-x[0])/(len(x)-1))
        if not T:
            print("""Need to provide a gas temperture (in Kelvin).
            Assuming 300 K""")
            T = 300

        if self.transition.emitter.symbol == "H":
            return self.hydrogen_profile(x,T,n)
        w = self.get_width(self.transition, n, T)
        y = lorentz_function(x,middle_wl,w)

        return y/max(y)


    def get_width(self, transition, n, T):
        """
        Interface for get_width2 using the transition properties.
        Energies in eV, alpha in m^2,
        lambda in nm, n in m^-3, T in K, m in amu.
        """
        
        xc = self.transition.wl
        m1 = self.transition.particle.m
        m2 = self.pert.m
        Eion = self.transition.particle.Ei
        Eu = self.transition.upperE
        El = self.transition.lowerE
        lu = self.transition.upperl
        ll = self.transition.lowerl
        alpha = self.pert.dipole_polarizability*const.value('Bohr radius')**3
                
        return self.get_width2(xc,m1,m2,T,n,Eion,Eu,El,lu,ll,alpha)


    def get_shift(self,x, n, T):
        w = self.get_width(self.transition, n, T)
        s = w*0.28 # citation needed
        return s

    
    def get_width2(self,xc,m1,m2,T,n,Eion,Eu,El,lu,ll,alpha):
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
        alpha = alpha*1e6 # to cm^3
        n = n/1e6 # to cm^-3
        xc = xc/1e7 # to cm
        w = 8.18e-12 * (xc**2 ) * ((alpha * R**2)**(2/5)) * ((T/mu)**(3/10))*n
        return w*1e7 # to nm


    def hydrogen_profile(self,x,T,n):
        """ Data from NIST: J. Phys. Chem. Ref. Data, Vol. 38, No. 3, 2009"""
        middle_wl = x[int(len(x)/2)] - 0.5*abs((x[-1]-x[0])/(len(x)-1))
        xc = self.transition.wl
        m1 = self.transition.particle.m
        m2 = self.pert.m
        Eion = self.transition.particle.Ei
        a = self.pert.dipole_polarizability*const.value('Bohr radius')**3

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
                w = self.get_width2(xc,m1,m2,T,n,Eion,*comp[2:],a)
                y = y + Aik*lorentz_function(x,wl-656.280+middle_wl,w)


        if round(self.transition.wl, 0) == 486:
            H = self.transition.particle
            #[Aik*gi, wl, Eu, El, lu, ll]
            components = [
            [1.7188e+07*4, 486.1278624, 12.74853800, 10.19880606, 2, 1],
            [2.0625e+07*6, 486.1361516, 12.74853989, 10.19885143, 2, 1],
            [3.4375e+06*4, 486.1365118, 12.74853800, 10.19885143, 2, 1],

            [9.6680e+06*4, 486.1286949, 12.74853801, 10.19881044, 1, 0],
            [9.6683e+06*2, 486.1297761, 12.74853234, 10.19881044, 1, 0],

            [8.5941e+05*2, 486.1288370, 12.74853289, 10.19880606, 0, 1],
            [1.7190e+06*2, 486.1374864, 12.74853289, 10.19885143, 0, 1]]

            y = np.zeros(len(x))
            for comp in components:
                Aik = comp[0]
                wl = comp[1]
                w = self.get_width2(xc,m1,m2,T,n,Eion,*comp[2:],a)
                y = y + Aik*lorentz_function(x,wl-486.133+middle_wl,w)


        return y/max(y)
