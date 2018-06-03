#!/bin/python3
import numpy as np
from scipy.signal import fftconvolve as convolve
from scipy import constants as const
from scipy import interpolate
#from . import transition

__all__ = [
    "interpol",
    "doppler_maxwell",
    "gauss_function",
    "doppler_thompson",
    "thompson_function",
    "thompson_side",
    "m_j",
    "zeeman",
    "get_spectroscopic_name",
    "parse_spectroscopic_name",
    ]

def interpol(new_x, y, old_x):
    f = interpolate.interp1d(new_x, y, fill_value=0.0, bounds_error=False)
    return f(old_x)

def doppler_maxwell(x, xc, T, m):
    w = xc/const.c * np.sqrt(8 * const.k * T * np.log(2) / (m*const.u))
    return gauss_function(x, xc, w)

def gauss_function(x, xc, w):
    "A is Area, w is FWHM, H = y0 + A/w / np.sqrt(0.25 * np.pi / np.log(2))"
    return (1/(w*np.sqrt(np.pi/(4*np.log(2)))) *
        np.exp(-4*np.log(2)*(x-xc)**2/w**2)) # gaussian function copied from origin, w is FWHM

def doppler_thompson(x, xc, Eb, m, side=False):
    vb = np.sqrt(2 * Eb * const.eV / (m*const.u))
    if side:
        return thompson_side(x, xc, vb)
    else:
        return thompson_function(x, xc, vb)

def thompson_function(x, xc, vb):
    "Thomspon function derived by myselfe, A is Area; H = y0 + 0.5 A"
    c = const.c
    lx = x[x<=xc]
    rx = x[x>xc]
    profile = np.append((16/np.pi * (xc*vb/c)**3 * ((lx-xc)**2 / ((lx-xc)**2 + (xc/c * vb)**2)**3)),rx*0)
    return profile

def thompson_side(x, xc, vb):
    """Thomspon function derived by e.g. Motohashi et al (Physica Scripta. Vol. T73, 329-331, 1997)
    A is Area. Normed by myselfe"""
    c = const.c
    norm_factor = 0.5 * (xc * vb / c)**2
    return norm_factor*( (x-xc)**2 + (vb*xc/c)**2 )**(-3/2)

def m_j(j):
    if j == round(j,0):
        return np.append(-np.arange(j+1)[::-1][:-1],np.arange(j+1))
    if j == round(j,0) + 0.5 or  j == round(j,0) - 0.5:
        return np.append(-np.arange(j+0.5)[::-1][:-1],np.arange(j+1.5))-0.5


def zeeman(x, cwl, B, upperJ, lowerJ, upperG, lowerG, side=False):
    bm =  const.physical_constants['Bohr magneton'][0]
    z_peaks = []
    y = np.zeros(len(x))
    J_u = upperJ
    J_l = lowerJ
    counter = 0 # for intensity scaling
    for mu in m_j(J_u):
        for ml in m_j(J_l):
            if abs(mu-ml) <= 1: # transition optical allowed?
                counter = counter + 1
                dE = bm * B * mu * upperG - bm * B * ml * lowerG
                E0 = const.h * const.c / (cwl * 1e-9)
                if dE != 0:
                    wl = cwl * (dE/E0 + 1)
                else:
                    wl = cwl
                # intensity acording to Condon,Shortley Page 387
                if J_u == J_l: # J -> J transition
                    if mu == ml: # pi transition with delta M = 0
                        intensity = mu**2
                        if side:
                            intensity = 0
                    else: # sigma transistion with delta M = +-1
                        intensity = 0.25 * (J_u + mu) * (J_u - mu + 1)
                if J_u == (J_l - 1): # J -> J + 1 transition
                    if mu == ml:
                        intensity = (J_u + 1)**2 - mu**2
                        if side:
                            intensity = 0
                    if mu == (ml + 1): # M -> M - 1
                        intensity = 0.25 * (J_u - mu + 1) * (J_u - mu + 2)
                    if mu == (ml - 1): # M -> M + 1
                        intensity = 0.25 * (J_u + mu + 1) * (J_u + mu + 2)
                if J_u == (J_l + 1): # J -> J -1 transition
                    if mu == ml:
                        intensity = J_u**2 - mu**2
                        if side:
                            intensity = 0
                    if mu == (ml + 1): # M -> M - 1
                        intensity = 0.25 * (J_u + mu) * (J_u + mu - 1)
                    if mu == (ml - 1): # M -> M + 1
                        intensity = 0.25 * (J_u - mu) * (J_u - mu - 1)
                index = np.argmin(np.abs(x-wl)) # index of line component
                y[index] = y[index] + intensity # in case two transitions hit the same pixel
    y = y/counter
    return y

def get_spectroscopic_name(element, charge):
    roman = ["I","II","III","IV","V","VI"]
    name = element.title() # ti -> Ti
    designation = roman[charge]
    return name + " " + designation

def parse_spectroscopic_name(name):
    roman = ["I","II","III","IV","V","VI"][::-1]
    for i,num in enumerate(roman):
        if num in name:
            name = name.replace(num, "").strip()
            return name.title(),(len(roman)-i-1)

