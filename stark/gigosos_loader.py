#!/usr/bin/python3

import numpy as np
from scipy import constants as const
from ..util import interpol
import os
from platformdirs import user_data_dir
import requests
from pathlib import Path
from zipfile import ZipFile
from io import BytesIO
data_folder = os.path.join(user_data_dir("owl-OES", "owl-OES"), "Gigosos2003H")

class gigosos_loader():
    def __init__(self, transition):
        self.transition = transition
        
        
    def download_profiles(redownload = False):
        if os.path.exists(data_folder+"/HalphaProfiles.zip") and redownload==False:
            return
        print("Downloading Stark broadening data tables for hydrogen.")
        URL = "https://ars.els-cdn.com/content/image/1-s2.0-S0584854703000971-mmc1.zip"
        try:
            response = requests.get(URL)
            zip_file = response.content
        except:
            print("ERROR: could not download hydrogen Stark broadening tables")
            
        try:
            Path(data_folder).mkdir(parents=True, exist_ok=True)
            with ZipFile(BytesIO(zip_file)) as zip_ref:
                zip_ref.extractall(data_folder)
        except:
            print("ERROR: could now save hydrogen Stark broadening tables to disk")


    def load(self, ne, Te, pert):
        r0 = (3/(4 * np.pi * ne))**(1/3)
        rD = ((const.epsilon_0 * const.eV * Te)/(ne * const.e**2))**(1/2)
        rho = r0/rD
        if pert:
           emitter_m = self.transition.particle.m
           reduced_m = (emitter_m * pert.m)/(emitter_m + pert.m)
           mu = reduced_m * Te/pert.T
        else:
           mu = 1
        return self.load_stark_profile(ne, mu, rho)


    def closest(self,val,array):
        return array[np.argmin(np.abs(array-val))]


    def slightly_smaller(self,val,array):
        array = array[array < val]
        return array[np.argmin(np.abs(array-val))]


    def slightly_bigger(self,val,array):
        array = array[array > val]
        return array[np.argmin(np.abs(array-val))]


    def mirror_profile(self,profile):
        return np.append(np.array(profile[::-1]), profile)


    def mirror_x(self,x):
        return np.append(np.array(x[::-1])*-1, x)


    def closest_ne(self,ne):
        ne_stock = 10**np.linspace(20,25,16)
        if ne in ne_stock:
            ne_low = ne_high = ne
        else:
            ne_low  = self.slightly_smaller(ne,ne_stock)
            ne_high = self.slightly_bigger(ne,ne_stock)
        return ne_low,ne_high


    def closest_mu(self,mu):
        mu_stock = np.array([0.50, 0.80, 0.90, 1.00, 1.25, 1.5, 1.75, 2.00,
            2.50, 3.00, 4.00, 5.00, 6.00, 7.00, 8.00, 9.00, 10.0])
        if mu in mu_stock:
            mu_low = mu_high = mu
        elif mu < mu_stock[0]:
            mu_low = mu_high = mu_stock[0]
        elif mu > mu_stock[-1]:
            mu_low = mu_high = mu_stock[-1]
        else:
            mu_low  = self.slightly_smaller(mu,mu_stock)
            mu_high = self.slightly_bigger(mu,mu_stock)
        return mu_low,mu_high


    def closest_rho(self,rho):
        rho_stock = np.linspace(0.1,0.6,11)
        if rho in rho_stock:
            rho_low = rho_high = rho
        elif rho < rho_stock[0]:
            rho_low = rho_high = rho_stock[0]
        elif rho > rho_stock[-1]:
            rho_low = rho_high = rho_stock[-1]
        else:
            rho_low  = self.slightly_smaller(rho,rho_stock)
            rho_high = self.slightly_bigger(rho,rho_stock)
        return rho_low,rho_high


    def rho_to_T(self,rho,ne):
        r0 = (3/(4 * np.pi * ne))**(1/3)
        T = ((ne*const.e**2)/(const.epsilon_0 * const.k))*(r0/rho)**2
        T = T*0.9999944 # Gigosos temperature calculations are strange...
        return T


    def load_file(self,ne,mu,rho,prefix,datazip):        
        ne_name = str(int(round(np.log10(ne)*100,0)))
        T = str(int(round(self.rho_to_T(rho,ne),0)))
        T_name = str(T).zfill(7)
        mu_name = str(int(mu*100)).zfill(4)
        filename = prefix + ne_name + "t" + T_name + "m" + mu_name + ".dlp"
        try:
            datafile = datazip.open(filename)
            x,y = np.loadtxt(datafile).T
        except:
            raise SystemExit('Error: Could not find Stark broadening data tables.')
        return x,y


    def interpolate_stark_profile(self,low_x, high_x, low_y, high_y, low_val, high_val, val):
        "Creates profile for arbitrary ne. Arrays must use same x!"
        low_y = interpol(low_x, low_y, high_x) # all to high_x
        if high_val != low_val:
            low_y = np.array(low_y)
            high_y = np.array(high_y)
            mix_factor = ((val - low_val)/(high_val - low_val))#**(1.5)
            return high_x,(mix_factor * high_y + (1-mix_factor) * low_y)
        else:
            return high_x,high_y


    def interpolate_stark_profile_ne(self,low_x, high_x, low_y, high_y, low_val, high_val, val):
        "Creates profile for arbitrary ne. Arrays must use same x!"
        low_y = interpol(low_x, low_y, high_x) # all to high_x
        if high_val != low_val:
            low_y = np.array(low_y)
            high_y = np.array(high_y)
            mix_factor = ((val - low_val)/(high_val - low_val))
            return high_x,(mix_factor * high_y + (1-mix_factor) * low_y)
        else:
            return high_x,high_y


    def load_stark_profile(self,ne_, mu_, rho_):
        
        folder = data_folder
        ele = self.transition.emitter.symbol
        if ele == "H" and round(self.transition.wl,0) == 656.0:
            datazip = os.path.join(folder,"HalphaProfiles.zip")
            prefix = "BAn"
        if ele == "H" and round(self.transition.wl,0) == 486.0:
            datazip = os.path.join(folder,"HbetaProfiles.zip")
            prefix = "BBn"
        if ele == "H" and round(self.transition.wl,0) == 434.0:
            datazip = os.path.join(folder,"HgammaProfiles.zip")
            prefix = "BGn"
        datazip = ZipFile(datazip, 'r')
        
        x = None
        profiles_ne = []
        vals_ne = []
        for ne in self.closest_ne(ne_):
            profiles_mu = []
            vals_mu = []
            for mu in self.closest_mu(mu_):
                rho_low,rho_high = self.closest_rho(rho_)
                low_x,y_low = self.load_file(ne,mu,rho_low,prefix,datazip)
                high_x,y_high = self.load_file(ne,mu,rho_high,prefix,datazip)
                x,y = self.interpolate_stark_profile(low_x,high_x,y_low,y_high,rho_low,rho_high,rho_)
                profiles_mu.append((x,y))
                vals_mu.append(mu)
            x,y = self.interpolate_stark_profile(
                profiles_mu[0][0],profiles_mu[1][0],profiles_mu[0][1],
                profiles_mu[1][1],vals_mu[0],vals_mu[1],mu_)

            profiles_ne.append((x,y))
            vals_ne.append(ne)

        x,y = self.interpolate_stark_profile_ne(
            profiles_ne[0][0],profiles_ne[1][0],profiles_ne[0][1],
            profiles_ne[1][1],vals_ne[0],vals_ne[1],ne_)

        y = self.mirror_profile(y)/1e9
        x = self.mirror_x(x)*1e9

        return x,y



