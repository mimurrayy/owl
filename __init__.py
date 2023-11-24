#!/usr/bin/python
from . import emitter
from . import spectrometer
from .util import *
from .emission_line import *
from .spectrum import *
from .stark import gigosos_loader, gigosos_he_loader

def download_stark_tables():
    print("Starting download of Stark broadening data tables for hydrogen.")
    gigosos_loader.download_profiles()
    print("Starting download of Stark broadening data tables for helium.")
    gigosos_he_loader.download_profiles()
    print("Done.")
