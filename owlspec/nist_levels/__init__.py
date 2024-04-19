# Licensed under a 3-clause BSD style license - see LICENSE
"""
Fetches level information from the NIST Atomic Spectra Database.
"""
from astropy import config as _config


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astroquery.nist_levels`.
    """
    server = _config.ConfigItem(
        ['https://physics.nist.gov/cgi-bin/ASD/energy1.pl'],
        'Name of the NIST URL to query.')
    timeout = _config.ConfigItem(
        30,
        'Time limit for connecting to NIST server.')


conf = Conf()

from .core import NistLevels, NistLevelsClass

__all__ = ['NistLevels', 'NistLevelsClass',
           'Conf', 'conf',
           ]
