#!/usr/bin/env python

"""
In this module python file you can add your favorit exoplanet. NOTE in the
following the astropy-units are required, however, the the choice of units
(e.g. seconds vs. hours) are optional. The unit for each parameter is
reinforced to match the calculations within the script "simVariability.py".
knowing the following parameter space:

Parameters
----------
t0 : float
    Time of inferior conjunction (i.e. central time of first transit) [unit required]
P : float
   Orbital period [astropy.units]
a : float
    Semi-major axis [unit required]
i : float
    Orbital inclination (90-0) [astropy.units]
e : float
    Eccentricity (0-1)
w : float
    Longitude of periastron (0-360) [astropy.units]
rp : float
    Planet radius [astropy.units]
fp : float
    Planet-to-star flux ratio

Optional SPIDERMAN
------------------
xi : float
    Ratio of radiative to advective timescale
Tn : float
    Temperature of nightside [astropy.units]
dT : float
   Day-night temperature contrast [astropy.units]
"""

from astropy import units as u

def load_exoplanet(source):

    #-------------------------------------------------------#
    #                      SOLAR SYSTEM                     #
    #-------------------------------------------------------#

        
    if source == 'Jupiter':
        # Parameters are drawn from astropy
        params = {'t0': 10 * u.d,
                  'P' : 100 * u.d,
                  'i' : 90 * u.deg,
                  'e' : 0,
                  'w' : 90 * u.deg,
                  'rp': 0.5 * u.R_jup,
                  'mp': 0.5 * u.M_jup,
                  'xi': 0.0,
                  'Tn': 300 * u.K,
                  'dT': 300 * u.K}

    #-------------------------------------------------------#
    #                      HOT-JUPITERS                     #
    #-------------------------------------------------------#

    if source == 'hotJupiter':
        # Parameters are drawn from astropy
        params = {'t0': 1 * u.d,
                  'P' : 2 * u.d,
                  'i' : 90 * u.deg,
                  'e' : 0,
                  'w' : 90 * u.deg,
                  'rp': 1 * u.R_jup,
                  'mp': 1 * u.M_jup,
                  'xi': 0.0,
                  'Tn': 1128 * u.K,
                  'dT': 942 * u.K}

    if source == 'CoRoT-1b':
        # http://exoplanet.eu/catalog/corot-1_b/
        params = {'t0': 1 * u.d,
                  'P' : 1.5089557 * u.d,
                  'e' : 0.0,
                  'i' : 83.96 * u.deg,
                  'w' : 90.0 * u.deg,
                  'rp': 1.49 * u.R_jup,
                  'mp': 1.03 * u.M_jup,
                  'xi': 0.1,
                  'Tn': 1757 * u.K,
                  'dT': (3144 - 1757) * u.K}

    if source == 'WASP-33b':  # A5 V
        # http://exoplanet.eu/catalog/wasp-33_b/
        params = {'t0': 1 * u.d,
                  'P' : 1.21986967 * u.d,
                  'e' : 0.0,
                  'i' : 87.7 * u.deg,
                  'w' : 130.0 * u.deg,
                  'rp': 1.603 * u.R_jup,
                  'mp': 2.8 * u.M_jup,
                  'xi': 0.1,
                  'Tn': 1757 * u.K,
                  'dT': (3144 - 1757) * u.K}

    if source == 'WASP-43b':  # K7 V
        # http://exoplanet.eu/catalog/wasp-43_b/
        params = {'t0': 1 * u.d,
                  'P' : 0.81347753 * u.d,
                  'e' : 0.0035,
                  'i' : 82.33 * u.deg,
                  'w' : 328 * u.deg,
                  'rp': 1.036 * u.R_jup,
                  'mp': 2.052 * u.M_jup,
                  'xi': 0.0,
                  'Tn': 1000 * u.K,
                  'dT': 1000 * u.K}

    #-------------------------------------------------------#
    #                      EARTH ANALOGS                    #
    #-------------------------------------------------------#

    if source == 'hotMars':
        # 
        params = {'t0': 1 * u.d,
                  'P' : 2 * u.d,
                  'e' : 0.,
                  'i' : 88.0 * u.deg,
                  'w' : 0. * u.deg,
                  'rp': 0.531 * u.R_earth,
                  'mp': 0.107 * u.M_earth,
                  'xi': 0.,
                  'Tn': 300. * u.K,
                  'dT': 0. * u.K}


    if source == 'hotEarth':
        # 
        params = {'t0': 1 * u.d,
                  'P' : 50 * u.d,
                  'e' : 0.,
                  'i' : 90. * u.deg,
                  'w' : 0. * u.deg,
                  'rp': 1. * u.R_earth,
                  'mp': 1. * u.M_earth,
                  'xi': 0.,
                  'Tn': 300. * u.K,
                  'dT': 0. * u.K}

    if source == 'Earth':
        # 
        params = {'t0': 10 * u.d,
                  'P' : 365.25 * u.d,
                  'e' : 0.0167,
                  'i' : 90.0 * u.deg,
                  'w' : 0. * u.deg,
                  'rp': 1. * u.R_earth,
                  'mp': 1. * u.M_earth,
                  'xi': 0.,
                  'Tn': 300. * u.K,
                  'dT': 50. * u.K}

    if source == 'Neptune':
        # 
        params = {'t0': 10 * u.d,
                  'P' : 365.25 * u.d,
                  'e' : 0.0167,
                  'i' : 90.0 * u.deg,
                  'w' : 0. * u.deg,
                  'rp': 3.9 * u.R_earth,
                  'mp': 17.15 * u.M_earth,
                  'xi': 0.,
                  'Tn': 300. * u.K,
                  'dT': 50. * u.K}

    if source == 'Kepler-21b':  # F6 IV
        # http://exoplanet.eu/catalog/kepler-21_b/
        params = {'t0': 1 * u.d,
                  'P' : 2.78578 * u.d,
                  'e' : 0.02,
                  'i' : 83.96 * u.deg,
                  'w' : -15. * u.deg,
                  'rp': 1.636 * u.R_earth,
                  'mp': 5.079 * u.M_earth,
                  'xi': 0.,
                  'Tn': 300. * u.K,
                  'dT': 50. * u.K}



        

    return params
