#!/usr/bin/env python

"""
This module is a placeholder for user defined stars for which the stellar spectrum
is created from the stellar parameters. Stars can either be created from a spectral
type or by adding a custom star to the list below.

PARAMETERS
----------
Teff : int, float
    Stellar effective temperature [astropy.units]
R : int, float
    Stellar radius [astropy.units]
M : int, float
    Stellar mass [astropy.units]
"""

import numpy as np
from astropy import units as u

#============================================================================#
#                             STAR FROM SPECTRAL TYPE                        #
#============================================================================#

# Function to create a star from its spectral type
# TODO not working yet!

def create_star(spectype, spectypes_datafile=None):


    if spectype in spectype_names: 
        spectype_names = list(np.genfromtxt(spectypes_datafile, dtype='str', delimiter='', usecols=(0), unpack=True))
        teffs, luminosities, radii,masses = np.genfromtxt(spectypes_datafile, delimiter='', usecols=(1,2,3,4), unpack=True)

    else:
        raise ValueError("SpType {} not one of {}".format(spectype,spectype_names))

    index = sptype_names.index(spectype)
    #-- then fill in the values
    star = {}
    star['teff'] = teffs[index],'K'
    star['radius'] = radii[index],'Rsol'
    star['mass'] = masses[index],'Msol'

    return star


# Based on tabulated data deduced from the spectral type

# if source == 'K5V':
#     M = create_star(source, spectypes_datafile='spectypes.dat')['mass']
#     R = create_star(source, spectypes_datafile='spectypes.dat')['radius']
#     g = (constants.GG_cgs * M[0] * constants.Msol_cgs / (R[0] * constants.Rsol_cgs)**2, 'cm s-2')
#     Teff = create_star(source, spectypes_datafile='spectypes.dat')['teff']
#     L = conversions.derive_luminosity(R, Teff, units='Lsol').as_tuple()  # "bolometric" luminosity


#============================================================================#
#                              USER DEFINED STARS                            #
#============================================================================#

# Based on manual settings for stellar parameters

def load_star(source):

    if source == 'GJ1214':  # M-dwarf
        M = 0.15  * u.M_sun
        R = 0.216 * u.R_sun
        Teff = 3026 * u.K
        logg = 4.5
        Z    = 0.0

    if source == 'WASP-43':  # K7 V
        # http://exoplanet.eu/catalog/wasp-33_b/
        M = 0.717 * u.M_sun
        R = 0.667 * u.R_sun
        Teff = 4520 * u.K
        logg = 4.5
        Z    = 0.0

    if source == 'CoRoT-1':  # G0 V
        # http://exoplanet.eu/catalog/wasp-33_b/
        M = 0.95 * u.M_sun
        R = 1.11 * u.R_sun
        Teff = 6298 * u.K
        logg = 4.5
        Z    = 0.0

    if source == "Sun":  # G0 V
        R = 1. * u.R_sun
        M = 1. * u.M_sun
        Teff = 5777. * u.K
        logg = 4.5
        Z    = 0.0

    if source == 'HD209458':
        M = 1.26 * u.M_sun
        R = 1.20 * u.R_sun
        Teff = 6071 * u.K
        logg = 4.5
        Z    = 0.0

    if source == 'WASP-33':  # A5 V
        # http://exoplanet.eu/catalog/wasp-33_b/
        M = 1.59 * u.M_sun
        R = 1.77 * u.R_sun
        Teff = 7430 * u.K
        logg = 4.5
        Z    = 0.0

    if source == 'Kepler-21':  # F6 IV
        # http://exoplanet.eu/catalog/wasp-33_b/
        M = 1.41 * u.M_sun
        R = 1.90 * u.R_sun
        Teff = 6305 * u.K
        logg = 4.5
        Z    = 0.0

    if source == 'dSct':
        M = 1.26 * u.M_sun
        R = 1.20 * u.R_sun
        Teff = 6071 * u.K
        logg = 4.0
        Z    = 0.0

    if source == 'gDor':
        M = 1.26 * u.M_sun
        R = 1.20 * u.R_sun
        Teff = 6071 * u.K
        logg = 4.0
        Z    = 0.0

        

    return M, R, Teff, logg, Z
