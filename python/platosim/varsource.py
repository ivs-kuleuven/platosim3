#!/usr/bin/env python

"""
This is a script holding all relevant secondary synthetic model classes.
For an explaination of the full parameter space used in each class, have
a look at the "source_exoplanets.py" module.
"""

import os
import math
import random
import zipfile
import urllib.request

import pathlib
import numpy as np
from numba import njit
from astropy.io import fits
from astropy import units as u
from astropy import constants as c
from PyAstronomy import funcFit, pyasl
from matplotlib import pyplot as plt
from scipy.stats import norm, truncnorm

# PlatoSim
from platosim.utilities import errorcode, downloadFromFTP


#==============================================================#
#                       MODELS PARAMETERS                      #
#==============================================================#


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





def load_star(source):
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




def load_exoplanet(source):

    """Module containing a few standard exoplanet paramters. 

    NOTE in the following the astropy-units are required, however,
    the the choice of units (e.g. seconds vs. hours) are optional. 
    The unit for each parameter is reinforced to match the calculations
    within the script "simVariability.py". knowing the following 
    parameter space:

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





#==============================================================#
#                          MODELS CLASS                        #
#==============================================================#


class LimbDarkening(funcFit.OneDFit):
    """
    Class for fitting the Limb Darkening Coefficients.
    Integrated into the BATMAN and SPIDERMAN packages.
    """

    def __init__(self):
        funcFit.OneDFit.__init__(self, ['u1', 'u2'])

    def evaluate(self, x):

        model = 'quadratic'

        if model == 'linear':
            coef = ['u1']
        else:
            coef = ['u1', 'u2']

        y = 1. - self['u1']*(1. - x)  # Linear model
        if model == 'quadratic':   y = y - self['u2']*(1. - x)**2
        if model == 'square-root': y = y - self['u2']*(1. - np.sqrt(x))
        if model == 'logarithmic': y = y - self['u2']*x*np.log(x)
        if model == 'exponential': y = y - self['u2']/(1 - np.exp(x))
        if model == 'nonlinear':   y = 1. - self['u1']*(1 - np.sqrt(x)) - self['u2']*(1. - x)

        return y



    def limb_darkening_coefficients(self, bandpass, star_source, plot=False):
 
        """Calculate the Limb Darkening (LD) coefficient.

        To compute our custom limb-darkening transit duration coefficients that meet
        the PLATO transmission response function. We used the angle-dependent ("MU")
        Specific Intensity Spectra (SIS) from PHOENIX (Goettingen 2018), which exactly
        as above is a library of the stellar effective temperature, surface gravity,
        and metallicity. The limb darkening are naturally calclated for the exact same
        stellar parameter as used for the granulation and oscialltions. 

        See links:
        Webpage : https://phoenix.astro.physik.uni-goettingen.de/
        Download: http://phoenix.astro.physik.uni-goettingen.de/data/
        """

        #  Convert input parameters
        wvl_tele  = self.wvl_tele.to('Å').value
        tran_tele = self.tra_tele

        # Interpolate (piecewise cubic) into higher resolution grid
        grid_no  = 1000
        wvl_int  = np.linspace(wvl_tele[0], wvl_tele[-1], grid_no)
        passband = make_interp_spline(wvl_tele, tran_tele, k=3)
        tran_int = passband(wvl_int)

        # Limb darkening model options:
        if args.ldm: limbDarkModel = args.lmd
        else: limbDarkModel = 'quadratic'

        # ANGLE DEPENDENT SIS FROM PHOENIX
        # TODO use proper PHOENIX LD model

        # Get low resolution angle dependent SIS spectra:
        # - Each row in data is a spectrum
        # - mu is then the corresponding values for each spectrum
        phoenix = Phoenix()
        wvl_data, mu, data = phoenix.getSpecIntFITS(self.Teff.value, self.logg, self.Z)

        # Convert data from Vacuum-To-Air (VTA)
        data_int     = np.zeros(grid_no*len(mu)).reshape(grid_no, len(mu))
        data_int_VTA = np.zeros(grid_no*len(mu)).reshape(grid_no, len(mu))
        for jj in range(len(mu)):
            data_int[:,jj] = np.interp(wvl_int, wvl_data, data[jj])
            data_int_VTA[:,jj], vind = pyasl.specAirVacConvert(wvl_int, data_int[:,jj], direction="vactoair")

        # Colvolve data with transmission to find the VTA intensity
        intensity_VTA = [np.nansum(data_int_VTA[:,i]*tran_int) for i in range(len(mu))]
        intensity_VTA = intensity_VTA/max(intensity_VTA)

        # Interpolate to a finer grid
        mu_int = np.linspace(0, 1, grid_no)
        intensity_VTA_int = np.interp(mu_int, mu, intensity_VTA)

        # Get rid of the lowest mu's to do 2nd order fitting
        # NOTE The stellar limb causes a jump for mu<0.15
        ind = np.where(mu_int >= 0.15)
        mu_trunc = mu_int[ind]
        intensity_VTA_trunc = intensity_VTA_int[ind]

        # Initialize the class that will do the fitting, and prepare for the fit
        model_ldc = LimbDarkening()
        model_ldc.thaw(['u1', 'u2'])
        model_ldc.assignValue({'u1':0.2, 'u2':0.2})
        model_ldc.fit(mu_trunc, intensity_VTA_trunc, xtol=1.e-7, ftol=1.e-7, disp=0)
        ldc = [model_ldc["u1"], model_ldc["u2"]]

        # Plot results
        if args.plot:
            pt.plot_passband_ldc(wvl_int, tran_int, grid_no, mu_trunc,
                                 intensity_VTA_trunc, model_ldc, ldc)

        # Finito!
        return ldc


    


class DopplerBeaming(funcFit.OneDFit):
    """
    Class for the calculation of the Beaming Effect (also called
    Doppler Beaming/Boosting).

    This utility models the beaming effect descriped in Shporer (2017).
    The effect is composed by: 1) Doppler shift, 2) Time dialation
    3) Light abberation. The calculated beaming amplitude asssumes that
    v_RV << c.

    Resources
    ---------
    Sphorer (2019)         : https://arxiv.org/abs/1703.00496
    Murray & Correia (2011): https://arxiv.org/abs/1009.1738v2
    Loeb & Gaudi (2003)    : https://arxiv.org/abs/astro-ph/0303212

    Parameters
    ----------
    Keplerian model from PyAstronomy using parameter space from the code
    Batman, hence, see documentation for input parameters.

    Return
    ------
    Beaming model (light curve) and Amplitude in relative flux [ppm].
    """

    def __init__(self):
        funcFit.OneDFit.__init__(self, ['wvl_c', 'Teff',
                                        't0', 'P', 'a',
                                        'e', 'i', 'w',
                                        'Mp', 'Ms'])

    def evaluate(self, x, nu):

        # Convert to SI units

        x     = x.to('s')
        nu    = nu * u.rad
        e     = self['e']
        i     = self['i'].to('rad')
        w     = self['w'].to('rad')
        t0    = self['t0'].to('s')
        P     = self['P'].to('s')
        a     = self['a'].to('m')
        wvl_c = self['wvl_c'].to('m')
        Ms    = self['Ms'].to('kg')
        Mp    = self['Mp'].to('kg')
        Teff  = self['Teff'].to('K')

        # Correction factor (alpha) between true bolmetric flux and finite flux:
        # We use Sphorer (2017) Eq.5 analytical expression obtained approximating a blackbody spectrum.
        # In bolometric light alpha = 1, but otherwise deviating due to finite bandpass measurement.

        xx = c.h * c.c / (wvl_c * c.k_B * Teff)
        alpha = 1/4. * xx*np.exp(xx)/(np.exp(xx) - 1)

        # Amplitude: Sphorer (2019) Eq. 3

        A = (0.0028 * alpha * P.to('d')**(-1/3) * (Ms.to('M_sun') + Mp.to('M_sun'))**(-2/3) *
             Mp.to('M_sun') * np.sin(i)) * 1e6

        # RV signal as function of nu: Murray & Correia (2011) Eq. 61, 65 and 66
        # NOTE "astar" in the following is the reduced semimajor axies due to the common
        # center-of-mass and "K" is the relative RV semi-amplitude of the star.

        phi = (x.value - t0.value) / P.value * u.rad
        astar = Mp/(Mp + Ms) * a
        K     = 2.*np.pi/P * astar * np.sin(i)/np.sqrt(1. - e**2)
        #RV    = K * (np.cos(w + nu) - e*np.cos(w))
        RV    = K * (np.cos(2*np.pi*phi + w) - e*np.cos(w))

        # Final beaming effect [ppm]: Second term in Eq.1 from Sphorer (2017) but normalized

        y = 4 * alpha * RV/c.c * 1e6

        return y.value, A.value





class EllipsoidalDistortion(funcFit.OneDFit):
    """
    Class for the calculation of the ellipsoidal effect.

    The model presented here is a fusion between the simple analytical description
    presented by Sphorer (2019) and adding a parametization solving the elliptical
    model. Higher order terms from Morris & Nafilan (1993) may be desired to have
    a look at for future improvemenets.

    Resources:
    ----------
    Sphorer (2019)         : https://arxiv.org/abs/1703.00496
    Murray & Correia (2011): https://arxiv.org/abs/1009.1738v2
    Morris & Nafilan (1993): http://cdsads.u-strasbg.fr/pdf/1993ApJ...419..344M
    Claret & Bloemen (2011): TODO this is for TESS but we don't have g fror PLATO
    NOTE Only w=90 produces the expected beaming effect with max and min between
    inferior and superior conjunction.
    """

    def __init__(self):
        funcFit.OneDFit.__init__(self, ['g', 'u',
                                        't0', 'P',
                                        'e', 'i', 'w',
                                        'a', 'Rs',
                                        'Mp', 'Ms'])

    def evaluate(self, x, nu):

        # Convert to SI units

        x  = x.to('s')
        nu = nu * u.rad
        e  = self['e']
        i  = self['i'].to('rad')
        w  = self['w'].to('rad')
        t0 = self['t0'].to('s')
        P  = self['P'].to('s')
        a  = self['a'].to('m')
        Rs = self['Rs'].to('m')
        Ms = self['Ms'].to('kg')
        Mp = self['Mp'].to('kg')

        # Orbital phase

        phi = (x.value - t0.value) / P.value

        # Coefficient accounting for the stellar LD and GD: Sphorer (2019) Eq.8

        alpha = 0.15 * (15. + self["u"]) * (1. + self["g"])/(3. - self["u"])

        # Amplitude of the ellipsoidal variation: Sphorer (2019) Eq.7

        A = alpha * Mp/Ms * (Rs/a)**3 * np.sin(i)**2 * 1e6

        # Add parameterization of elliptical orbits: Murray & Correia (2011) Eq.20
        # This is just 1 for e=0

        a1 = (1 + e*np.cos(nu)) / (1. - e**2)

        # Add parameterization of the ellipsoidal distortion itself

        a2 = - np.cos(4 * np.pi * phi)

        # Final ellipsoidal distortion [ppm]

        y = A * a1 * a2

        return y.value, A.value







class StellarFlares(funcFit.OneDFit):  # TODO test
    """
    A simplistic analytical description of stellar flares described by an sudden
    flux increase followed by an exponential decay. Given the time of the time
    series the corresponding flux is returned including the wanted flares.

    PARAMETERS
    ----------
    t_flares : float, narray
        Times for the maxima of each flare occurances [days]
    t_scale : float, narray
        The full time width at half-maximum-flux of each flare [days]
    asymmetry : float, narray
        Higher order asymmetry factor for each flare.

    RETURN
    ------
    flux : narray
    """
    def __init__(self):
        funcFit.OneDFit.__init__(self, ["sampling", "t_flares", "t_scale", "asymmetry"])

    def evaluate(self, time):

        t_step   = np.diff(time[:2])[0]
        n_flares = len(self["t_flares"])
        flux     = np.zeros_like(np.arange(time[0], time[-1], t_step))

        for m in range(n_flares):

            t0 = time[0]  - self["t_flares"][m]
            t1 = time[-1] - self["t_flares"][m]
            tn = np.arange(t0, t1, t_step)
            t  = tn / self["t_scale"]

            B = self["asymmetry"]
            C = 1./B
            b = - 1.941 - 0.175 + 2.246 + 1
            c = 1 - 0.689

            # Loop over every time-step in the time interval that is defined relative to this flares maxima
            # and put in units of the time-scale

            for i in range(len(t)):
                if (t[i]*B) > -1 and (t[i]*B) <= 0:
                    flux[i] += 1 + 1.941 * t[i]*B - 0.175 * np.power(t[i]*B,2) - 2.246 * np.power(t[i]*B,3) - b * np.power(t[i]*B,4)

                elif t[i]*C > 0:
                    flux[i] += 0.689 * np.exp(-1.6 * t[i]*C) + c * np.exp(-0.2783 * t[i]*C)

                else:
                    flux[i] += 0


        time = tn + self["t_flares"][n_flares-1]

        return flux






class GravityOscillations(funcFit.OneDFit):  # TODO test
    """
    time_start and time_end are given in days and define the time interval over which to simulate the flare
    the sampling should be given in exposures or data-points per day
    the period and amplitude define the allowed range in periods and amplitudes to include
    the number_modes defines how many different modes or periods/amplitudes to include
    power is the power to which the sum of every mode is raised, it introduces an asymmetry in the signal

    """
    def __init__(self):
        funcFit.OneDFit.__init__(self, ["P_range", "A_range", "N_modes", "power", "t_step"])

    def evaluate(self, time):

        seed = random.seed(seed)

        # convert the times in days to times in secons (because the time_step is given in seconds)
        
        self["t_step"] /= 24*60*60.

        t = np.arange(time[0], time[-1], self["t_step"], dtype=float)
        flux = np.zeros_like(t)

        # loop over the number of modes and each time pick a period, an amplitude and a phi
        #  at random out of the appropriate range
        
        for i in range(self["N_modes"]):

            phi = random.uniform(0, 2 * np.pi)
            P   = random.uniform(self["P_range"][0], self["P_range"][1])
            A   = random.uniform(self["A_range"][0], self["A_range"][1])

            flux += A * np.sin(2 * np.pi * (1 / P) * t + phi)  # sum of every mode

        # normalize the flux so its values lie in [-1, 1] (so roots are not undefined)
        
        A     = np.amax(np.absolute(flux))
        flux /= A
        flux  = A * (((flux + 1)**self["power"]) - 1)

        # That's it!
        
        return flux








# class SolarLikeOscillations(object):

#     def __init__(self):
#         """
#         PURPOSE: Initialize and prepare data structure
#         """
#         # Create data directories if they do not exist

@njit
def pulsations(time, freq, eta, Ntime, Nmode, amplsin, amplcos,
               kicktimestep, kick_amplitude, last_kicktime, next_kicktime):

    # Prepare for loop
    signal = np.zeros(Ntime)
    for j in range(Ntime):

        # Compute the contribution of each mode separatly
        for i in range(Nmode):

            # Let the oscillator evolve until right before 'time[j]'
            while (next_kicktime[i] <= time[j]):

                deltatime = next_kicktime[i] - last_kicktime[i]
                damp = np.exp(-eta[i] * deltatime)
                amplsin[i] = damp * amplsin[i] + kick_amplitude[i] * np.random.normal(0.,1.)
                amplcos[i] = damp * amplcos[i] + kick_amplitude[i] * np.random.normal(0.,1.)
                last_kicktime[i] = next_kicktime[i]
                next_kicktime[i] = next_kicktime[i] + kicktimestep

            # Now make the last small step until 'time[j]'
            deltatime = time[j] - last_kicktime[i]
            damp = np.exp(-eta[i] * deltatime)
            signal[j] = signal[j] + damp * (amplsin[i] * np.sin(2*np.pi*freq[i]*time[j])    \
                                  + amplcos[i] * np.cos(2*np.pi*freq[i]*time[j]))

    return(signal)







def solarosc(time, freq, ampl, eta, verbose):
    """
    Compute time series of stochastically excited damped modes

    See also De Ridder et al., 2006, MNRAS 365, pp. 595-605.

    Example:

    >>> time = np.linspace(0, 40, 100)      # in Ms
    >>> freq = np.array([23.0, 23.5])       # in microHz
    >>> ampl = np.array([100.0, 110.0])     # in ppm
    >>> eta = np.array([1.e-6, 3.e-6])      # in 1/Ms
    >>> oscsignal = solarosc(time, freq, ampl, eta)
    >>> flux = 1000000.0                      # average flux level
    >>> signal = flux * (1.0 + oscsignal)
    >>> # The same with a logger
    >>> import sys, logging, logging.handlers
    >>> myLogger = logging.getLogger("solarosc")
    >>> myLogger.addHandler(logging.StreamHandler(sys.stdout))
    >>> myLogger.setLevel(logging.INFO)
    >>> oscsignal = solarosc(time, freq, ampl, eta, myLogger)
    Simulating 2 modes
    Oscillation kicktimestep: 3333.333333
    300 kicks for warm up for oscillation signal
    Simulating stochastic oscillations

    @param time: time points [0..Ntime-1] (unit: e.g. Ms)
    @type time: ndarray
    @param freq: oscillation freqs [0..Nmodes-1] (unit: e.g. microHz)
    @type freq: ndarray
    @param ampl: amplitude of each oscillation mode
                 rms amplitude = ampl / sqrt(2.)
    @type ampl: ndarray
    @param eta: damping rates (unit: e.g. (Ms)^{-1})
    @type eta: ndarray
    @return: signal[0..Ntime-1]
    @rtype: ndarray
    """

    Ntime = len(time)
    Nmode = len(freq)

    if verbose:
        print("Simulating %d modes" % Nmode)

    # Set the kick (= reexcitation) timestep to be one 100th of the
    # shortest damping time. (i.e. kick often enough).
    kicktimestep = (1.0 / max(eta)) / 100.0

    if verbose:
        print("Oscillation kicktimestep: %f" % kicktimestep)
        
    # Init start values of amplitudes, and the kicking amplitude
    # so that the amplitude of the oscillator will be on average be
    # constant and equal to the user given amplitude

    amplcos = 0.0
    amplsin = 0.0
    kick_amplitude = ampl * np.sqrt(kicktimestep * eta)

    # Warm up the stochastic excitation simulator to forget the
    # initial conditions. Do this during the longest damping time.
    # But put a maximum on the number of kicks, as there might
    # be almost-stable modes with damping time = infinity

    damp = np.exp(-eta * kicktimestep)
    Nwarmup = min(20000, int(math.floor(1.0 / np.min(eta) / kicktimestep)))

    if verbose:
        print("Simulating stochastic oscillations")
        print("%d kicks for warm up for oscillation signal" % Nwarmup)

    for i in range(Nwarmup):
        amplsin = damp * amplsin + np.random.normal(np.zeros(Nmode), kick_amplitude)
        amplcos = damp * amplcos + np.random.normal(np.zeros(Nmode), kick_amplitude)

    # Initialize the last kick times for each mode to be randomly chosen
    # a little before the first user time point. This is to avoid that
    # the kicking time is always exactly the same for all of the modes.

    last_kicktime = np.random.uniform(time[0] - kicktimestep, time[0], Nmode)
    next_kicktime = last_kicktime + kicktimestep

    signal = pulsations(time, freq, eta, Ntime, Nmode, amplsin, amplcos,
                        kicktimestep, kick_amplitude, last_kicktime, next_kicktime)

    # Finito!

    return signal










def add_region(nday, ic, lon, lat, k, bsiz1, phase):
# Add one active region of a particular size:
    w_org = 0.4 * bsiz1                            # original width (degrees)
    width = 4.0                                    # final width (degrees)
    bmax = 250. * (w_org / width)**2               # final peak flux density (G) 
    bsizr = np.pi * bsiz1 / 180                    # pole separation in radians
    width = np.pi * width / 180                    # final width in radians
    # For tilt angles, see Wang and Sheeley, Sol. Phys. 124, 81 (1989)
    #                      Wang and Sheeley, Ap. J. 375, 761 (1991)
    #                      Howard, Sol. Phys. 137, 205 (1992)
    x = np.random.normal()
    while abs(x) > 1.6:
        x = np.random.normal()
    y = np.random.normal()
    while abs(y) >= 1.6:
        y = np.random.normal()
    z = np.random.uniform()
    if z > 0.14:
        ang = 0.5 * lat + 2.0 + 27. * x * y # tilt angle (degrees)
    else:
        while z > 0.5:
            z = np.random.normal()
        ang =  z * np.pi / 180
    lat = np.pi * lat / 180                               # latitude (radians)
    ang = np.pi * ang / 180                               # tilt angle (radians)
    dph = ic * 0.5 * bsizr * np.cos(ang) / np.cos(lat)    # delta phi (radians)
    dth = ic * 0.5 * bsizr * np.sin(ang)                  # delta theta (radians)
    phcen = np.pi * lon / 180                             # longitude (radians)
    if k == 0:                       # Insert on N hemisphere
        thcen = 0.5 * np.pi - lat
        phpos = phcen + dph
        phneg = phcen - dph
    else:                            # Insert on S hemisphere
        thcen = 0.5 * np.pi + lat
        phpos = phcen - dph
        phneg = phcen + dph

    thpos = thcen + dth
    thneg = thcen - dth
    str_ = '{:5d} {:8.5f} {:8.5f} {:8.5f} {:8.5f} {:8.5f} {:7.1f} {:8.5f}'.format(int(nday),thpos,phpos,thneg,phneg,width,250,ang)
    return str_          

def regions(seed = None, activityrate = 1, cycle_period = 1, cycle_overlap = 0, maxlat = 70,
            minlat = 0, tsim = 1000, tstart = 0, randspots = False, odir = None, verbose  = True):
# ; inputs
# ; activityrate - number of bipoles (1= solar)
# ; cyclelength - length of cycle in years
# ; cycleoverlap - cycleoverlap time in years
# ; tsim - length of simulation in days 
# ; tstart - first day to start outputting bipoles 
# ; minlat - minimum latitude of spot emergence
# ; maxlat - maximum latitude of spot emergence 
# ; randspots - if True, no Butterfly pattern
# ; odir - output directory 
# ; verbose - if True, output printed to stdout

# ; output: list of regions with parameters:
# ;
# ;    nday = day of emergence
# ;    thpos= theta of positive pole (radians)
# ;    phpos= phi   of positive pole (radians)
# ;    thneg= theta of negative pole (radians)
# ;    phneg= phi   of negative pole (radians)
# ;    width= width of each pole (radians)
# ;    bmax = maximum flux density (Gauss)
# ;
# ;  According to Schrijver and Harvey (1994), the number of active regions
# ;  emerging with areas in the range [A,A+dA] in a time dt is given by 
# ;
# ;    n(A,t) dA dt = a(t) A^(-2) dA dt ,
# ;
# ;  where A is the "initial" area of a bipole in square degrees, and t is
# ;  the time in days; a(t) varies from 1.23 at cycle minimum to 10 at cycle
# ;  maximum.
# ;
# ;  The bipole area is the area within the 25-Gauss contour in the
# ;  "initial" state, i.e. time of maximum development of the active region.
# ;  The assumed peak flux density in the initial sate is 1100 G, and
# ;  width = 0.2*bsiz (see disp_region). The parameters written onto the
# ;  file are corrected for further diffusion and correspond to the time
# ;  when width = 4 deg, the smallest width that can be resolved with lmax=63.
# ;
# ;  In our simulation we use a lower value of a(t) to account for "correlated"
# ;  regions.
# ;

    nbin=5                              # number of area bins
    delt=0.5                            # delta ln(A)
    amax=100.                           # orig. area of largest bipoles (deg^2)
    dcon = np.exp(0.5*delt)-np.exp(-0.5*delt)   # contant from integ. over bin

    if verbose:
        print('Creating regions with the following parameters:')
        print('Acivity rate: {} x Solar rate.'.format(activityrate))
        print('Activity cycle period: {} years.'.format(cycle_period))
        print('Duration of overlap between consecutive activity cycles: {} years.'.format(cycle_overlap))
        print('Maximum spot latitude: {} degrees.'.format(maxlat))
        print('Minimum spot latitude: {} degrees.'.format(minlat))
        print('Duration of simulation: {} days.'.format(tsim))
        print('Time at start of simulation: {} days.'.format(tstart))
        
    latrmsd = 5
    atm = 10 * activityrate    
    # a(t) at cycle maximum (deg^2/day)
    # cycle period (days)
    # cycle duration (days)
     
    ncycle = int(cycle_period * 365)         # cycle length in days   
    nclen = int((cycle_period + cycle_overlap) * 365)
    fact = np.exp(delt*np.arange(nbin))     # array of area reduction factors
    ftot = fact.sum()                         # sum of reduction factors
    bsiz = np.sqrt(amax/fact)               # array of bipole separations (deg)
    tau1 = 5                                  # first and last times (in days) for
    tau2 = 15                                 #   emergence of "correlated" regions
    prob = 0.001                              # total probability for "correlation"
    nlon = 36                                 # number of longitude bins
    nlat = 16                                 # number of latitude bins       
    nday1 = 0                                 # first day to be simulated
    ndays = int(tsim)                              # number of days to be simulated
    dt = 1

    # Initialize random number generator
    np.random.seed(seed)
                
    # Initialize time since last emergence of a large region, as function
    # of longitude, latitude and hemisphere:
    tau = np.zeros((nlon,nlat,2),'int') + tau2
    dlon = 360. / nlon
    dlat = maxlat / nlat
                  
    # Open file for results
    if odir is None: 
        odir = os.getcwd() + '/output'
    ofile = os.path.join(odir, 'regions.txt')
        
    with open(ofile, 'w') as flo:
        ncnt = 0
        # Loop over time (in days):
        ncur = 0
        start_day = 0
        
        for nd in range(ndays):
#            print('nd:', nd)
            nday = nd + nday1
            
            # Compute index of most recently started cycle:
            ncur_now = int(nday / ncycle)
            ncur_prev = int((nday-1) / ncycle)
            if ncur_now > ncur_prev:
                ncur = ncur + 1
#            print('ncur:', ncur)
            #  Initialize rate of emergence for largest regions, and add 1 day
            #  to time of last emergence:

            tau = tau + 1
            rc0 = np.zeros((nlon,nlat,2))
            l = (tau > tau1) & (tau <= tau2)
            if l.any():
                rc0[l] = prob / (tau2 - tau1)
 
            #  Loop over current and previous cycle:
            for icycle in [0,1]:
                nc = ncur-icycle # index of cycle
#                print('icycle:', icycle)
#                print('nc:', nc)
                if ncur == 0:
                    start_day = nc * ncycle
                else:  
                    if ncur == 1:
                        if icycle == 0:
                            start_day = ncycle * nc
                        elif icycle == 1:
                            start_day = 0
                    else:
                        start_day = ncycle * nc
#                print('start_day', start_day)
           
                nstart = start_day        # start date of cycle
                if (nday-nstart) < nclen:  
#                    print('inif')
                    ic = 1 - 2 * ((nc + 2) % 2) # +1 for even, -1 for odd cycle
                    phase = float(nday-nstart) / nclen # phase within the cycle
#                    print('Cycle phase:', phase)
                    # Emergence rate of largest "uncorrelated" regions (number per day,
                    # both hemispheres), from Schrijver and Harvey (1994):
            
                    ru0_tot = atm * np.sin(np.pi * phase)**2 * (1.0 * dcon) / amax
            
                    # Emergence rate of largest "uncorrelated" regions per latitude/longitude
                    # bin (number per day), as function of latitude:
                    if randspots:
                        latavg = (maxlat - minlat) / 2. 
                        latrms = maxlat - minlat
                        nlat1 = np.floor(minlat / dlat).astype(int)
                        nlat2 = np.floor(maxlat / dlat).astype(int)
                        nlat2 = min([nlat2, nlat - 1])
                    else:
                        latavg = maxlat + (minlat - maxlat)*phase #+ 5.*phase**2
                        # latavg=70.0-68.*phase+5.*phase**2 # average latitude (degrees)
                        latrms = (maxlat/5.) - latrmsd * phase # rms latitude (degrees)
                        nlat1 = np.floor(max([maxlat * 0.9 - 1.2 * maxlat * phase, 0.0]) / dlat).astype(int) # first and last index
                        nlat2 = np.floor(min([maxlat + 15. - maxlat * phase, maxlat]) / dlat).astype(int)
                        nlat2 = min([nlat2, nlat - 1])

                    js = np.arange(nlat2 - nlat1).astype(int)

                    p = np.zeros(nlat)
                    for j in np.arange(nlat2-nlat1+1).astype(int) + nlat1:
                        p[j] = np.exp( - ((dlat * (0.5 + j) - latavg) / latrms)**2)
                    ru0 = ru0_tot * p / (p.sum() * nlon * 2)
            
                    # Loops over hemisphere and latitude:
                    for k in [0,1]:
                        for j in np.arange(nlat2-nlat1+1).astype(int) + nlat1:
                            # Emergence rates of largest regions per longitude/latitude bin (number
                            # per day):
                            r0 = ru0[j] + rc0[:,j,k]
                            rtot = r0.sum()
                            ssum = rtot * ftot
                            x = np.random.uniform()
                            if x <= ssum:
                                nb = 0
                                sumb = rtot * fact[0]
                                while x > sumb:
                                    nb = nb + 1
                                    sumb = sumb + rtot * fact[nb]
                                i = 0
                                sumb = sumb + (r0[0] - rtot) * fact[nb]
                                while x > sumb:
                                    i = i + 1
                                    sumb = sumb + r0[i] * fact[nb]
                                lon = dlon * (np.random.uniform() + float(i))
                                lat = dlat * (np.random.uniform() + float(j))
                                if (nday > tstart):
                                    str_ = add_region(nday/dt,ic,lon,lat,k,bsiz[nb],phase)
                                    if verbose:
                                        print(str_)
                                    flo.write('{}\n'.format(str_))
                                ncnt = ncnt + 1
                                if nb < 1:
                                    tau[i,j,k] = 0
    flo.close()
    if verbose:
        print('Total number of regions:  ',ncnt)
    return

def butterfly(odir=None):
    if not odir:
        odir = os.getcwd()
    regions = np.genfromtxt(os.path.join(odir, 'regions.txt')).T
    print(regions[7])
    nreg = len(regions)
    angle_regions = 0.5 * (regions[1] + regions[3])
    lats_regions = np.pi/2 - angle_regions
    l = angle_regions < 0.0
    lats_regions[l] = angle_regions[l] - np.pi/2
    lats_regions = lats_regions * 180 / np.pi
    plt.figure(figsize=(12,4))
    plt.plot(regions[0], lats_regions, 'ko')
    plt.ylim(-90, 90)
    plt.ylabel('Latitude (deg)')
    plt.xlabel('Time (days)')
    plt.savefig(os.path.join(odir,'butterfly.png'))

if __name__ == "__main__":
    regions(cycle_overlap=0.2,verbose=False)
    butterfly()
    plt.show()











class StellarActivity(object):
    """
    A simplistic sinusoidal model of stellar activity.
    """
    def __init__(self):
        
        self.Prot_sun  = 27.0
    

    def run_model(self, time, timedur, cadence, nsim=1,
                  odir=None, seed=None, verbose=True, plot=False):
        """
        PURPOSE: Quicktool to to create a noise-less time series.
        """
        # Run script to produce the spot cycles
        self.make_starspot_regions(tsim=timedur)
        # Run script to produce light curves
        self.make_starspot_lightcurve(time=time, dur=timedur, cadence_hours=cadence,
                                      verbose=verbose, plot=plot)
        # Plot if defined by user
        #self.make_plot()
        # Return data
        return self.dF.sum(0)

    
    def make_plot(self):
        fig, axes = plt.subplots(3,1, figsize=(12,9), sharex=True)
        X = np.genfromtxt(self.reg_file).T
        i = 0
        #string = \
        # """
        # Stellar activity cycles
        # AR = {:5.3f}, CL = {:6.3f}, R = {:1d} 
        # sini = {:6.2f}, PEQ = {:6.2f}, PPL = {:6.2f} 
        # tau  = {:5.2f}, NS = {:5d}
        # """.format(self.ar[i],
        #            self.clen[i],
        #            int(self.rr[i]),
        #            np.sin(self.incl[i]),
        #            self.period_eq[i],
        #            self.period_pole[i],
        #            self.tau_evol[i],
        #            self.s.nspot)
        axes[0].set_title('Stellar Acticity Cycle')
        #axes[1].text(self.time[-1]-self.time[-1]/10, 0, string)
        for j in range(self.s.nspot):
            
            axes[0].plot(self.s.t0[j], self.s.lat[j]*180/np.pi, 'ko',
                         ms=self.s.amax[j]*(1./3e-4)*5, alpha=0.5)
            axes[0].set_ylim(-90,90)
            axes[0].set_ylabel('Spot latitude [deg]')

            axes[1].plot(self.time, self.area.sum(0)*100, 'k-')
            axes[1].set_ylabel('Spot coverage [%]')
            
            axes[2].plot(self.time, self.dF.sum(0), 'k-')
            axes[2].set_ylabel('Relative flux')
            axes[2].set_xlim(self.time.min(), self.time.max())
            axes[2].set_xlabel('Time [days]')
            
            #plt.savefig(os.path.join(odir, 'lightcurve_{:04d}.png'.format(i)))
        #plt.tight_layout()
        plt.show()
        plt.close('all')




        
        
    def make_starspot_regions(self, tsim, odir=None, nsim=1, seed=None, verbose=True):

        # Handle output directory
        if odir is None:
            odir = os.getcwd() + '/output'
        ofile = os.path.join(odir, 'regions_par.txt')

        # Verbosity
        header = '{:4s} {:5s} {:6s} {:6s} {:6s} {:6s} {:4s}'.format('NO', 'AR', 'CLEN', 'COVER', 'LMIN', 'LMAX', 'RAND')
        if verbose:
            print(header)
        flo = open(os.path.join(odir, 'regions_par.txt'), 'w')
        flo.write('{}\n'.format(header))

        # Draw global parameters from random distributions
        np.random.seed(seed)

        # Activity cycle period, uniform between 1 and 10 years
        clen = np.random.uniform(1, 10, nsim) 

        # Overlap between consecutive cycles, uniform between 0 and 0.1 cycle period
        coverlap = np.random.uniform(0, 0.1, nsim) * clen 

        # Rate of emergence of magnetic bipoles:
        # log uniform between 0.3 (10^-0.5) lower and 3x (10^0.5) times solar
        ar = 10.0**(np.random.uniform(-0.5, 0.5, nsim))

        # minimum spot latitute, uniform between 0 and 40 degrees
        minlat = np.random.uniform(0, 40, nsim)

        # maximum spot latitute, ranging from minlat to 80 degrees
        maxlat = np.random.uniform(0,1,nsim)**0.3 * (80 - minlat) + minlat

        # A fifth of the LCs will have no butterfly pattern, but the rate of emergence
        # of active regions still fluctuates periodically with period clen
        randspots = np.zeros(nsim, 'bool')
        randspots[0:int(nsim/5)] = True
        np.random.shuffle(randspots)


        for i in range(nsim):
            lin = '{:4d} {:5.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:1d}'.format(i, ar[i], clen[i], coverlap[i], minlat[i], maxlat[i], randspots[i])
            if verbose:
                print(lin)
            flo.write('{}\n'.format(lin))

            # Simulate regions
            regions(activityrate = ar[i], cycle_period = clen[i], cycle_overlap = coverlap[i],
                    tsim = max(tsim, 1.2 * clen[i] * 365.25), minlat = minlat[i], maxlat = maxlat[i],
                    randspots = randspots[i], verbose=False) 

            # Extract random subset 
            rfile_tmp = os.path.join(odir, 'regions.txt')
            rfile_fin = os.path.join(odir, 'regions_{:04d}.txt'.format(i))

            X = np.genfromtxt(rfile_tmp).T
            t = X[0]
            tlen = int(t.max())
            if tlen > tsim:
                tstart = np.random.randint(tlen - tsim)
                l = (t >= tstart) & (t < (tstart + tsim))
                X = X[:,l]
                X[0] -= tstart
            np.savetxt(rfile_fin,X.T)

        flo.close()


    def make_starspot_lightcurve(self, time, seed=None, dur=1100, cadence_hours=2.0, 
                                 verbose=True, odir=None, plot=True):

        cad = cadence_hours / 24.0 # in days
        if odir is None:
            odir = os.getcwd() + '/output'

        # read in global activity cycle parameters
        reg_par_file = os.path.join(odir, 'regions_par.txt')
        X = np.genfromtxt(reg_par_file, skip_header=1).T
        no, ar, clen, cover, lmin, lmax, rr = X

        # Avoid a crash when only one star is requested
        try:
            nsim = len(ar)
        except TypeError:
            nsim = 1
            ar = [ar]
            clen = [clen]
            cover = [cover]
            lmin = [lmin]
            lmax = [lmax]
            rr = [rr]

        # draw stellar parameters from random distributions
        np.random.seed(seed)
        # random orientations -> uniform in sin(i)
        incl = np.arcsin(np.random.uniform(0, 1, nsim))
        # (equatorial) rotation period uniform between 3 and 90 days
        period_eq = np.random.uniform(3, 90.0, nsim)
        omega_eq = self.Prot_sun / period_eq # equatorial rotation rate (in radians)
        # relative differential rotation rate uniform between 0 and 0.25 (solar is about 0.15)
        delta_omega_rel = np.random.uniform(0, 0.25, nsim)
        delta_omega = delta_omega_rel * omega_eq # absolute differential rotation in radians
        period_pole = period_eq / (1 - delta_omega_rel) # polar period
        # spot half life ranges of 1 to 10 rotation periods (log uniform)
        tau_evol = 10.0**(np.random.uniform(0, 1, nsim))
        # average size of spots (using this scaling approximately matches the amplitude of variability for the Sun if ar=1)
        alpha_med = np.sqrt(ar) * 3e-4

        spot_par_file = os.path.join(odir, 'spot_par.txt')
        flo = open(spot_par_file, 'w')
        str_ = '#  N    AR   CLEN  COVER   LMIN   LMAX R   SINI    PEQ   PPOL  A_MED  TAU  NSPOT'
        flo.write('{}\n'.format(str_))

        if verbose: print(str_)
            
        for i in range(nsim):
            reg_file = os.path.join(odir, 'regions_{:04d}.txt'.format(i))

            
            s = StarSpots(incl = incl[i], omega = omega_eq[i],
                          delta_omega = delta_omega_rel[i], alpha_med = alpha_med[i],
                          tau_evol = tau_evol[i], regions = reg_file)

            str_ = '{:4d} {:5.3f} {:6.3f} {:6.3f} {:6.3f} {:6.3f} {:1d} {:6.2f} {:6.2f} {:6.2} {:5.2f} {:5.2f} {}'.format(i, ar[i], clen[i], cover[i], lmin[i], lmax[i], int(rr[i]), np.sin(incl[i]), period_eq[i], period_pole[i], np.log10(alpha_med[i]), tau_evol[i], s.nspot)
            flo.write('{}\n'.format(str_))

            if verbose:print(str_)

            tstart = s.dur - dur
            #time   = np.r_[tstart:s.dur:cad]
            area, ome, beta, dF = s.calc(time)
        
            X = np.zeros((2,len(time)))
            X[0,:] = time
            X[1,:] = dF.sum(0)
            #np.savetxt(os.path.join(odir, 'lightcurve_{:04d}.txt'.format(i)), X.T)
            # Create plot if defined by user
            self.reg_file = reg_file
            self.ar = ar
            self.clen = clen
            self.rr = rr
            self.incl = incl
            self.period_eq = period_eq
            self.period_pole = period_pole
            self.tau_evol = tau_evol
            self.s = s
            self.dF = dF
            self.time = time
            self.area = area
            self.regions = regions
            if plot:
                self.make_plot()

        flo.close()
        return






    


        
    
class StarSpots():
    """Holds parameters for spots on a given star"""
    def __init__(self, dur = None, alpha_med = 0.0001, incl = None,
                 omega = 2.0, delta_omega = 0.3, 
                 tau_evol = 5.0, threshold = 0.1, 
                 regions = '/Users/aigrain/Soft/idl/diffrot/regions.txt'):
        '''Generate initial parameter set for spots (emergence times
        and initial locations are read from regions file)'''
        # set global stellar parameters which are the same for all spots
        # inclination
        if incl == None:
            self.incl = np.arcsin(np.random.uniform())
        else:
            self.incl = incl
        # rotation and differential rotation (supplied in solar units)
        self.omega_sun = 2 * np.pi / (27.0 * 86400)
        self.omega = omega * self.omega_sun                                  # [rad]
        self.delta_omega = delta_omega * self.omega_sun 
        self.per_eq = 2 * np.pi / self.omega / 86400                         # [day]
        self.per_pole = 2 * np.pi / (self.omega - self.delta_omega) / 86400  # [day]
        #self.diffrot_func = diffrot_func
        # spot emergence and decay timescales
        self.tau_em = min(2.0, self.per_eq * tau_evol / 10.0)
        self.tau_decay = self.per_eq * tau_evol
        # read in regions file
        X = np.genfromtxt(regions).T
        t0 = X[0].astype(float)
        lat = 0.5*(X[1]+X[3])
        l = lat < 0
        lat = np.pi/2. - lat
        lat[l] = -lat[l]
        lon = 0.5*(X[2]+X[4])
        Bem = X[7]
        # keep only spots emerging within specified time-span, with peak B-field > threshold
        if dur == None:
            self.dur = t0.max() 
        else:
            self.dur = dur
        l = (t0 < self.dur) * (Bem > threshold)
        self.nspot = l.sum()
        self.t0 = t0[l]
        self.lat = lat[l]
        self.lon = lon[l]
        # scale to achieve desired median alpha, where alpha = spot contrast * spot area
        self.amax = Bem[l] * alpha_med / np.median(Bem[l])
                                               

    def diffrot_func(self, name, omega_0, delta_omega, lat):
        """
        # at pole lat = 90 deg then sin(lat) = 1 and omega_90 = omega_0 - delta_omega, 
        # therefore period_90 = 2pi/omega_90 = 2 pi / (omega_0 - delta_omega)
        # but omega_0 = 2 pi / period_0 and delta_omega = delta_omega_rel * omega_0 so
        # period_pole = 2 pi / (2 pi / period_0 - delta_omega_rel * 2 pi / period_0)
        #             = 2 pi / (2 pi / period_0 (1 - delta_omega_rel))
        #             = period_0 / (1 - delta_omega_rel)
        """
        if name == 'sin2': 
            return omega_0 - delta_omega * np.sin(lat)**2
                                               
    def calci(self, time, i):
        '''Evolve one spot and calculate its impact on the stellar flux'''
        '''NB: Currently there is no spot drift or shear'''
        # Spot area
        area = np.ones(len(time)) * self.amax[i]
        tt = time - self.t0[i]
        l = tt<0
        area[l] *= np.exp(-tt[l]**2 / 2. / self.tau_em**2) # emergence
        l = tt>0
        area[l] *= np.exp(-tt[l]**2 / 2. / self.tau_decay**2) # decay
        # Rotation rate
        ome = self.diffrot_func('sin2', self.omega, self.delta_omega, self.lat[i])
        # Fore-shortening 
        phase = ome * time * 86400 + self.lon[i]
        beta = np.cos(self.incl) * np.sin(self.lat[i]) + \
            np.sin(self.incl) * np.cos(self.lat[i]) * np.cos(phase)
        # Differential effect on stellar flux
        dF = - area * beta
        dF[beta < 0] = 0
        return area, ome, beta, dF

    def calc(self, time):
        '''Calculations for all spots'''
        N = len(time)
        M = self.nspot
        area = np.zeros((M, N))
        ome = np.zeros(M)
        beta = np.zeros((M, N))
        dF = np.zeros((M, N))
        for i in np.arange(M):
            area_i, omega_i, beta_i, dF_i = self.calci(time, i)
            area[i,:] = area_i
            ome[i] = omega_i
            beta[i,:] = beta_i
            dF[i,:] = dF_i
        return area, ome, beta, dF




    
class PlanetMRforecast():
    """
    Class to forecast the mass from a planets radius.
    """

    def __init__(self):
        
        # constant
        mearth2mjup = 317.828
        mearth2msun = 333060.4
        rearth2rjup = 11.21
        rearth2rsun = 109.2

        # Boundary
        mlower = 3e-4
        mupper = 3e5

        # Number of different populations
        n_pop = 4

        # read parameter file
        filepath = 'inputfiles/data_varsim/varsim_exomass_fitting_parameters.h5' 
        hyper_file = Path(os.getenv("PLATO_PROJECT_HOME")) / filepath

        # Fetch PIC catalogue from FTP server
        try:
            h5 = h5py.File(hyper_file, 'r')
        except:
            errorcode('message', 'Inuaguration: Welcome to the PLATO variability simulator!')
            print(f"Downloading mass-radius parameterisation file...")
            downloadFromFTP(hyper_file.name, hyper_file.parents[0], server='plato')

        # Open file
        h5 = h5py.File(hyper_file, 'r')
        all_hyper = h5['hyper_posterior'][:]
        h5.close()


    def indicate(M, trans, i):
        '''
        indicate which M belongs to population i given transition parameter
        '''
        ts = np.insert(np.insert(trans, n_pop-1, np.inf), 0, -np.inf)
        ind = (M>=ts[i]) & (M<ts[i+1])
        return ind


    def split_hyper_linear(hyper):
        """
        split hyper and derive c
        """
        c0, slope,sigma, trans = \
        hyper[0], hyper[1:1+n_pop], hyper[1+n_pop:1+2*n_pop], hyper[1+2*n_pop:]
        
        c = np.zeros_like(slope)
        c[0] = c0
        for i in range(1,n_pop):
                c[i] = c[i-1] + trans[i-1]*(slope[i-1]-slope[i])

        return c, slope, sigma, trans


    
    def piece_linear(hyper, M, prob_R):
        '''
        model: straight line
        '''
        c, slope, sigma, trans = split_hyper_linear(hyper)
        R = np.zeros_like(M)
        for i in range(4):
                ind = indicate(M, trans, i)
                mu = c[i] + M[ind]*slope[i]
                R[ind] = norm.ppf(prob_R[ind], mu, sigma[i])

        return R


    def ProbRGivenM(radii, M, hyper):
        '''
        p(radii|M)
        '''
        c, slope, sigma, trans = split_hyper_linear(hyper)
        prob = np.zeros_like(M)

        for i in range(4):
                ind = indicate(M, trans, i)
                mu = c[i] + M[ind]*slope[i]
                sig = sigma[i]
                prob[ind] = norm.pdf(radii, mu, sig)

        prob = prob/np.sum(prob)

        return prob


    def classification( logm, trans ):
        '''
        classify as four worlds
        '''
        count = np.zeros(4)
        sample_size = len(logm)

        for iclass in range(4):
                for isample in range(sample_size):
                        ind = indicate( logm[isample], trans[isample], iclass)
                        count[iclass] = count[iclass] + ind

        prob = count / np.sum(count) * 100.
        print('Terran %(T).1f %%, Neptunian %(N).1f %%, Jovian %(J).1f %%, Star %(S).1f %%' \
                        % {'T': prob[0], 'N': prob[1], 'J': prob[2], 'S': prob[3]})
        return None


    
    def Mpost2R(mass, unit='Earth', classify='No'):
        """
        Forecast the Radius distribution given the mass distribution.

        Parameters
        ---------------
        mass: one dimensional array
                The mass distribution.
        unit: string (optional)
                Unit of the mass. 
                Options are 'Earth' and 'Jupiter'. Default is 'Earth'.
        classify: string (optional)
                If you want the object to be classifed. 
                Options are 'Yes' and 'No'. Default is 'No'.
                Result will be printed, not returned.

        Returns
        ---------------
        radius: one dimensional array
                Predicted radius distribution in the input unit.
        """

        # mass input
        mass = np.array(mass)
        assert len(mass.shape) == 1, "Input mass must be 1-D."

        # unit input
        if unit == 'Earth':
                pass
        elif unit == 'Jupiter':
                mass = mass * mearth2mjup
        else:
                print("Input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default.")

        # mass range
        if np.min(mass) < 3e-4 or np.max(mass) > 3e5:
                print('Mass range out of model expectation. Returning None.')
                return None

        ## convert to radius
        sample_size = len(mass)
        logm = np.log10(mass)
        prob = np.random.random(sample_size)
        logr = np.ones_like(logm)

        hyper_ind = np.random.randint(low = 0, high = np.shape(all_hyper)[0], size = sample_size)	
        hyper = all_hyper[hyper_ind,:]

        if classify == 'Yes':
                classification(logm, hyper[:,-3:])


        for i in range(sample_size):
                logr[i] = piece_linear(hyper[i], logm[i], prob[i])

        radius_sample = 10.** logr

        ## convert to right unit
        if unit == 'Jupiter':
                radius = radius_sample / rearth2rjup
        else:
                radius = radius_sample 

        return radius



    def Mstat2R(mean, std, unit='Earth', sample_size=1000, classify = 'No'):	
        """
        Forecast the mean and standard deviation of radius given the mena and standard deviation of the mass.
        Assuming normal distribution with the mean and standard deviation truncated at the mass range limit of the model.

        Parameters
        ---------------
        mean: float
                Mean (average) of mass.
        std: float
                Standard deviation of mass.
        unit: string (optional)
                Unit of the mass. Options are 'Earth' and 'Jupiter'.
        sample_size: int (optional)
                Number of mass samples to draw with the mean and std provided.
        Returns
        ---------------
        mean: float
                Predicted mean of radius in the input unit.
        std: float
                Predicted standard deviation of radius.
        """

        # unit
        if unit == 'Earth':
                pass
        elif unit == 'Jupiter':
                mean = mean * mearth2mjup
                std = std * mearth2mjup
        else:
                print("Input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default.")

        # draw samples
        mass = truncnorm.rvs( (mlower-mean)/std, (mupper-mean)/std, loc=mean, scale=std, size=sample_size)	
        if classify == 'Yes':	
                radius = Mpost2R(mass, unit='Earth', classify='Yes')
        else:
                radius = Mpost2R(mass, unit='Earth')

        if unit == 'Jupiter':
                radius = radius / rearth2rjup

        r_med = np.median(radius)
        onesigma = 34.1
        r_up = np.percentile(radius, 50.+onesigma, interpolation='nearest')
        r_down = np.percentile(radius, 50.-onesigma, interpolation='nearest')

        return r_med, r_up - r_med, r_med - r_down



    def Rpost2M(radius, unit='Earth', grid_size = 1e3, classify = 'No'):
        """
        Forecast the mass distribution given the radius distribution.

        Parameters
        ---------------
        radius: one dimensional array
                The radius distribution.
        unit: string (optional)
                Unit of the mass. Options are 'Earth' and 'Jupiter'.
        grid_size: int (optional)
                Number of grid in the mass axis when sampling mass from radius.
                The more the better results, but slower process.
        classify: string (optional)
                If you want the object to be classifed. 
                Options are 'Yes' and 'No'. Default is 'No'.
                Result will be printed, not returned.

        Returns
        ---------------
        mass: one dimensional array
                Predicted mass distribution in the input unit.
        """

        # unit
        if unit == 'Earth':
                pass
        elif unit == 'Jupiter':
                radius = radius * rearth2rjup
        else:
                print("Input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default.")


        # radius range
        if np.min(radius) < 1e-1 or np.max(radius) > 1e2:
                print('Radius range out of model expectation. Returning None.')
                return None



        # sample_grid
        if grid_size < 10:
                print('The sample grid is too sparse. Using 10 sample grid instead.')
                grid_size = 10

        ## convert to mass
        sample_size = len(radius)
        logr = np.log10(radius)
        logm = np.ones_like(logr)

        hyper_ind = np.random.randint(low = 0, high = np.shape(all_hyper)[0], size = sample_size)	
        hyper = all_hyper[hyper_ind,:]

        logm_grid = np.linspace(-3.522, 5.477, 1000)

        for i in range(sample_size):
                prob = ProbRGivenM(logr[i], logm_grid, hyper[i,:])
                logm[i] = np.random.choice(logm_grid, size=1, p = prob)

        mass_sample = 10.** logm

        if classify == 'Yes':
                classification(logm, hyper[:,-3:])

        ## convert to right unit
        if unit == 'Jupiter':
                mass = mass_sample / mearth2mjup
        else:
                mass = mass_sample

        return mass



    def Rstat2M(mean, std, unit='Earth', sample_size=1e3, grid_size=1e3, classify = 'No'):	
        """
        Forecast the mean and standard deviation of mass given the mean and standard deviation of the radius.

        Parameters
        ---------------
        mean: float
                Mean (average) of radius.
        std: float
                Standard deviation of radius.
        unit: string (optional)
                Unit of the radius. Options are 'Earth' and 'Jupiter'.
        sample_size: int (optional)
                Number of radius samples to draw with the mean and std provided.
        grid_size: int (optional)
                Number of grid in the mass axis when sampling mass from radius.
                The more the better results, but slower process.
        Returns
        ---------------
        mean: float
                Predicted mean of mass in the input unit.
        std: float
                Predicted standard deviation of mass.
        """
        # unit
        if unit == 'Earth':
                pass
        elif unit == 'Jupiter':
                mean = mean * rearth2rjup
                std = std * rearth2rjup
        else:
                print("Input unit must be 'Earth' or 'Jupiter'. Using 'Earth' as default.")

        # draw samples
        radius = truncnorm.rvs( (0.-mean)/std, np.inf, loc=mean, scale=std, size=sample_size)	
        if classify == 'Yes':
                mass = Rpost2M(radius, 'Earth', grid_size, classify='Yes')
        else:
                mass = Rpost2M(radius, 'Earth', grid_size)

        if mass is None:
                return None

        if unit=='Jupiter':
                mass = mass / mearth2mjup

        m_med = np.median(mass)
        onesigma = 34.1
        m_up = np.percentile(mass, 50.+onesigma, interpolation='nearest')
        m_down = np.percentile(mass, 50.-onesigma, interpolation='nearest')

        return m_med, m_up - m_med, m_med - m_down
        
