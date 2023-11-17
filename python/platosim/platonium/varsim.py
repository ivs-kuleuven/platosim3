#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script is an integrated part of PlatoSim's toolkit PLATOnium. Given a star and exoplent
this script creates a synthetic stellar and exoplanet variability model that can be used
directly as input for PlatoSim. A bolometric PLATO passband correction is applied to each 
photometric amplitude signal using synthetic high-resolution PHOENIX spectra (FGKM stars) or 
medium-resolution ATLAS9 spectra (OBAF stars). The current implemented models are:

Solar-like stars (F5-K7 dwarf and subgiants):
  - Granulation noise 
  - Convection driven oscillations (p-modes)
  - Stellar activity (spot modulations)

Other types of stars:
  - Cepheids/RR Lyrae  (fromFile)
  - gamma-Doradus star (analytic, fromFile)
  - delta-Scuti star   (analytic)

Exoplanet phase curve variations:
  - Transits               (using BATMAN)
  - Occultations           (using SPIDERMAN)
  - Limb Darkening         (using LDTk)
  - Doppler beaming        (using PyAstronomy)
  - Ellipsoidal distortion (using PyAstronomy)

User examples:
  $ varsim --star Sun --planet hotJupiter --quarter 1-8 -p
  $ varsim --star_params 1 1 5800 4.5 0 --planet_params 1 10 0 90 0 1 1 -o </path/to/varsource.txt>
"""

# Built-in
import os
import sys
import glob
import random
import argparse
import datetime
import warnings
from pathlib import Path

# PlatoSim standard
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.ndimage import median_filter
from scipy.interpolate import make_interp_spline
from astropy.io import fits
from astropy import constants as c
from astropy import units as u

# PLATOnium standard
import natsort
import batman
from ldtk import LDPSetCreator, TabulatedFilter

# PlatoSim functions
import platosim.plot      as pt
import platosim.utilities as ut
from platosim.utilities import errorcode
from platosim.spectrum  import Spectrum
from platosim.varsource import (StellarFlares,
                                StellarSpots,
                                SolarLikeOscillator,
                                GravityOscillator,
                                SurfaceModulations,
                                SMBHB,
                                PlanetMRforecast,
                                DopplerBeaming,
                                EllipsoidalDistortion)


#==============================================================#
#                         BEGIN CLASS                          #
#==============================================================#


class VarSim(object):

    """Class to generate noise-less light curves.

    NOTE Input parameters with a physical unit need an astropy.units attached.
    """
    
    def __init__(self, args):

        # CONSTANTS
        
        self.Teff_sun = 5777.  # [K]
        
        
        # I/O SETTINGS

        # Parameters in {True, False, None}
        self.plot  = args.plot
        self.seed  = args.seed
        self.star  = args.star
        self.ofile = args.ofile
        self.mocka = args.mocka
        self.binary = args.binary
        self.planet = args.planet
        #self.phase_curve = args.phase_curve TODO
        self.starID = None

        # Prepare pandas series for parameters
        self.df = pd.Series()
        
        # Random number BitGenerator (PCG64)
        if not self.seed:
            self.rng = np.random.default_rng()
        else:
            self.rng = np.random.default_rng(self.seed)

        # Verbosity (a.k.a log level) -> Identical to PlatoSim usage
        # verbose = 0: Cluster mode: Disabling print and warnings, and no log files are saved
        # verbose = 1: Default mode: Print details to bash but do not save log files
        # verbose = 3: Debug mode  : Print details to bash and saves all log files
        if args.verbose == 0:
            self.verbose = 0
            # Turn off warnings
            warnings.filterwarnings("ignore")
        elif args.verbose is None or args.verbose == 1:
            self.verbose = 1
            # Turn off warnings
            warnings.filterwarnings("ignore")
        else:
            self.verbose = 3

        # Print software name
        if self.verbose > 0:
            errorcode('software', '\nVariable Source Simulator\n')
            
        # Data path for VarSim
        self.idir = str(Path(os.getenv("PLATO_PROJECT_HOME")) / 'inputfiles/data_varsim')

        # Output file
        if self.ofile is not None:
            self.ofile = Path(self.ofile).resolve()
            
        # Add latex font if catalogue is saved
        if self.ofile is None:
            from platosim.matplotlibrc import setup
            setup()
        else:
            from platosim.matplotlibrc import latex
            latex()

        # Data (download) for usage of varsim
        if not Path(self.idir + '/passband_plato.txt').is_file():
            errorcode('message', 'Inuaguration: Welcome to the VarSim!')
            print(f'Downloading a few prerequisite files')
            ut.downloadFromFTP('passband_plato.txt',       self.idir)
            ut.downloadFromFTP('passband_cheops.txt',      self.idir)
            ut.downloadFromFTP('passband_tess.txt',        self.idir)
            ut.downloadFromFTP('passband_kepler.txt',      self.idir)
            ut.downloadFromFTP('varsim_meunier19a_t1.txt', self.idir)
            ut.downloadFromFTP('varsim_mainFitsBiSON.txt', self.idir)

            
        # OBS PARAMETERS

        # Cadence [d]
        if args.cadence:
            cadence = args.candece / 86400.
        else:
            cadence = 25 / 86400.

        # Parsing "quarter" overwrites "time"
        # NOTE Default is Q1 (t=0 and 90 days duration)
        if args.quarter:
            Q = ut.convertQuarterRange(args.quarter)
            if len(Q) == 1:
                numQuarters = 1
            elif len(Q) > 1:
                numQuarters = Q[1] - Q[0] + 1
            else:
                errorcode('error', 'Wrong input format of "quarter"!')        
            timeStart = round(ut.quarter() * (Q[0]-1) )
            timeDur   = round(ut.quarter() * numQuarters)
        elif args.time:
            timeStart = 0
            timeDur   = args.time
        else:
            timeStart = 0
            timeDur   = 90
        
        # Time points (ensure even number of time points)
        time = np.arange(0, timeDur, cadence)
        if len(time) % 2 != 0:
            time = np.arange(0, timeDur + cadence, cadence)

        # Store parameters
        self.time      = time * u.d
        self.cadence   = cadence * u.d
        self.timeDur   = timeDur * u.d
        self.timeStart = timeStart * u.d

        # Load the instrument passband
        if args.inst: self.intrument = args.inst
        else: self.instrument = 'plato'
        passband = pd.read_csv(self.idir + f'/passband_{self.instrument}.txt', comment='#')
        self.wvl_tele = passband.wavelength.to_numpy() * u.nm  # Wavelengths [nm]
        self.tra_tele = passband.passband.to_numpy()           # Normalized transmission
        
        # Define pandas data frame to store all signals
        self.lc = pd.DataFrame(data=self.time.to('s').value, columns=['time'])
        
        if self.verbose > 0:
            print(f'Simulating time series for : ' +
                  f'{len(self.time)} x {self.cadence.to("s")} ({self.timeDur.value} days)')
            print(f'Simulating {self.instrument} bandpass  : ' +
                  f'{self.wvl_tele[0]} - {self.wvl_tele[-1]}')




            
    #--------------------------------------------------------------#
    #                     BENCHMARK STARS/PLANETS                  #
    #--------------------------------------------------------------#

    
    def load_star(self, source):

        """Function to load benchmark stars used by varsim.

        NOTE in the following the astropy-units are required, however,
        the the choice of units (e.g. seconds vs. hours) are optional. 

        Parameters
        ----------
        source : str
            Stellar effective temperature [astropy.units]

        Returns
        -------
        Stellar parameters {M, R, Teff, logg, Z}
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

        if source == 'Ceph':
            M = 2.0 * u.M_sun
            R = 2.0 * u.R_sun
            Teff = 7500 * u.K
            logg = 4.0
            Z    = 0.0

            
        return M, R, Teff, logg, Z




    def load_exoplanet(self, source):

        """Module to load benchmark exoplanets to be used by varsim.

        NOTE in the following the astropy-units are required, however,
        the the choice of units (e.g. seconds vs. hours) are optional. 

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
        



    
    #--------------------------------------------------------------#
    #                       STELLAR PARAMETERS                     #
    #--------------------------------------------------------------#

    
    def stellar_source(self):

        """Select the stellar paramters.
        """

        if self.verbose > 0:
            errorcode('module', '\nStellar parameters\n')
        
        # Check star source or use Sun as default
        if args.star:
            self.star_source = args.star
        else:
            self.star_source = 'Sun'
                    
        # If stellar parameters are parsed
        if args.star_params:            
            star = args.star_params[0]
            self.M    = star[0] * u.M_sun
            self.R    = star[1] * u.R_sun
            self.Teff = star[2] * u.K
            self.logg = star[3]
            self.Z    = star[4]            
            
        # If Gaia ID is parsed
        elif self.mocka:
            self.star_source = self.df.spec
            self.Teff = self.df.Teff * u.K
            self.logg = self.df.logg
            self.Z    = self.df.Z
            self.R    = self.df.R * u.R_sun
            self.M    = self.df.M * u.M_sun
            self.L    = self.df.L * u.L_sun

        # Else use a benchmark star is parsed
        else:
            star = self.load_star(self.star_source)
            self.M    = star[0] 
            self.R    = star[1]
            self.Teff = star[2]
            self.logg = star[3]
            self.Z    = star[4]

        # Use theoretical logg if user has set logg=0
        if self.logg == 0:
            g = c.G.cgs * self.M.cgs / self.R.cgs**2
            self.logg = np.log10(g.value)
        
        # Derive theoretical luminosity
        if not self.mocka:
            self.L = self.R.value**2 * (self.Teff.value/self.Teff_sun)**4 * u.L_sun

        # Return parameters
        self.df['L_Lsun'] = self.L.to('L_sun').value
        self.df['M_Msun'] = self.M.to('M_sun').value
        self.df['R_Rsun'] = self.R.to('R_sun').value
        self.df['Teff_K'] = self.Teff.value
        self.df['logg']   = self.logg
        self.df['Z']      = self.Z
        self.df['alpha']  = 0.0
            
        # Print available stellar model parameters
        if self.verbose > 0:
            print(f"Spectral type   : {self.star_source}")
            print(f"Stellar Teff    : {self.df.Teff_K:4.0f}")
            print(f"Surface gravity : {self.df.logg:.2f} dex")
            print(f"Stellar [M/Fe]  : {self.df.Z:.3f} dex")
            print(f"Stellar radius  : {self.df.R_Rsun:.2f}")
            print(f"Stellar mass    : {self.df.M_Msun:.2f}")
            print(f"Luminosity      : {self.df.L_Lsun:.2f}")


            

        
    def binary_source(self):

        """Select the stellar paramters of a binary system.
        """

        if self.verbose > 0:
            errorcode('module', '\nBinary sources\n')

        


            
    def stellar_spectrum(self):

        """Calculates the bolometric correction from high-res spectra.
        
        This function uses the model grid method described in Sarkar+2018:
        https://academic.oup.com/mnras/article/481/3/2871/5092616

        # NOTE to compare theo L while using PhoenixAtmos, divide with np.pi
        """

        if self.verbose > 0:
            errorcode('module', '\nStellar spectrum\n')
        
       # Load parameters
        wvl_tele  = self.wvl_tele.to('AA').value
        tran_tele = self.tra_tele
        Teff = self.Teff.value
        logg = self.logg
        Z    = self.Z
        R    = self.R
        L    = self.L
        
        # Initialise class for synthetic spectra
        spec = Spectrum(verbose=self.verbose)

        # Select temeperature ranges
        if Teff <= 12200:
            library = 'PhoenixHiRes'
            if Teff < 7000: dT = 50
            else: dT = 100
        else:
            library = 'Atlas9'
            if Teff < 13000: dT = 125
            else: dT = 500

        # Make sure parameters are within limits
        Teff_lower, logg, Z, alpha = spec.nearest_parameters(Teff-dT, logg, Z, 0, library)
        Teff_upper, logg, Z, alpha = spec.nearest_parameters(Teff+dT, logg, Z, 0, library)

        # Fetch high resolution spectra
        if Teff <= 12200:
            wvl1_in, flux1_in = spec.getPhoenixHiResFITS(Teff_lower, logg, Z, alpha=0)
            wvl2_in, flux2_in = spec.getPhoenixHiResFITS(Teff_upper, logg, Z, alpha=0)
        else:
            wvl1_in, flux1_in = spec.getAtlasFITS(Teff_lower, logg, Z, alpha=0)
            wvl2_in, flux2_in = spec.getAtlasFITS(Teff_upper, logg, Z, alpha=0)
        
        # Cut off IR part TODO delete?
        # wvl_max  = 13000
        # dex      = np.where(wvl1_in > wvl_max)
        # wvl1_in  = np.delete(wvl1_in, dex)      # [AA]
        # wvl2_in  = np.delete(wvl2_in, dex)      # [AA]
        # flux1_in = np.delete(flux1_in, dex)     # [ergs/s/cm2/AA]
        # flux2_in = np.delete(flux2_in, dex)     # [ergs/s/cm2/AA]
        
        # Check that the interpolation between the two SEDs can be done
        if len(wvl1_in) != len(wvl2_in):
            errorcode('error', 'Spectra are not of same size! Check interpolation')

        # Create SED for star by interpolating nearby absolutely calibrated spectra
        wvl  = (wvl1_in + wvl2_in) / 2.
        flux = (flux1_in + (Teff-Teff_lower) * (flux2_in-flux1_in) *
                float(Teff_upper-Teff_lower)**-1)

        # Measure bolometric luminosity from SED [ergs/s]
        L1_bolometric = 4*np.pi*(R.cgs.value)**2 * np.trapz(flux1_in, wvl1_in)
        L2_bolometric = 4*np.pi*(R.cgs.value)**2 * np.trapz(flux2_in, wvl2_in)
            
        # Luminosity amplitude gradient in passband
        dex_wvl_min = ut.findNearestIndex(wvl, wvl_tele[0])
        dex_wvl_max = ut.findNearestIndex(wvl, wvl_tele[-1])
        if dex_wvl_max - dex_wvl_min == 1:
            L1_passband = (4*np.pi * (R.cgs.value)**2 *
                           ( wvl1_in[dex_wvl_max] -  wvl1_in[dex_wvl_min]) *
                           (flux1_in[dex_wvl_max] + flux1_in[dex_wvl_min]) / 2.)
            L2_passband = (4*np.pi * (R.cgs.value)**2 *
                           ( wvl2_in[dex_wvl_max] -  wvl2_in[dex_wvl_min]) *
                           (flux2_in[dex_wvl_max] + flux2_in[dex_wvl_min]) / 2.)
        else:
            L1_passband = (4*np.pi * (R.cgs.value)**2 *
                           np.trapz(flux1_in[dex_wvl_min:dex_wvl_max],
                                     wvl1_in[dex_wvl_min:dex_wvl_max]))
            L2_passband = (4*np.pi * (R.cgs.value)**2 *
                           np.trapz(flux2_in[dex_wvl_min:dex_wvl_max],
                                     wvl2_in[dex_wvl_min:dex_wvl_max]))
        
        # Bolometric cofficient
        self.bol_coeff = (ut.diff(L2_passband,   L1_passband) /
                          ut.diff(L2_bolometric, L1_bolometric))
        
        # Consistnecy check
        if self.verbose > 0:
            Lum = 4*np.pi*(R.cgs.value)**2 * np.trapz(flux, wvl)
            print(f'Theoretical luminosity : {L.to("erg/s"):.3e}')
            print(f'Synthetic   luminosity : {Lum * u.erg/u.s:.3e}')            
            print(f'Bolometric coefficient : {self.bol_coeff:.4f}')

        # Return parameters
        self.df['BC'] = self.bol_coeff
            
        # Rebinned spectrum TODO delete?
        wvl_equi = np.arange(wvl_tele[0], wvl_tele[-1], 1)
        wvl_equi, flux_equi = ut.rebin3(wvl_equi, wvl, flux)

        # # Absolute bolometric magnitude of star TODO delete?
        # M_bolometric = round(ut.diff(L2_bolometric, L1_bolometric) /
        #                      ut.diff(Teff_upper, Teff_lower), 2)
        # M_passband   = round(ut.diff(L2_passband, L1_passband) /
        #                      ut.diff(Teff_upper, Teff_lower), 2)
        # # Consistency check [mag]
        # if self.verbose > 0:
        #     print(f'Absolute magnitude bol : {round(M_bolometric, 4)}')
        #     print(f'Absolute magnitude lam : {round(M_passband,   4)}')            
        
        # PROLOGUE

        # Plot interpolation
        if args.plot:
            pt.plotSED(wvl, wvl1_in, wvl2_in, wvl_equi,
                       flux, flux1_in, flux2_in, flux_equi,
                       Teff, Teff_upper, Teff_lower)




            
    #--------------------------------------------------------------#
    #                   MODELS OF SOLAR-LIKE STARS                 #
    #--------------------------------------------------------------#


    def solar_granosc(self):

        """Model convection driven oscillations.
        """

        if self.verbose > 0:
            errorcode('module', '\nSolar-like oscillations\n')

        # Initialize and prepare model input
        params = [self.Teff, self.R, self.M, self.L]
        model  = SolarLikeOscillator(self.time, params, self.idir, seed=False)

        # Default scaling relations
        if args.gran is None:
            args.gran = 'Kallinger2014'
        if args.puls is None:
            args.puls = 'Corsaro2013'
        
        # Model granulation
        if not args.gran == None:
            params_gran = model.init_granulation(scaling=args.gran)
            self.lc['gran'] = model.eval_granulation() * self.bol_coeff
            self.df['A_gran_ppm'] = ut.rootMeanSquare(self.lc.gran)
            if self.verbose > 0:
                print(f'Granulation amplitude : {self.df.A_gran_ppm:.2f} ppm')

        # Model stochastic oscillations
        if not args.puls == None:
            params_puls = model.init_oscillations(scaling=args.puls)
            self.params_puls = params_puls[:3]
            self.lc['puls'] = model.eval_oscillations() * self.bol_coeff
            self.df['A_puls_ppm']   = ut.rootMeanSquare(self.lc.puls)
            self.df['numax_muHz']   = self.params_puls[0]
            self.df['deltanu_muHz'] = self.params_puls[0]
            if self.verbose > 0:
                print(f'Pulsation   amplitude : {self.df.A_puls_ppm:.2f} ppm')
                print(f'Frequency   nu_max    : {params_puls[0]:.2f} microHz')
                print(f'Splitting   Delta_nu  : {params_puls[1]:.2f} microHz')

        # Fix that the granulation time series is some time 1 time point shorter
        if args.gran and args.puls:
            if len(self.lc.gran) != len(self.lc.puls):
                self.lc.gran = np.append(self.lc.gran, self.lc.gran.mean())
                
        # Plot rsults
        if self.plot:
            pt.plot_amplitude_spectrum(self.lc, numax=params_puls[0])


    

        
    def solar_spots(self):

        """Model stellar spot modulations.
        """

        if self.verbose > 0:
            errorcode('module', '\nSolar-like spot modulations\n')

        # Use a random uniform distribution
        # NOTE secure a lower misalignment for planetary systems
        # Else use cos(i_star) uniform between 0-90 deg
        if self.planet or self.planet_params or self.kul20:
            incl = self.rng.uniform(85, 90)
        else:
            incl = None

        # Initialise model
        model = StellarSpots(seed=self.seed)
            
        # Generate model
        lc, params = model.evaluate(teff=self.Teff.value,
                                    time=self.time.to('d').value,
                                    dur=self.timeDur.to('d').value,
                                    cadence_hours=self.cadence.to('h').value,
                                    incl=incl)
        
        # print them to screen          
        if self.verbose:
            print(f'B-V           : {params[0]:.3f} mag')
            print(f"log R'_HK     : {params[1]:.3f}")
            print(f'Activity rate : {params[2]:.3f} solar')
            print(f'Rot. period   : {params[3]:.3f} days')
            print(f'Min. period   : {params[4]:.3f} days')
            print(f'Max. period   : {params[5]:.3f} days')
            print(f'Cycle period  : {params[6]:.3f} years')
            print(f'Cycle overlap : {params[7]:.3f} years')
            print(f'Max. latitude : {params[8]:.3f} deg')
            print(f'Inclination   : {params[9]:.3f} deg')

        # Store global variables
        self.lc['spot'] = lc * 1e6
        self.df['B_V']      = params[0]
        self.df['logR_HK']  = params[1]
        self.df['AR_ARsun'] = params[2]
        self.df['Prot_day'] = params[3]
        self.df['Pmin_day'] = params[4]
        self.df['Pmax_day'] = params[5]
        self.df['Pcyc_day'] = params[6]
        self.df['Povl_day'] = params[7]
        self.df['Lmax_deg'] = params[8]
        self.df['I_deg']    = params[9]

        # Plot model
        if self.plot: model.plot()
        


        
        
    def solar_flares(self): # TODO under construction

        """Model solar flares.
        """

        # Initialise model
        time  = self.time.to('d').value
        model = StellarFlares(time, seed=self.seed)

        # Run simple model for now
        model.initToyModelBeta0()

        # Return model [mag -> ppm]
        mag = model.evaluate(plot=self.plot)
        self.lc['flux'] = ut.fromMagToFlux(mag) * self.bol_coeff




    
    #--------------------------------------------------------------#
    #                         OTHER PULSATORS                      #
    #--------------------------------------------------------------#
    

    def star_roap(self):

        """Generate light curve for roAp stars.

        This class is of BAF stars with so-called surface spots.
        """
        
        # Start script
        if self.verbose > 0:
            errorcode('module', '\nRotational variable (roAp)\n')

        # Initialize class
        time  = self.time.to('d').value
        model = SurfaceModulations(time, seed=False)

        # Prepare model parameters
        params = model.initToyModel()
        if self.verbose > 0:
            print(f'Rotational period   : {round(params[0],3)} days')
            print(f'Random phase offset : {round(np.rad2deg(params[1]),1)} deg')
            print(f'Relative amplitude  : {round(params[2],3)}')
            print(f'Scaled amplitude    : {round(params[3],3)}')
        
        # Return model
        self.lc['roap'] = model.evaluate(plot=args.plot)
        self.df['Prot_day'] = params[0]
        self.df['dphi_rad'] = params[1]
        self.df['Arel']     = params[2]
        self.df['scale']    = params[3]



        
    
    def star_gdor(self):

        """Generate light curves for gamma-Dor stars.

        Notes 
        -----
        This function uses the "varsouce.GravityOscillator" class.
        This class provide two model generations:
        1) Toy model using a characteristic power of 2.2
        2) Draw (freq, ampl, phase) from Kepler observations by
           using the flag "--puls gang2020".
        """

        # Start script
        if self.verbose > 0:
            errorcode('module', '\nPulsator: gamma-Dor (g-modes)\n')

        # Initialize and prepare model input
        time  = self.time.to('d').value
        model = GravityOscillator(time, power=2.2, seed=None) # TODO code is wrong?
        
        # Check if a file with pulsations are parsed
        if args.puls == 'gang2020':
            model.initGang2020(self.idir, starID=self.starID)
        else:
            model.initToyModel([0.5, 3], [0.5, 2.5])

        # Return model [mag -> ppm]
        mag = model.evaluate(plot=args.plot)
        self.lc['flux'] = ut.fromMagToFlux(mag) * self.bol_coeff



            
        
    def star_dsct(self): # TODO

        """Generate light curves for delta-Scuti stars.

        Notes 
        -----
        This function used the "gravity_oscill" utility with a characteristic
        power of 1.0 for these g-mode pulsators.
        """

        # Start script
        if self.verbose > 0:
            errorcode('module', '\nPulsator: delta-Scuti (g-modes)\n')

        # Initialize and prepare model input
        time  = self.time.to('d').value
        model = GravityOscillator(time, power=1.0, seed=self.seed)
        
        # Check if a file with pulsations are parsed
        model.initToyModel([1, 30], [10, 30])

        # Return model [mag -> ppm]
        mag = model.evaluate(plot=args.plot)
        self.lc['flux'] = ut.fromMagToFlux(mag) * self.bol_coeff





    def star_ceph(self): # TODO

        """Generate ligth curve for Cepheid and RR Lyrae stars.

        This function uses precomputed models of Cepheids and RR Lyrae
        stars to generate the light curve from their harmonics. For now
        this function will randomly select a star.
        """

        if self.verbose > 0:
            errorcode('module', '\nClassical pulsators (RR Lyrae & Cepheids)\n')

        # Select a random object from the list and load Fourier data
        filename = 'varsource_RRLyrae_Cepheid_bodi2023'
        try:
            filenames = glob.glob(f'{self.idir}/{filename}/*.fou')
            starfile = random.choice(filenames)
        except:
            zipfile = f'{filename}.zip'
            errorcode('message', 'Classic, I like your style!')
            print(f'Downloading {zipfile} files..')
            ut.downloadFromFTP(filename=zipfile, outputDir=self.idir, server='plato')
            os.system(f'unzip {self.idir}/{zipfile} -d {self.idir}')
            os.system(f'rm {self.idir}/{zipfile}')
            print('')

        # Load file with harmonics
        filenames = glob.glob(f'{self.idir}/{filename}/*.fou')
        starfile  = random.choice(filenames)
        fourier   = np.loadtxt(starfile)    
        if self.verbose > 0:
            print(f'Using file {starfile} with frequencies:')
            print(fourier)
            
        # Convert units of input parameters
        time = self.time.value
        flux = np.zeros_like(time)


        # Generate the lightcurve in the dataframe
        components = len(np.array(fourier))
        for i in range(components):
            flux += fourier[i,1] * np.sin((2*np.pi*fourier[i,0] * time) + fourier[i,2])

        # Save magnitude to list
        lc = pd.DataFrame(data = {'time':time, 'flux':flux})
        mag = - 2.5 * np.log10(lc.flux + 1)

        #self.lc['flux'] = ut.fromMagToFlux(mag) * self.bol_coeff
        

        # plot light curve
        if args.plot:
            plt.figure(figsize=(10, 5))
            plt.plot(time, mag*1e3, 'm-')
            plt.xlabel('Time [d]')
            plt.ylabel(r'$\delta m$ [mmag]')
            plt.xlim(np.min(time), np.max(time))
            plt.tight_layout()
            plt.show()
        



        
    #--------------------------------------------------------------#
    #                          BINARY SYSTEMS                      #
    #--------------------------------------------------------------#
    
    
    def binary_smbh(self):

        """Function to generate a SMBH binary light curve.

        A Super Massive Black Hole (SMBH) binary system consist of several 
        components, for which this model includes two of the effects:
          1) The doppler boosting
          2) The gravitational lensing effect
          3) The stochastic quasar variability
        """

        if self.verbose > 0:
            errorcode('module', '\nSMBH binary\n')

        # Set the stellar source entry
        self.star_source = 'SMBHB'
            
        # Fetch time array
        time  = self.time.to('d').value
        model = SMBHB(time, seed=self.seed)

        # Fetch model parameters
        P, A_beam, A_lens, phi, tmax, tdur = model.initToyModel()

        if self.verbose:
            print(f'Model parameters of toy model:')
            print(f'Orbital period     : {P:.3f} day')
            print(f'Beaming amplitude  : {A_beam*1e3:.3f} mmag')
            print(f'Lensing amplitude  : {A_lens*1e3:.3f} mmag')
            print(f'Lens time maximum  : {tmax:.3f} day')
            print(f'Lens time duration : {tdur:.3f} day')
        
        # Get model
        flux, flux_beam, flux_lens = model.evalToyModel(P, A_beam, A_lens, phi, tmax, tdur) 

        # plot light curve
        if self.plot: model.plot()
                
        # Return light curve
        self.lc['flux']      = flux
        self.lc['flux_beam'] = flux_beam
        self.lc['flux_lens'] = flux_lens

        # Return parameters
        self.df['P_day']         = P
        self.df['A_beam_mag']    = A_beam
        self.df['A_lens_mag']    = A_lens
        self.df['tmax_lens_day'] = tmax
        self.df['tdur_lens_day'] = tdur        



        
    #--------------------------------------------------------------#
    #                          PLANET MODELS                       #
    #--------------------------------------------------------------#
        
        
    def ldc(self):

        """Compute the Limb Darkening (LD) coefficients.
        """
        
        #  Convert input parameters
        wvl_tele  = self.wvl_tele.value  # [nm]
        tran_tele = self.tra_tele        # [0-1]

        # Interpolate (piecewise cubic) into higher resolution grid
        grid_no  = 1000
        wvl_int  = np.linspace(wvl_tele[0], wvl_tele[-1], grid_no)
        passband = make_interp_spline(wvl_tele, tran_tele, k=3)
        tran_int = passband(wvl_int)

        # Create passband object
        filters = [TabulatedFilter('plato', wvl_int, tran_int)]

        # Create instance of class
        # NOTE uncertainties are vital for the software to work!
        sc = LDPSetCreator(teff=(self.Teff.value, 50),
                           logg=(self.logg, 0.20),
                           z=(self.Z, 0.05),
                           filters=filters)

        # Create the limb darkening profiles
        ps = sc.create_profiles()

        # Estimate quadratic law coefficients
        # Take care of occations when LDTk fails
        try:
            u, _ = ps.coeffs_qd(do_mc=True)
        except:
            self.ldc = [0.430, 0.170]
            errorcode('warning', 'LD coefficients failed for ' +
                      f'(Teff, logg, Z) = ({self.Teff}, {self.logg}, {self.Z}')
        else:
            self.ldc = u[0]
                    

            
    def planet_model(self): # TODO Include into class

        """Calculation of exoplanet model parameters.

        The following equations are from chapters in Seager et al. (2010): "Exoplanets":
        Winn (2014)             : https://arxiv.org/pdf/1001.2010.pdf
        Murray & Correia (2011) : https://arxiv.org/abs/1009.1738v2
        NOTE Both "a" and "Rs" are given in units of Rstar.

        Assumptions
        -----------
        - The calculations in this code block are under the following assumptions that
        that eclipses are centered around conjunction. This is not valid for extremely
        eccentric and close-in orbits with grazing eclipses. However, non-grazing close
        in orbits are still valid.
        - The time seperation between transit and occultation in the following is a
        first order approximation in "e" by integrating "dt/dF".
        """

        # Stellar parameters
        time = self.time.to('d')
        Ms   = self.M.to('kg')
        Rs   = self.R.to('m')
        Teff = self.Teff.to('K')

        
        # LOAD PLANET MODEL
        
        if args.planet == 'random': # TODO implement as part of KUL20 mode!

            # Select benchmark planets
            name_benchmark = ['Earth-like', 'Neptune-like', 'Jupiter-like']
            Rp_benchmark = [1.0, 3.9, 11.2]
            Rp = np.random.choice(Rp_benchmark) * u.R_earth
            # Log-normal distribution of planet radii 
            #Rp = np.random.lognormal(5., 1.) * u.R_earth
            #if Rp.to('R_jup').value < 0.05: Rp = 0.05 * u.R_jup
            #if Rp.to('R_jup').value > 8.50: Rp = 1.80 * u.R_jup

            # Draw mass from M-R forecaster
            # Paper: Chen, J., & Kipping, D. M. 2018, MNRAS, 473, 2753
            # Code: https://github.com/chenjj2/forecaster
            # NOTE forecasting only valid between (0.05 < Rp/R_jup < 8.5) 
            if self.verbose > 0: classifier = 'yes'
            else: classifier = 'no'
            Mp, _, _ = mr_forecast.Rstat2M(mean=Rp.to('R_jup').value, unit='Jupiter',
                                           std=0.01, sample_size=1000, grid_size=1000,
                                           classify=classifier)
            Mp = Mp * u.M_jup

            # Select random uniform period
            P  = np.random.uniform(0, 1) * u.yr

            # Make sure at least one transit fall within the duration of the observation
            t0 = np.random.uniform(0, self.timeDur.to('yr').value) * u.yr
            if P <= self.timeDur.to('yr'):
                t0 = t0 - self.timeDur.to('yr')

            # Simple for now
            e = 0

            # Unbiased uniform distribution 
            #i = np.arccos(np.random.uniform(0, 90/85-1)) * 180/np.pi * u.deg  # between 85-90 deg
            i = 90 * u.deg
            
            # Uniform orientation
            w = np.random.uniform(0, 360) * u.deg

            # Convert units 
            t0 = t0.to('s')
            P  = P.to('s')
            i  = i.to('rad')
            w  = w.to('rad')
            Rp = Rp.to('m')
            Mp = Mp.to('kg')
        
        elif args.planet_params is None:
            # Load exoplanet parameters [SI units]
            try:
                params = self.load_exoplanet(args.planet)
            except UnboundLocalError:
                errorcode('error', 'No planet with that name! Check --planet entry')
            else:
                t0 = params['t0'].to('s')
                P  = params['P'].to('s')
                e  = params['e']
                i  = params['i'].to('rad')
                w  = params['w'].to('rad')
                Rp = params['rp'].to('m')
                Mp = params['mp'].to('kg')
                xi = params['xi']
                Tn = params['Tn'].to('K').value
                dT = params['dT'].to('K').value
            
        else:
            # Load exoplanet parameters [SI units]
            params = args.planet_params[0]
            t0 = (params[0] * u.d).to('s')
            P  = (params[1] * u.d).to('s')
            e  = params[2]
            i  = (params[3] * u.deg).to('rad')
            w  = (params[4] * u.deg).to('rad')
            Rp = (params[5] * u.R_earth).to('m')
            Mp = (params[6] * u.M_earth).to('kg')
            xi = 0.
            Tn = 0.
            dT = 0.

            
        # ORBITAL DYNAMICS

        # Handy definitions
        x = np.sqrt(1 - e**2)   # Optimization constant
        k = Rp/Rs               # Radius constant

        # Semi-major axis (K3)
        a = (c.G * P**2 * (Ms + Mp) / (4*np.pi**2))**(1/3.)

        # Impact parameter: Winn (2014) Eq. 7 & 8
        b_tra = a*np.cos(i)/Rs * (1 - e**2)/(1 + e*np.sin(w))
        b_occ = a*np.cos(i)/Rs * (1 - e**2)/(1 - e*np.sin(w))

        # Transit and occultation times: Winn (2014) Eq. 14, 15 & 16
        ep = x/(1 + e*np.sin(w)) / u.rad
        em = x/(1 - e*np.sin(w)) / u.rad
        t_tra_tot = P/np.pi * np.arcsin( Rs/a * np.sqrt((1 + k)**2 - b_tra**2)/np.sin(i) ) * ep
        t_tra_ful = P/np.pi * np.arcsin( Rs/a * np.sqrt((1 - k)**2 - b_tra**2)/np.sin(i) ) * ep
        t_occ_tot = P/np.pi * np.arcsin( Rs/a * np.sqrt((1 + k)**2 - b_occ**2)/np.sin(i) ) * em
        t_occ_ful = P/np.pi * np.arcsin( Rs/a * np.sqrt((1 - k)**2 - b_occ**2)/np.sin(i) ) * em
        tau_tra = (t_tra_tot - t_tra_ful)/2.
        tau_occ = (t_occ_tot - t_occ_ful)/2.

        # Transits to occultations time seperation: Winn (2014) Eq. 33
        dt_c = P/2 * (1 + 4*e*np.cos(w)/np.pi)

        # Time of first (t0) central passage transit and occultation
        t0_tra_cen = t0
        t0_occ_cen = t0 + dt_c

        # Check for model applicability
        if b_tra > 1 - Rp/Rs and b_tra <= 1 + Rp/Rs:
            errorcode('warning', 'Planetary model consist of grazing eclipses!')
        elif b_tra > 1 + Rp/Rs:
            errorcode('warning', 'Planet model do not have any physical eclipses!')
            
        # Show parameters apce
        if self.verbose > 0:
            errorcode('module', '\nPlanet eclipse model')
            print('')
            print("Planet name        : {}".format(args.planet))
            print("Planet mass        : {:.2f}".format(Mp.to('M_earth')))
            print("Planet radius      : {:.2f}".format(Rp.to('R_earth')))
            print("Semimajor axis     : {:.2f} starRad".format(a.to('m')/Rs.to('m')))
            print("Eccentricity       : {:.3f}".format(e))
            print("Inclination        : {:.2f}".format(i.to('deg')))
            print("Arg. of periastron : {:.2f}".format(w.to('deg')))
            print("Orbital Period     : {:.2f}".format(P.to('d')))
            print("Time of emphemeris : {:.2f}".format(t0.to('d')))
            print("Total tra duration : {:.3f}".format(t_tra_tot.to('h')))
            print("Full  tra duration : {:.3f}".format(t_tra_ful.to('h')))
            print("In/Egress duration : {:.3f}".format(tau_tra.to('min')))
            print("Impact parameter   : {:.3f}".format(b_tra))
            # if self.phase_curve: TODO
            #     print('')
            #     print("\nTransit-to-Occultation time  : {:.3f}".format(dt_c.to('d')))
            #     print("Total occultation duration   : {:.3f}".format(t_occ_tot.to('h')))
            #     print("Full  occultation duration   : {:.3f}".format(t_occ_ful.to('h')))
            #     print("In/Egress occult. duration   : {:.3f}".format(tau_occ.to('min')))
            #     print("Impact parameter Occultation : {:.3f}\n".format(b_occ))
            
        # Store parameters
        self.df['Mp_Mearth'] = Mp.to('M_earth').value
        self.df['Rp_Rearth'] = Rp.to('R_earth').value
        self.df['a_Rstar']   = (a.to('R_sun')/Rs).value
        self.df['P_day']     = P.to('d').value
        self.df['t0_day']    = t0.to('d').value
        self.df['e']         = e
        self.df['i_deg']     = i.to('deg').value
        self.df['w_deg']     = w.to('deg').value
        self.df['u1']        = self.ldc[0]
        self.df['u2']        = self.ldc[1]

        # Store parameters
        self.Mp = Mp
        self.Rp = Rp
        self.t0 = t0
        self.P  = P
        self.a  = a
        self.e  = e
        self.i  = i
        self.w  = w
        self.t_tra_tot  = t_tra_tot
        self.t_occ_tot  = t_tra_tot
        self.t0_tra_cen = t0_tra_cen
        self.t0_occ_cen = t0_occ_cen
        self.dt_c = dt_c

            



    def planet_transit(self):

        """Model exoplanet transits.

        In the following the exoplanet transits are being modelled with Batman:
        https://lweb.cfa.harvard.edu/~lkreidberg/batman/quickstart.html

        NOTE: t0 and P can principly be anything as long as they are consistant. 
        Here we make sure to use consistent reference time unit.
        """

        # Limb darkening model options:
        if args.ldm: limbDarkModel = args.lmd
        else: limbDarkModel = 'quadratic'

        # Initialize batman model
        batman_params = batman.TransitParams()
        batman_params.t0        = self.t0.to('d').value
        batman_params.per       = self.P.to('d').value
        batman_params.a         = (self.a.to('m')/self.R.to('m')).value
        batman_params.ecc       = self.e
        batman_params.inc       = self.i.to('deg').value
        batman_params.w         = self.w.to('deg').value
        batman_params.rp        = (self.Rp.to('m')/self.R.to('m')).value
        batman_params.u         = self.ldc
        batman_params.limb_dark = limbDarkModel

        # Initializes transit model and extract light curve [ppm]
        model = batman.TransitModel(batman_params, self.time.value)
        self.lc['tran'] = (model.light_curve(batman_params) - 1) * 1e6

        # True anomaly at each time: This will be used in our custom models later
        self.nu = model.get_true_anomaly()

        # Time of periastron passage (calculated from t0)
        self.tau = model.get_t_periastron(batman_params)

        # Print to bash
        if self.verbose:
            print(f"Mid-transit depth  : {np.abs(np.min(self.lc.tran)):.1f} ppm")
            print(f"LD coefficients    : {self.ldc[0]:.3f}, {self.ldc[1]:.3f}")





    def planet_occultation(self): # TODO Test again for student project!
        """
        In the following the exoplanet transits are being modelled with Batman:
        https://spiderman.readthedocs.io/en/latest/index.html

        Spiderman models the secondary eclipses plus phase curves, however, not the
        transit hence Batman is needed still for this. The two codes are compatible
        in terms of parameter definitions, as they are developed in collaboration.

        NOTE the unit of "t0" and "P" can in principly be anything as long as they
        are consistant. Here we make sure to use consistent reference time unit.
        NOTE "n_layer" is the number of layers in the 2D "web" defining the integration
        grid. A value of 20 give an error less than 0.1 ppm.
        """

        # Initialize model and assign parameters
        # TODO brightness model should be loaded from file by user

        # import spiderman
        # exo_brightness_model = 'zhang'

        # # TODO The phase curve is now calculated using a blackbody model but it is also
        # # possible to use PHOENIX med-res spectra. This not work now - fix it!

        # spider_params = spiderman.ModelParams(brightness_model=exo_brightness_model, stellar_model='PHOENIX')
        #                                       #stellar_model= datapath + 'stellar_model.txt')

        # spider_params.t0    = t0.to('d').value
        # spider_params.per   = P.to('d').value
        # spider_params.a     = (a/Rs).value
        # spider_params.a_abs = a.to('AU').value
        # spider_params.ecc   = e
        # spider_params.inc   = i.to('deg').value
        # spider_params.w     = w.to('deg').value
        # spider_params.rp    = (Rp/Rs).value
        # spider_params.p_u1  = ldc[0]
        # spider_params.p_u2  = ldc[1]
        # spider_params.T_s   = Teff.value
        # spider_params.l1    = wvl_tele[0].to('m').value
        # spider_params.l2    = wvl_tele[-1].to('m').value
        # spider_params.n_layer = 100  # The number of layers in the 2D "web"

        # # TODO convolve with the acual bandpass!
        # # Use spiderman's weithing scheme with the instrument response

        # spider_params.filter = os.getcwd() + '/data/Passbands/response_plato_spiderman.txt'

        # # USE USER DEFINED BRIGHTNESS MODELS

        # # TODO Verify all the available brigthness models

        # # Spherical hamonics: Non-physical model
        # # NOTE there is nothing in the implementation to prevent negative surface fluxes!

        # if exo_brightness_model == 'spherical':
        #     spider_params.degree = degree
        #     spider_params.la0    = la0
        #     spider_params.lo0    = lo0
        #     spider_params.sph    = sph

        # # Zhang and Showman (2017): http://adsabs.harvard.edu/abs/2017ApJ…836…73Z
        # # A temperature map model based on semi-physical reproducing the main features
        # # of hot Jupiter phase-curves - offset hotspots. Called with “zhang”

        # if exo_brightness_model == 'zhang':
        #     spider_params.xi      = xi   # Ratio of radiative to advective timescale
        #     spider_params.T_n     = Tn   # Temperature of nightside
        #     spider_params.delta_T = dT   # Day-night temperature contrast

        # # Offset hotspot and Two sided planet:

        # if (exo_brightness_model == 'hotspot_b' or exo_brightness_model == 'hotspot_t'):
        #     spider_params.la0  = la0
        #     spider_params.lo0  = lo0
        #     spider_params.size = size
        #     spider_params.grid_size = grid_size

        #     if exo_brightness_model == 'hotspot_b':

        #         if spot_b is not None and p_b is not None:
        #             spider_params.spot_b = spot_b
        #             spider_params.p_b    = p_b
        #         if pb_d is not None and pb_n is not None:
        #             spider_params.pb_d = pb_d
        #             spider_params.pb_n = pb_n

        #     if exo_brightness_model == 'hotspot_t':
        #         spider_params.spot_T = spot_T
        #         spider_params.p_T    = p_T

        # # Two sided planet:

        # if exo_brightness_model == 'spherica':
        #     spider_params.degree = degree
        #     spider_params.la0    = la0
        #     spider_params.lo0    = lo0
        #     spider_params.sph    = sph

        # Contruct model light curve [ppm]
        #lc_occ = (spider_params.lightcurve(time.value) - 1) * 1e6
        lc_occu = np.zeros(len(self.time))
        self.lc['occu'] = lc_occu.tolist()




        
    def planet_beaming(self): # TODO Implement into class
        """
        Doppler beaming model.
        """

        # Central wavelength of PLATO bandpass [m]
        wvl_c = self.wvl_tele[np.argmax(self.tra_tele)]

        # Initialize and prepare model input
        model_beam  = DopplerBeaming()
        params_beam = {'wvl_c': wvl_c,
                       'Teff': self.Teff,
                       't0': self.t0,
                       'P': self.P,
                       'a': self.a,
                       'e': self.e,
                       'i': self.i,
                       'w': self.i,
                       'Mp': self.Mp,
                       'Ms': self.M}

        # Assign and calculate model
        model_beam.assignValue(params_beam)
        lc_beam, self.A_beam = model_beam.evaluate(self.time, self.nu)

        # Combine models
        self.lc['beam'] = lc_beam.tolist()






    def planet_ellipsoidal(self): # TODO Implement into class

        """Model ellipsoidal distortion.
        """

        # Initialize and prepare model input
        model_elli  = EllipsoidalDistortion()
        params_elli = {'g': 0.05,
                       'u': self.ldc[0],
                       't0': self.t0,
                       'P': self.P,
                       'e': self.e,
                       'i': self.i,
                       'w': self.w,
                       'a': self.a,
                       'Rs': self.R,
                       'Mp': self.Mp,
                       'Ms': self.M}

        # Assign and calculate model
        model_elli.assignValue(params_elli)
        lc_elli, self.A_elli = model_elli.evaluate(self.time, self.nu)

        # Combine models
        self.lc['elli'] = lc_elli.tolist()




    def plot_phase_curve(self):

        # Plot exoplanet model
        if (self.time[-1] >= self.P.to('d') + self.t0.to('d')):
            fig = plt.figure(figsize=(13, 10))
            # Adjust times
            lc_exo = self.lc['tran'] + self.lc['occu'] + self.lc['beam'] + self.lc['elli']
            pt.plot_orbital_phase_curve(fig, self.time.value,
                                        self.lc['tran'].to_numpy(),
                                        self.lc['occu'].to_numpy(),
                                        self.lc['beam'].to_numpy(),
                                        self.lc['elli'].to_numpy(),
                                        lc_exo.to_numpy(),
                                        self.t0.to('d').value,
                                        self.P.to('d').value,
                                        self.dt_c.to('d').value,
                                        self.t0_tra_cen.to('d').value,
                                        self.t_tra_tot.to('d').value,
                                        self.t0_occ_cen.to('d').value,
                                        self.t_occ_tot.to('d').value,
                                        self.A_beam, self.A_elli)
        elif (self.time[-1] < self.P.to('d') + self.t0.to('d')):
            errorcode('warning',
                      'No phase plot, time series is shorter than the orbital period!')





    #--------------------------------------------------------------#
    #                      PROLOGUE AND SAVING                     #
    #--------------------------------------------------------------#
        

    def run_prolog(self):
        
        if self.verbose > 0:
            errorcode('module', '\nPrologue\n')

        # SORT LIGHT CURVE

        # Convert time unit [d -> s]
        if args.quarter:
            self.lc['time'] += self.timeStart * 86400
        
        # Variability classes
        stars    = ['std', 'dSct', 'gDor', 'Ceph']
        binaries = ['SMBH']
            
        # Combine all signals for solar-like stars
        if (not self.star in stars and
            not self.binary in binaries and
            self.mocka is None):

            # Granulation and pulsation are additive
            self.lc['flux'] = np.zeros(len(self.lc.time))
            if 'gran' in self.lc:
                self.lc['flux'] += self.lc.gran
            if 'puls' in self.lc:
                self.lc['flux'] += self.lc.puls
            if 'spot' in self.lc:
                self.lc['flux'] += self.lc.spot

            # Convert to relative flux to multiply with transits
            self.lc['flux'] = self.lc['flux'] / 1e6 + 1 
                
            # Spots and transits are multiplicative
            if 'tran' in self.lc:
                self.lc['flux'] *= (self.lc.tran / 1e6 + 1)
                
            # Plot combined light curve
            if self.plot:
                fig, ax = pt.plot_final_lc(self.lc)
                plt.show()
                                
            
        # SAVE DATA
        
        if self.ofile:

            # Filenames
            ofile_parameters = self.ofile.parents[0] / f'{self.ofile.stem}_parameters.ftr'
            ofile_components = self.ofile.parents[0] / f'{self.ofile.stem}_components.ftr'
            
            if self.verbose:
                print(f'Saving file : {self.ofile}')
                print(f'Saving file : {ofile_parameters}')
                print(f'Saving file : {ofile_components}')

            # Convert to magnitude [mag]
            df = self.lc.flux.to_numpy() 
            dm = - 2.5 * np.log10(df)            
                
            # Save light curve
            data = np.transpose([self.lc['time'], dm])
            np.savetxt(self.ofile, data, fmt=['%.1f', '%.8f'])

            # Save parameter space
            self.df = self.df.to_frame().T
            self.df.to_feather(ofile_parameters)
            self.lc.to_feather(ofile_components)
            




    #--------------------------------------------------------------#
    #                         SOFTWARE MODES                       #
    #--------------------------------------------------------------#


    def mode_single(self):

        """Given stellar properties asign variable signal.
        """
        
        # Select star
        self.stellar_source()

        # Bolometric correction
        self.stellar_spectrum()

        # Activate spot modulation by default
        if args.spot is True or args.spot is None:
            args.spot = True
        
        # Include stellar variability
        if args.star == 'roAp':
            v.star_roap()
            
        elif args.star == 'gDor':
            v.star_gdor()

        elif args.star == 'dSct':
            v.star_dsct()

        elif args.star == 'Ceph':
            v.star_ceph()

        else:
            # Solar-like stars
            if args.star or args.star_params:
                if args.spot is True:
                    v.solar_spots()
                if not args.gran or not args.puls:
                    v.solar_granosc()

            # Include exoplanet
            if args.planet or args.planet_params or args.planet == 'random':
                v.ldc()
                v.planet_model()
                v.planet_transit()

                # For hot-Jupiters include phase curve TODO
                # if args.phase_curve:
                #     v.planet_occultation()
                #     v.planet_beaming()
                #     v.planet_ellipsoidal()
                #     if not args.kul20 and args.plot:
                #         v.plot_phase_curve()

        # Combine and save
        self.run_prolog()




        
    def mode_binary(self):

        """Given stellar properties asign variable signal.
        """
        
        # Select binary system
        self.binary_source()

        # Bolometric correction
        #self.stellar_spectrum()

        if args.star == 'EB':
            v.binary_system()
        
        elif args.binary == 'SMBH':
            v.binary_smbh()

        # Combine and save
        self.run_prolog()

        


        
    def mode_kul20(self):

        """Given stellar properties asign variable signal.
        """
        
        # Notes on flag "--kul20" -> used for KUL20
        # 0 -> Std/constant (2 hamonics)
        # 1 -> Gran, Puls
        # 2 -> Gran, puls, Spot
        # 3 -> Gran, Puls, Spots, Exo
        # x -> Constant stars is any other number x
        if args.kul20 == 0:
            args.star = 'roAp'
        if args.kul20 in (1, 2, 3):
            args.star          = 'Sun'
            args.planet_params = False
        if args.kul20 == 2:
            args.spot          = True
            args.planet        = False
        if args.kul20 == 3:
            args.spot          = True
            args.planet        = 'random'
        if not args.kul20 in (0, 1, 2, 3):
            args.kul20 = False

        # Add steps from default mode
        self.mode_default()




        
    def mode_mocka(self):

        """Given stellar properties asign variable signal.
        """

        # I/O EXTRA

        project, starType, starID, conFlag, odir = args.mocka[0]
        idir = Path(os.getenv('PLATO_WORKDIR')) / project / 'input'        
        odir = Path(odir).resolve()
        self.starID = int(starID)
        
        # Load feather
        df0 = pd.read_feather(idir /  'starcat_GaiaDR3_PlatoCS.ftr')
        ds0 = pd.read_feather(idir / f'starcat_GaiaDR3_PlatoCS_{starType}.ftr')
        
        # Output directory
        starDir = f'{self.starID}'.zfill(9)
        self.odir = odir / starDir
        self.odir.mkdir(parents=True, exist_ok=True)


        # QUERY STARS IN SUBFIELD

        # Select target star
        df_i = ds0.iloc[self.starID-1]

        # Fetch smaller region around target
        x = 45/3600.
        dc_i = df0[(df0.ra  > df_i.ra  - x) & (df0.ra  < df_i.ra  + x) &
                   (df0.dec > df_i.dec - x) & (df0.dec < df_i.dec + x)]

        # Find radial distance [arcsec] 
        dc_i['dis'] = ut.radialDistance(df_i.ra, df_i.dec, dc_i.ra, dc_i.dec) * 3600.
        dc_i = dc_i.sort_values(by=['dis'])
        dc_i = dc_i.reset_index(drop=True)

        # If target distance is NaN we secure it is placed as first row
        target_row = dc_i[dc_i.gaiaDR3 == df_i.gaiaDR3].index[0]
        df = ut.pdMoveRowToFirst(dc_i, target_row, reset_index=True)
        df.dis.iloc[0] = 0.0
        

        # GENERATE LIGHT CURVES

        # Check contaminant variability
        if conFlag == 'no':
            nstar = 1
        elif conFlag == 'yes':
            nstar = df.shape[0]
        else:
            errorcode('error', 'Not valid mocka.CFLAG value! Use [yes, no]')

            
        # Loop over each star in subfield

        varSourceFiles = []
        for i in range(nstar):

            self.df = df.iloc[i]
            
            
            # FETCH STELLAR PARAMETERS

            # Case 1) Only SpecType
            # Case 2) Only SpecType, BP-RP
            # Case 3) Only SpecType, BP-RP, Teff, logg
            # Case 4) All parameters exist
            
            # Case 1: Draw from spectral type distributions
            if pd.isna(self.df.Ag) and pd.isna(self.df.BP_RP):
                dx   = pd.read_feather(f'{idir}/starcat_GaiaDR3_Teff_SpecType_{df0.spec}.ftr')
                Teff = random.choices(dx.Teff, weights=dx.density, k=1)[0]
                dx   = df0.iloc[ut.findNearestIndex(df0.Teff, Teff)]
                self.df.Pmag = ut.passbandConversionG2P(self.df.Gmag, dx.BP_RP)
                self.df.Teff = dx.Teff
                self.df.logg = dx.logg
                self.df.Z = dx.Z
                self.df.R = dx.R
                self.df.M = dx.M
                self.df.L = dx.L
                
            # Case 2: Draw from nearest match to input catalogue df0
            elif pd.isna(self.df.Teff):
                dx = df0[df0.spec == df0.spec]
                dx = dx.iloc[ut.findNearestIndex(dx.BP_RP, self.df.BP_RP)]
                self.df.Teff = dx.Teff
                self.df.logg = dx.logg
                self.df.Z = dx.Z
                self.df.R = dx.R
                self.df.M = dx.M
                self.df.L = dx.L
                
            # Case 3: Try small parameter space around Teff and logg 
            elif pd.isna(self.df.L):
                try:
                    dx = df0[(self.df.Teff > self.df.Teff_low) &
                             (self.df.Teff < self.df.Teff_upp) &
                             (self.df.logg > self.df.logg_low) &
                             (self.df.logg < self.df.logg_upp)]
                except: dx = df0
                dx = dx.iloc[ut.findNearestIndex(dx.Teff, self.df.Teff)]
                self.df.R = dx.R 
                self.df.M = dx.M
                self.df.L = dx.L
                

            # GENRIC STEPS

            self.stellar_source()
            self.stellar_spectrum()

            
            # SELECT VARIABLE SIGNAL

            # beta Cepheid Bowman et al. 2020 -> see Burssens et al. 2020, Fig. 3
            # Mira stars: Cunha+2020
            
            # Seperate dwarf (MS) and sub-gaint (post MS) stars
            #ds = df[self.R.value < ut.getMainSequenceLimit(df0.Teff)]
            #sg = df[self.R.value > ut.getMainSequenceLimit(df0.Teff)]
            
            # Solar-like oscillator
            
            #self.stellar_granosc()
            #self.stellar_activity()
            #self.stellar_flares()
            #self.star_roap()
            self.star_gdor()
            #self.star_dsct()
            #self.star_ceph()
            
            # if df0.spec == 'O':
            #     self.stellar_gran_osc()
                
            # elif df0.spec == 'B':
            #     self.photometric_standard()

            # elif df0.spec == 'A':
            #     self.photometric_standard()

            # elif df0.spec == 'F':
            #     self.stellar_gran_osc()
                
            # elif df0.spec == 'G':
            #     self.stellar_gran_osc()
            #     self.stellar_activity()

            # elif df0.spec == 'K':
            #     self.stellar_gran_osc()
            #     self.stellar_activity()
            #     self.stellar_flares()

            # elif df0.spec == 'M':
            #     self.stellar_activity()
            #     self.stellar_flares()
                
            # print(df0)
            # exit()
        
            # GENERATE LIGHT CURVE

            # Save each varsource to file
            sfile = 'varsource_'+f'{i+1}'.zfill(3)+'.txt'
            self.ofile = self.odir.joinpath(sfile)
            self.run_prolog()

            # Use cluster name for PLATOnium
            clusterDir = f'$VSC_SCRATCH/platosim/mocka/{starType}/{starDir}/'
            varSourceFiles.append(clusterDir + sfile)
            
        # GENERATE VARIABLE CATALOG FILE
        
        starIDs = np.arange(1, nstar+1).astype(str)
        varSourceList = self.odir / 'varSourceList.txt'
        
        if isinstance(varSourceFiles, str):
            varSourceFiles = [varSourceFiles]
            
        with open(varSourceList, 'w') as f:
            for i in range(nstar):
                f.write(f'{starIDs[i]} {varSourceFiles[i]}\n')



            
#--------------------------------------------------------------#
#                PARSING COMMAND-LINE ARGUMENTS                #
#--------------------------------------------------------------#

parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)


parser.add_argument('-p', '--plot',     action='store_true',     help='Flag to plot the synthetic models')
parser.add_argument('-o', '--ofile',    metavar='STR', type=str, help='Output filename [<path/to/ofile.txt>]')
parser.add_argument('-v', '--verbose',  metavar='INT', type=int, help='Verbosity level [0, 1, 3] (Default: 1)')

obs_group = parser.add_argument_group('OBS PARAMETERS')
obs_group.add_argument('--time',    metavar='DAY',  type=int, help='Duration of simulation (Default: 90 days)')
obs_group.add_argument('--quarter', metavar='NUM',  type=str, help='Quarter number or range of quaters to simulate (Default: 1)')
obs_group.add_argument('--cadence', metavar='SEC',  type=int, help='Cadence of observation (Default: 25 seconds)')
obs_group.add_argument('--inst',    metavar='NAME', type=str, help='Photometric instrument (Default: plato)')
obs_group.add_argument('--seed',    metavar='INT',  type=int, help='Option to bootstrap seed to reproduce results')

star_group = parser.add_argument_group('STAR PARAMETERS')
star_group.add_argument('--star', metavar='NAME', type=str, help='Benchmark star [None, Sun, gDor, dSct, roAp, <Object>]')
star_group.add_argument('--star_params', action='append', type=float, nargs=5, metavar=('M', 'R', 'Teff', 'logg', 'Z'),
                        help='Stellar model parameters [M/Msun, R/Rsun, Teff/K, logg/dex, Z/dex]')
star_group.add_argument('--gran',     metavar='RELATION', type=str, help='Scaling relation of Granulation [Kallinger2014, None]')
star_group.add_argument('--puls',     metavar='RELATION', type=str, help='Scaling relation of Pulsations [Corsaro2013, None]')
star_group.add_argument('--spot',     metavar='BOOL',     type=str, help='Inclusion of stellar spots [True, False] (Default: True)')
star_group.add_argument('--pulslist', metavar='FILE',     type=str, help='Use file with pulsations [periods, amplitudes, phases]')

star_group = parser.add_argument_group('BINARY PARAMETERS')
star_group.add_argument('--binary', metavar='NAME', type=str, help='Benchmark eclipsing binary [None, EB, SMBH, <Object>]')
#star_group.add_argument('--binary_params', action='append', type=float, nargs=5, metavar=('M', 'R', 'Teff', 'logg', 'Z'),
#                        help='Stellar model parameters with units [M/Msun, R/Rsun, Teff/K, logg/rel, Z/rel]')


planet_group = parser.add_argument_group('PLANET PARAMETERS')
planet_group.add_argument('--planet', metavar='NAME', type=str, help='Benchmark planet [None, Earth, hotJupiter, <object>]')
planet_group.add_argument('--planet_params', action='append', type=float, nargs=7, metavar=('t0', 'P', 'e', 'i', 'w', 'Rp', 'Mp'),
                          help='Planet model parameters [t0/days, P/days, i/deg, w/deg, Rp/Rearth, Mp/Mearth]')
#planet_group.add_argument('--phase_curve', action='store_true', help='Flag orbital phase curve (occultation, beaming, ellipsoidal)')
planet_group.add_argument('--ldm',   metavar='MODEL', type=str, help='Limb darkening model [quadratic]')

dis_group = parser.add_argument_group('DISTRIBUTION MODES')
dis_group.add_argument('--kul20', metavar='INT',   type=int, help='Option designed for KUL-TN-20 [0, 1, 2, 3]')
dis_group.add_argument('--mocka', action='append', type=str, nargs=5, metavar=('PROJECT', 'CLASS', 'ID', 'CFLAG', 'ODIR'), help='Option designed for MOCKA')

args = parser.parse_args()

#--------------------------------------------------------------#
#                            WORKFLOW                          #
#--------------------------------------------------------------#

# Monitor script speed
tic = datetime.datetime.now()

# Initialize instance of class
v = VarSim(args)

# Mode for PLATO-CS
if args.mocka:
    v.mode_mocka()

# Mode for KUL-TN-20
elif args.kul20:
    v.mode_kul20()

# Default mode for binaries
elif args.binary:
    v.mode_binary()
    
# Default mode for single stars
else:
    v.mode_single()

# Print run time
if (args.verbose is None) or (args.verbose > 0):
    toc = datetime.datetime.now()
    print(f'\nTotal execution time : {toc-tic} [hh:mm:ss]\n')

    
