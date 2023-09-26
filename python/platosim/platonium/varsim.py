#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
This script is an integrated part of PlatoSim's toolkit PLATOnium.
Given a star and exoplent this script creates a synthetic stellar
and exoplanet variability model that can be used directly as input
for PlatoSim. A bolometric correction to the PLATO passband is 
determined from a best fit interpolation between the passband and
the synthetic high-res PHOENIX spectra. The script can model:

Solar-like stars (F5-K7 dwarf and subgiants):
  - Activity modulations
  - Granulation noise 
  - Convection driven oscillations (p-modes)

Other types of stars:
  - gamma-Doradus star (analytic model)
  - delta-Scuti star   (analytic model)
  - Cepheids/RR Lyrae  (from files)

Exoplanet phase curve variations:
  - Transits               (using BATMAN)
  - Occultations           (using SPIDERMAN)
  - Limb Darkening         (using LDTk)
  - Doppler beaming        (using PyAstronomy)
  - Ellipsoidal distortion (using PyAstronomy)

User examples:
  $ varsim --star Sun --planet hotJupiter -p
  $ varsim --star CoRoT-1 --planet CoRoT-1b --time 10 -p -o </path/to/varsouce.txt>
  $ varsim --star_params 1 1 5800 4.5 0 --planet_params 1 10 0 90 0 1 1 --quarter 1-8 -o </path/to/varsource.txt>
"""

# Built-in
import os
import sys
import glob
import random
import argparse
import datetime

# External packages
import natsort
import numpy as np
import pandas as pd
from pathlib import Path
from matplotlib import pyplot as plt
from scipy.ndimage import median_filter
from scipy.interpolate import make_interp_spline
from astropy.io import fits
from astropy import constants as c
from astropy import units as u

# External for exoplanet
from ldtk import LDPSetCreator, TabulatedFilter
import batman

# PlatoSim
import platosim.plot      as pt
import platosim.utilities as ut
from platosim.utilities import errorcode, downloadFromFTP
from platosim.spectrum  import Spectrum
from platosim.varsource import load_star, load_exoplanet
from platosim.varsource import (LimbDarkening,
                                DopplerBeaming,
                                EllipsoidalDistortion,
                                PlanetMRforecast,
                                StellarActivity,
                                StellarFlares,
                                GravityOscillations,
                                solarosc)


#-------------------------------------------------------#
#                     CONFIGURATION                     #
#-------------------------------------------------------#
        
# Construct a new Generator with the default BitGenerator (PCG64)
from numpy.random import default_rng
rng = default_rng()

# Monitor script speed
tic = datetime.datetime.now()

# Constants
t_quater = 90.     # [d]
t_day    = 86400.  # [s]


#==============================================================#
#                         BEGIN CLASS                          #
#==============================================================#

class VarSim(object):

    """Class to generate noise-less light curves.

    NOTE All input parameters with a physical unit needs a astropy.units attached.
    """
    
    def __init__(self, args):
        
        # VERBOSITY (a.k.a log level) -> Identical to PlatoSim usage
        # verbose = 0: Cluster mode: Disabling print and warnings, and no log files are saved
        # verbose = 1: Default mode: Print details to bash but do not save log files
        # verbose = 3: Debug mode  : Print details to bash and saves all log files
        if args.verbose == 0:
            self.verbose = 0
            # Turn off warnings
            import warnings
            warnings.filterwarnings("ignore")
        elif args.verbose is None or args.verbose == 1:
            self.verbose = 1
            # Turn off warnings
            import warnings
            warnings.filterwarnings("ignore")
        else:
            self.verbose = 3

        # SETTINGS

        # Add latex font if catalogue is saved
        if args.outfile is None:
            from platosim.matplotlibrc import setup
            setup()
        else:
            from platosim.matplotlibrc import latex
            latex()

        # Use LaTeX when plotting
        if args.plot:
            plt.rcParams['text.usetex'] = True
        
        # Absolute cwd path
        self.path     = str(Path(os.getenv("PLATO_PROJECT_HOME"))) + "/inputfiles"
        self.datapath = self.path + '/data_varsim'

        
        # DOWNLOAD DATA
        if not Path(self.datapath + '/response_plato.txt').is_file():
            errorcode('message', 'Inuaguration: Welcome to the Variable Source Simulator!')
            print(f'Downloading a few prerequisite files..')
            ut.downloadFromFTP(filename='response_plato.txt',
                               outputDir=self.datapath, server='plato')
            ut.downloadFromFTP(filename='response_tess.txt',
                               outputDir=self.datapath, server='plato')
            ut.downloadFromFTP(filename='response_kepler.txt',
                               outputDir=self.datapath, server='plato')
            ut.downloadFromFTP(filename='meunier_19a_t1.txt',
                               outputDir=self.datapath, server='plato')
            ut.downloadFromFTP(filename='Main_Fits_BiSON_8640d_lbest_UseInSolarCycle.txt',
                               outputDir=self.datapath, server='plato')

            
        # OBSERVATIONAL PARAMETERS

        # Cadence/sampling [s -> d]
        if args.samp:
            cadence = args.samp / 86400.
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
            timeStart = round(90. * (Q[0]-1) )
            timeDur   = round(90. * numQuarters)
        elif args.time:
            timeStart = 0
            timeDur   = args.time
        else:
            timeStart = 0
            timeDur   = 90
        
        # Time series points
        # NOTE Ensure that the time series has an even number of time points
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
        passband = np.loadtxt(self.datapath + '/response_{0}.txt'.format(self.instrument))
        self.wvl_tele = passband[:,0] * u.nm  # Wavelengths [nm]
        self.tra_tele = passband[:,1]         # Normalized transmission
        
        # Define pandas data frame to store all signals
        self.lc = pd.DataFrame(data=self.time.to('s').value, columns=['time'])
        
        if self.verbose > 0:
            print(f'Simulating time series for : ' +
                  f'{len(self.time)} x {self.cadence.to("s")} ({self.timeDur.value} days)')
            print(f'Simulating {self.instrument} bandpass : ' +
                  f'{self.wvl_tele[0]} - {self.wvl_tele[-1]}')

            
    #--------------------------------------------------------------#
    #                         MODEL OF STAR                        #
    #--------------------------------------------------------------#

    
    def stellar_source(self):

        """Select the stellar paramters.
        """

        # Check star source or use Sun as default
        if args.star:
            self.star_source = args.star
        else:
            self.star_source = 'Sun'
        
        if args.star_params:
            # If star parameters are given
            star = args.star_params[0]
            self.M    = star[0] * u.M_sun
            self.R    = star[1] * u.R_sun
            self.Teff = star[2] * u.K
        else:
            # If star source from selection
            star = load_star(self.star_source)
            self.M    = star[0] 
            self.R    = star[1]
            self.Teff = star[2]
        # Same parameters
        self.logg = star[3]
        self.Z    = star[4]
        
        # Use theoretical logg if user has set logg=0
        self.g = c.G.cgs * self.M.cgs / self.R.cgs**2
        if self.logg == 0:
            self.logg = np.log10(self.g.value)
        
        # Derive theoretical luminosity
        self.Teff_sun = 5777.  # [K]
        self.L = self.R.value**2 * (self.Teff.value/self.Teff_sun)**4 * u.L_sun
 
        # Print available stellar model parameters
        if self.verbose > 0:
            errorcode('module', '\nStellar Spectrum\n')
            print("Spectral Type         : {}".format(self.star_source))
            print("Stellar Mass          : {:.2f}".format(self.M.to('M_sun')))
            print("Photosphere Radius    : {:.2f}".format(self.R.to('R_sun')))
            print("Surface Gravity       : {:.2f}".format(self.g.si))
            print("Log10 Surface Gravity : {:.2f} dex".format(self.logg))
            print("Effective Temperature : {:4.0f}".format(self.Teff.to('K')))
            print("Bolometric Luminosity : {:.2f}\n".format(self.L.to('L_sun')))

        # Return parameters
        self.star_params = [self.M.to('M_sun').value, self.R.to('R_sun').value,
                            self.Teff.value, self.logg, self.Z, 0.0]


            
        
        
    def stellar_spectrum(self):

        """Calculates the bolometric correction from high-res spectra.
        
        This function uses the model grid method described in Sarkar+2018:
        https://academic.oup.com/mnras/article/481/3/2871/5092616
        """

        # Load parameters
        wvl_tele  = self.wvl_tele.to('AA').value
        tran_tele = self.tra_tele
        Teff = self.Teff.value
        logg = self.logg
        Z    = self.Z
        R    = self.R
        L    = self.L
        
        # Initialise class for synthetic spectra
        spec = Spectrum()

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
        
        # Cut off IR part -> above 2 microns and check
        wvl_max  = 20000
        dex      = np.where(wvl1_in > wvl_max)
        wvl1_in  = np.delete(wvl1_in, dex)      # [AA]
        wvl2_in  = np.delete(wvl2_in, dex)      # [AA]
        flux1_in = np.delete(flux1_in, dex)     # [ergs/s/cm2/AA]
        flux2_in = np.delete(flux2_in, dex)     # [ergs/s/cm2/AA]
        
        # Check that the interpolation between the two SEDs can be done
        if len(wvl1_in) != len(wvl2_in):
            errorcode('error', 'Spectra are not of same size! Check interpolation')

        # Create SED for star by interpolating nearby absolutely calibrated spectra
        wvl  = (wvl1_in + wvl2_in) / 2.
        flux = flux1_in + (Teff-Teff_lower) * (flux2_in-flux1_in) * float(Teff_upper-Teff_lower)**-1

        # Measure total luminosity from SED [ergs/s]
        Total_Luminosity1 = 4*np.pi*(R.cgs.value)**2 * np.trapz(flux1_in, wvl1_in)
        Total_Luminosity2 = 4*np.pi*(R.cgs.value)**2 * np.trapz(flux2_in, wvl2_in)

        # Consistency check [ergs/sec]
        # NOTE to compare theo L while using PhoenixAtmos, divide with np.pi
        if self.verbose > 0:
            Total_Luminosity  = 4*np.pi*(R.cgs.value)**2 * np.trapz(flux, wvl) * u.erg/u.s
            print('Theoretical Luminosity [ergs/sec] : {:.3e}'.format(L.to('erg/s')))
            print('Synthetic   Luminosity [ergs/sec] : {:.3e} \n'.format(Total_Luminosity))

        # Measure bolometric luminosity amplitude gradient (vs. temperature variation)
        # This will be used instead of the approximate relation provided by Kjeldsen & Bedding (1995)
        # bol_coeff = round(ut.diff(Total_Luminosity2, Total_Luminosity1) /
        #                   ut.diff(Teff_upper, Teff_lower), 2)

        # Measure passband luminosity amplitude gradient (vs. bolometric luminosity amplitude)
        dex_wavlMin = ut.findNearestIndex(wvl, wvl_tele[0])
        dex_wavlMax = ut.findNearestIndex(wvl, wvl_tele[-1])
        if dex_wavlMax-dex_wavlMin == 1:
            Luminosity1_lambda = (4*np.pi * (R.cgs.value)**2 *
                                  ( wvl1_in[dex_wavlMax] -  wvl1_in[dex_wavlMin]) *
                                  (flux1_in[dex_wavlMax] + flux1_in[dex_wavlMin]) / 2.)
            Luminosity2_lambda = (4*np.pi * (R.cgs.value)**2 *
                                  ( wvl2_in[dex_wavlMax] -  wvl2_in[dex_wavlMin]) *
                                  (flux2_in[dex_wavlMax] + flux2_in[dex_wavlMin]) / 2.)
        else:
            Luminosity1_lambda = (4*np.pi * (R.cgs.value)**2 *
                                  np.trapz(flux1_in[dex_wavlMin:dex_wavlMax],
                                           wvl1_in[dex_wavlMin:dex_wavlMax]))
            Luminosity2_lambda = (4*np.pi * (R.cgs.value)**2 *
                                  np.trapz(flux2_in[dex_wavlMin:dex_wavlMax],
                                           wvl2_in[dex_wavlMin:dex_wavlMax]))

        # Bolometric correction cofficient
        self.bol_coeff = (ut.diff(Luminosity2_lambda, Luminosity1_lambda) /
                          ut.diff(Total_Luminosity2, Total_Luminosity1))
        if self.verbose > 0:
            print(f'Bolometric to {self.instrument} passband correction ' +
                  f'coefficient: {round(self.bol_coeff, 4)}')

        # Rebinned spectrum TODO delete?
        wvl_equi = np.arange(wvl_tele[0], wvl_tele[-1], 1)
        wvl_equi, flux_equi = ut.rebin3(wvl_equi, wvl, flux)

        # PROLOGUE

        # Plot the above calculations
        bb_flux = 1
        if args.plot:
            pt.plot_sed(wvl, wvl1_in, wvl2_in, wvl_equi,
                        flux, flux1_in, flux2_in, flux_equi,
                        Teff, Teff_upper, Teff_lower)




            
    #--------------------------------------------------------------#
    #                   MODELS OF SOLAR-LIKE STARS                 #
    #--------------------------------------------------------------#


    def stellar_gran_osc(self):
        """
        Function to simulate stellar granulation and stochastic oscillation (p-modes).
        Each synthetic component of stellar variability is generated in the time domain
        and from asteroseismic scaling relations, for which a vararity of relations are
        are made available to choose from for the user. This function apply the bolometric
        correction when observing in the PLATO passband to each signal compponent.

        Resources
        ---------
        Kjeldsen & Bedding (1995) : https://arxiv.org/abs/astro-ph/9403015
        Michel et al.      (2005) : https://arxiv.org/abs/0809.1078
        Chaplin et al.     (2009) : https://arxiv.org/abs/0905.1722
        Kjeldsen & Bedding (2011) : https://arxiv.org/abs/1104.1659
        Ballot et al.      (2011) : https://arxiv.org/abs/1105.4557
        Kallinger et al.   (2014) : https://arxiv.org/abs/1408.0817

        Code & Data (p-modes)
        ---------------
        De Ridder et al.   (2006) : https://academic.oup.com/mnras/article/365/2/595/976827
        Broomhall et al.   (2009) : Tables in paper 
        """

        # Start script
        if self.verbose > 0:
            errorcode('module', '\nStochastic Oscillator\n')

        # Default scaling relations
        if args.gran is None:
            args.gran = 'Kallinger2014'
        if args.puls is None:
            args.puls = 'Corsaro2013'

        # We use the solar frequency spectrum as a template for the pulsations
        # But scale the frequencies and amplitudes according to the scaling relations
        # (we do not scale the mode lifetimes)
        Teff_sun       = 5777.  # [K]
        numax_sun      = 3140.  # [muHz]
        deltanu_sun    = 134.9  # [muHz]
        tau_gran_sun   = 375.   # [s]
        tau_puls_sun   = (2.88 * u.d).to('Ms').value
        A_gran_bol_sun = 41.    # [ppm] Kallinger et al. (2014) 
        A_puls_bol_sun = 3.6    # [ppm] Michel et al. (2009)
        
        # Convert units [Ms => microHz in frequency]
        # Load stellar parameters [solar]
        M       = self.M.to('M_sun').value
        R       = self.R.to('R_sun').value
        L       = self.L.to('L_sun').value
        Teff    = self.Teff.to('K').value/self.Teff_sun
        time    = self.time.to('Ms').value
        cadence = self.cadence.to('Ms').value

        # Solar frequency spectrum: mode line-widths and frequencies
        data = np.loadtxt(self.datapath + '/Main_Fits_BiSON_8640d_lbest_UseInSolarCycle.txt')
        eta  = np.exp(data[:,6])*np.pi  # From Chaplin et al. (2009) Eq. 1
        freq = data[:,2]                # frequencies of 96 modes

        # Frequency of maximum power and the primary frequency splitting
        numax   = M / R**2 / np.sqrt(Teff) * numax_sun  # [Sun] Keldsen & Bedding (1995) eq. 10
        deltanu = np.sqrt(M) * R**-1.5 * deltanu_sun    # [Sun] Keldsen & Bedding (1995) eq. 9
        if self.verbose > 0:
            print('Frequency of maximum amplitude : {} microHz'.format(round(numax,2)))
            print('Primary frequency splitting    : {} microHz\n'.format(round(deltanu,2)))

        # Scaling solar distinct pulsation frequencies of the 96 modes of the Sun [microHz]
        freq = (freq-numax_sun)/deltanu_sun * deltanu + numax

        # TODO Some scaling relations are not documented
        #intfreqs = np.linspace(np.min(freq), np.max(freq), 10000)
        #tmin = np.max(1./freq * u.Ms).to('min')
        #tmax = np.min(1./freq * u.Ms).to('min')

        # GRANULATION

        if args.gran == None:
            A_gran_bol = 0.

        elif args.gran == 'KjeldsenBedding2011':  # TODO incomplete model?

            # Timescale from Keldsen & Bedding (2011) Eq. 9 [Sun]
            tau_gran = L * M**-1 * Teff**-3.5

            # Amplitude from Kallinger et al. 2014 table 3 TODO where is this equation from?
            A_gran_bol = L * tau_gran**0.5 * M**-1.5 * Teff**-2.75 * 41.

        elif args.gran == 'Kallinger2014':
            """
            The idea to construct the PSD from Kallinger et al. (2014) Eq. 2 using only the
            granulation component modelled by the contribution of 1-3 super-Lorentzian
            functions. Here we only use 2 super-Lorentzian functions.
            """

            # Bolometric correction that scales with the Kepler bandpass: Sec. 5.4
            # Originally from Ballot et al. (2011); Michel et al. (2005)
            T0    = 5934.  # [K]
            alpha = 0.8
            C_bol = (Teff*Teff_sun/T0)**alpha

            # Parameters of power law fits to granulation: Sec. 5.2, from Tab. 2:
            # a**2 correspond to the area under the super-Loretzian in PSD, hence,
            # the variance in the time series, and bolmetric correction account for the bandpass.
            # {b1, b2} are the characteristic frequencies defined by the chracteristic timescale tau
            a  = C_bol * 3710. * numax**-0.613 * M**-0.26  # [ppm]
            b1 = 0.317 * numax**0.970                      # [microHz]
            b2 = 0.948 * numax**0.992                      # [microHz]

            # Prepare PSD model
            Nfreq = int(len(time)/2 + 1)
            frequencies = np.arange(float(Nfreq)) / (Nfreq-1) * 0.5 / cadence

            # Add Gaussian noise TODO add new numpy rng!
            realPart = np.random.normal(0., .5, Nfreq)
            imagPart = np.random.normal(0., .5, Nfreq)

            # Construct PSD [ppm^2/muHz] and the full fourier spectrum
            psd1 = ut.superLorentzian(frequencies, b1, a)
            psd2 = ut.superLorentzian(frequencies, b2, a)
            psd  = psd1 + psd2
            fourier = np.sqrt(2. * psd * (Nfreq-1) / cadence) * (realPart + imagPart * 1j)
            signal_gran = np.real(np.fft.irfft(fourier))

        # OSCILLATIONS

        if args.puls == None:
            A_puls_bol = 0.

        elif args.puls == 'KB1995Brown1991':
            # According to Corsaro et al. (2013) Eq. 6 [ppm]
            A_puls_bol = (numax/numax_sun)**-1 * Teff**1.5 * A_puls_bol_sun

        elif args.puls == 'KjeldsenBedding1995':
            # According to Kjeldsen and Bedding (1995), Eq. 4 usinf Eq. 8 [ppm]
            A_puls_bol = L * Teff**-1 * M**-1 * 4.7 * 550/623

        elif args.puls == 'Huber2011': # TODO
            # According to the relation by Huber et al. 2011b [ppm]
            s = 0.886; t = 1.89; r = 2.0
            A_puls_bol = L**s * M**-t * Teff**(1-r) * A_puls_bol_sun
            # Re-estimated by Corsaro et al. 2013 using Bayesian inference
            #s = 0.984; t = 1.66; r = 2.79
            # According to the relation by Huber et al. 2011b,
            # and the fit by Corsaro et al. 2013, equ. (19), [ppm]
            #ampl_puls_bol = ((numax/numax_sun)**(2*s-3*t) * (deltanu0/deltanu_sun)**(-4*s+4*t) *
            #                (Teff[0]/5777.)**(5*s-1.5*t-r+0.2) * 3.6)

        elif args.puls == 'Mosser2010': # TODO
            # According to Corsaro et al. 2013, equ. (24)
            tau_puls = conversions.convert('d', 'Ms', np.exp((5777. - Teff[0])/601.) * 2.65)
            r = 1.5  # free parameter
            # According to Kjeldsen and Bedding 2011 equ. (6), [ppm]
            A_puls_bol = (L[0] * (tau_puls/conversions.convert('d','Ms',2.65 ))**0.5 *
                          M[0]**-1.5 * (Teff[0]/5777.)**-(1.25+r) * 3.6)

        elif args.puls == 'KjeldsenBedding2011':  # TODO
            # According to Corsaro et al. 2013, equ. (24)
            tau_puls = conversions.convert('d', 'Ms', np.exp((5777. - Teff[0])/601.) * 2.65)
            r = 2.0
            # According to Kjeldsen and Bedding 2011 equ. (6), [ppm]
            A_puls_bol = (L[0] * (tau_puls/conversions.convert('d','Ms',2.65))**0.5 * M[0]**-1.5 *
                          (Teff[0]/5777.)**-(1.25+r) * 3.6)

        elif args.puls == 'Corsaro2013':

            # According to Corsaro et al. 2013, Eq. 24:
            # NOTE The value of T0 was calibrated using Kepler RGs in the open clusters
            # NGC 6791 and NGC 6819, and a sample of MS and subgiant Kepler field stars.
            T0   = 601.                         # [K]
            tau0 = (2.65 * u.d).to('Ms').value  # [Ms]
            tau_puls = np.exp((Teff_sun - Teff*Teff_sun) / T0) * tau0  # [Ms]

            # According to Corsaro et al. (2013) Eq. 26 (Model 6) [ppm]
            r = -2.8
            t = 1.56
            A_puls_bol = ( (numax/numax_sun)**(2-3*t) * (deltanu/deltanu_sun)**(4*t-4) *
                           (tau_puls/tau_puls_sun)**0.5 * Teff**(4.55-r-1.5*t) * A_puls_bol_sun )

        # Amplitude [ppm]
        A = A_puls_bol * np.exp(-(freq-numax)**2/(2*(1.5*deltanu)**2))

        if self.verbose > 0:
            #print('Bolometric granulation amplitude : {} ppm'.format(round(A_puls_bol,2)))
            print('Bolometric pulsation amplitude : {} ppm\n'.format(round(A_puls_bol,2)))

        # Timeseries of oscillations [ppm]
        signal_puls = solarosc(time, freq, A, eta, self.verbose)

        # PROLOGUE

        # Fix that the granulation time series is some time 1 time point shorter
        if len(signal_gran) != len(signal_puls):
            signal_gran = np.append(signal_gran, np.mean(signal_gran))

        # Combine model and correct for bandpass
        lc_gran = self.bol_coeff * signal_gran
        lc_puls = self.bol_coeff * signal_puls

        # print(np.sqrt(np.mean(signal_gran**2)))
        # print(np.sqrt(np.mean(signal_puls**2)))
        # print(np.sqrt(np.mean(lc_gran**2)))
        # print(np.sqrt(np.mean(lc_puls**2)))
        
        self.lc['gran'] = lc_gran.tolist()
        self.lc['puls'] = lc_puls.tolist()
        self.puls_params = [numax, deltanu, self.bol_coeff]

        # Plot rsults (time [Ms] -> freq [muHz])
        if args.plot:
            lc_tot = lc_gran + lc_puls  
            pt.plot_amplitude_time_series(time, lc_gran, lc_puls, lc_tot, self.star_source)





            
    def stellar_activity(self):

        """Model stellar spot modulations.

        This function simulate a synthetic noise-less light curve of main-sequence
        stars that include stellar activity in the form of cyclic spot modulations.

        Resources
        ---------
        Noyes et al.            (1984) : https://adsabs.harvard.edu/pdf/1984ApJ...279..763N
        Pillet et al.           (1993) : https://www.cambridge.org/core/journals/international-
                                         astronomical-union-colloquium/article/distribution-of-
                                         sunspot-decay-rates/9D2174592A1C0CD0E9DD8B1DF347B7FF#
        Baumann and Solanki     (2005) : https://www.aanda.org/articles/aa/abs/2005/45/aa3415-05/aa3415-05.html
        Mamajek and Hillenbrand (2008) : https://iopscience.iop.org/article/10.1086/591785/meta
        Llama et al.            (2012) : https://academic.oup.com/mnrasl/article/422/1/L72/971190?login=true
        Aigrain et al.          (2012) : https://academic.oup.com/mnras/article/419/4/3147/2908053?login=true
        Meunier et al.          (2019) : https://arxiv.org/abs/1911.05319
        
        Assumptions
        -----------
        - Model of spots only (missing e.g. faculae)
        - Model only valid for main-sequence stars
        
        Code courtesy
        -------------
        Suzanne Aigrain : Aigrain et al. (2012)
        """

        from platosim.starspot import simulate_lc
        
        # Start script
        if self.verbose > 0:
            errorcode('module', '\nStellar activity\n')

        # Use a random uniform distribution
        # NOTE secure a lower misalignment for planetary systems
        if args.planet or args.planet_params or args.xsource:
            incl = np.random.uniform(85, 90)
        else:
            # None means cos(i_star) uniform between 0-90 deg
            incl = None
            
        # Re-run random activity cycle until it fits in memory!
        lc = None
        while lc is None:
            lc, params = simulate_lc(teff=self.Teff.value,
                                     time=self.time.to('d').value,
                                     dur=self.timeDur.to('d').value,
                                     cadence_hours=self.cadence.to('h').value,
                                     incl=incl,
                                     verbose=self.verbose,
                                     doplot=args.plot)

        # Store global variables
        lc = lc * 1e6
        self.lc['spot'] = lc.tolist()
        self.spot_params = params








    def stellar_flare(self, tscale=False, tmax=False, amplitude=30, asymmetry=1):

        """Function to model flares.
        
        Parameters
        ----------
        tscale : float, ndarray
            Time scale duration of the flare(s) [days]
        tmax : float, ndarray
            Full time-width at half-maximum-flux of the flare(s) maximum intensity [days]
        amplitude : float, ndarray
            Amplitude of the flare(s) [mmag]
        asymmetry : float, ndarray
            Asymmetry factor of the flare(s)
        """

        # Convert units of input parameters
        time = self.time.to('d').value
        flux = np.zeros_like(time)

        # Secure that single flare works
        try: len(tmax)
        except: tmax = [tmax]
        
        # Loop over each flare event
        
        for m in range(len(tmax)):

            # Start and end of flare event
            t0 = (time[0]  - tmax[m])
            t1 = (time[-1] - tmax[m])
            dt = np.diff(time)[0]

            # Time array during flare event
            tn = np.arange(t0, t1, dt)
            t = tn/tscale

            # Model parameters of flare
            B = asymmetry
            C = 1/B
            b = -1.941 - 0.175 + 2.246 + 1
            c = 1 - 0.689

            # Loop over every time-step in the flare time interval
            # NOTE: this is defined relative to this flares maxima and put in units
            # of the time-scale here the analytic expressions for the rise and decay
            # are used to determine the flux of this flare

            for i in range(len(t)):

                # Rise of flare
                if t[i]*B > -1 and t[i]*B <= 0:
                    flux[i] += (1
                                + 1.941 * (t[i]*B)
                                - 0.175 * (t[i]*B)**2
                                - 2.246 * (t[i]*B)**3
                                - b     * (t[i]*B)**4)

                # Decay of flare
                elif t[i]*C > 0:
                    flux[i] += 0.689 * np.exp(-1.6 * t[i]*C) + c * np.exp(-0.2783 * t[i]*C)

                # No flare
                else:
                    flux[i] += 0

        # Convert to magnitude [mmag]
        mag = flux * amplitude
                    
        # plot light curve
        # if args.plot:
        #     plt.figure(figsize=(10, 5))
        #     plt.plot(time, mag , 'm-')
        #     plt.xlabel('Time [d]')
        #     plt.ylabel(r'$\delta m$ [mmag]')
        #     plt.xlim(np.min(time), np.max(time))
        #     plt.tight_layout()
        #     plt.show()
            

        return flux * amplitude + 1

        


            
    #--------------------------------------------------------------#
    #                     NON SOLAR-LIKE STARS                     #
    #--------------------------------------------------------------#
    

    def photometric_standard(self):

        """Function to model rotationally modulated chemcically perculiar A-type (Ap) stars.

        roAp stars are common and readily available objects that efficiently can
        be used as photometric calibrators in larger surveys. Rotational periods
        can be drawn randomly from a uniform distribution between 1 and 3 days. 
        The shape of light curves of these stars is typically sinusoidal with 
        frequent contribution of the second harmonic. The amplitude is typically
        10-30 mmag in the Kepler passband. We here assume that Kepler passband is
        representative for the PLATO passband.

        Author: Oleg Kochukhov (oleg.kochukhov@physics.uu.se)
        """
        # Start script
        if self.verbose > 0:
            errorcode('module', '\nPhotometric standard\n')

        # Convert units of input parameters
        time = self.time.to('d').value

        # Random value of rotational period in [1, 3] d range
        P = np.random.uniform(low=1, high=3)

        # Relative amplitudes of the main frequency and harmonic
        a1 = 1.
        a2 = np.random.uniform(low=0, high=a1)

        # Random phase offset [0, 2*pi] between main frequency and harmonic
        dphi = np.random.uniform(low=0, high=2*np.pi)

        # Create light curve and scale amplitude in [10,30] mmag range
        mag   = a1 * np.cos(2*np.pi * (time/P)) + a2 * np.cos(4*np.pi * (time/P) + dphi)
        scale = np.random.uniform(low=10, high=30) / max(abs(mag))
        mag  *= scale

        # plot light curve
        if args.plot:
            plt.figure(figsize=(10,3.5))
            plt.plot(time, mag, 'b-')
            plt.xlabel('Time [d]')
            plt.ylabel('Magnitude [mmag]')
            plt.tight_layout()
            plt.show()

        # Return [delta mag]
        mag *= 1e-3 
        self.lc['mag'] = mag.tolist()
        self.std_params = [P, dphi*180/np.pi, a1, a2, scale]



        


    def gravity_oscill(self, period_range, amplitude_range, nmodes, power=2.2, seed=0):

        """Function to generate g-mode pulsations

        Parameters
        ----------
        time_start : 
        """
        
        # time_start and time_end are given in days and define the time interval over which to simulate the flare
        # the sampling should be given in exposures or data-points per day
        # the period and amplitude define the allowed range in periods and amplitudes to include
        # the number_modes defines how many different modes or periods/amplitudes to include
        # power is the power to which the sum of every mode is raised, it introduces an asymmetry in the signal

        # Convert units of input parameters
        time = self.time.to('d').value
        mag  = np.zeros_like(time)

        # Check if a file with pulsations are parsed
        if args.pulslist:
            f = Path(args.pulslist).resolve()
            # Catch if the file doesn't exist
            if not f.is_file():
                errorcode('error', 'File do not exist, check filepath again!')
            else:
                data = np.loadtxt(f)
            # Else load the data
            period    = data[:,0]
            amplitude = data[:,1]
            phase     = data[:,2]
            nmodes = len(period)
        else:
            # Randomly generate pulsations if not from file
            period    = np.zeros(nmodes)
            amplitude = np.zeros(nmodes)
            phase     = np.zeros(nmodes)
            for i in range(nmodes):
                phase[i]     = random.uniform(0, 2*np.pi)
                period[i]    = random.uniform(period_range[0], period_range[1])
                amplitude[i] = random.uniform((amplitude_range[0]), amplitude_range[1])

        # Loop over the number of modes
        # Flux is the sum of every mode
        for i in range(nmodes):
            mag += amplitude[i] * np.sin(2 * np.pi * (1 / period[i]) * time + phase[i])

        # Normalize the flux so its values lie in [-1, 1] (so roots are not undefined)
        # Then add 1, raise the power and substract 1
        A = np.amax(np.absolute(mag))
        mag /= A
        mag = A * ( ( (1 + mag)**power ) - 1)

        # Create a table with value
        df = pd.DataFrame()
        df['period [days]']    = period 
        df['amplitude [mmag]'] = amplitude
        df['phase [rad]']      = phase
        
        # plot light curve
        if args.plot:
            print(df)
            plt.figure(figsize=(10, 5))
            plt.plot(time, mag, 'm-')
            plt.xlabel('Time [d]')
            plt.ylabel(r'$\delta m$ [mmag]')
            plt.xlim(np.min(time), np.max(time))
            plt.tight_layout()
            plt.show()
        
        return mag * 1e-3 




    
    def gamma_doradus(self):

        """Function to generate pulsations of a gamma-Dor star.

        Notes 
        -----
        This function used the "gravity_oscill" utility with a characteristic power
        of 2.2 for these g-mode pulsators.
        """

        # Start script
        if self.verbose > 0:
            errorcode('module', '\ngamma-Dor g-modes\n')

        # Draw number of pulsation from uniform distribution
        npuls = int(rng.normal(50, 5))

        # Model gravity modes
        flux = self.gravity_oscill([0.8, 3], [0.5, 2.5], nmodes=npuls, power=2.2, seed=0)
        self.lc['mag'] = flux.tolist()
        #self.std_params = [P, dphi*180/np.pi, a1, a2, scale]


    
        

    def delta_scuti(self):

        """Function to generate pulsations of delta scuti stars.

        Notes 
        -----
        This function used the "gravity_oscill" utility with a characteristic power
        of 1.0 for these g-mode pulsators.
        """

        # Start script
        if self.verbose > 0:
            errorcode('module', 'delta-Scuti g-modes\n')

        # Draw number of pulsation from uniform distribution
        npuls = int(rng.normal(50, 5))

        # Model pulsations
        flux = self.gravity_oscill([1,30], [10,30], nmodes=npuls, power=1.0, seed=0)
        self.lc['mag'] = flux.tolist()
        #self.std_params = [P, dphi*180/np.pi, a1, a2, scale]





    def classical_pulsator(self):

        """Function to generate pulsations of Cepheids and RR Lyrae stars.

        This function uses precomputed models of Cepheids and RR Lyrae
        stars to generate the light curve from their harmonics. For now
        this function will randomly select a star.
        """

        # Start script
        if self.verbose > 0:
            errorcode('module', '\nClassical pulsator\n')

        # Select a random object from the list and load Fourier data
        try:
            filenames = glob.glob(f'{self.datapath}/RRL-CEP/*.fou')
            starfile = random.choice(filenames)
        except:
            zipfile = 'RRL-CEP.zip'
            errorcode('message', 'Classic, I like your style!')
            print(f'Downloading {zipfile} files..')
            downloadFromFTP(filename=zipfile, outputDir=self.datapath, server='plato')
            os.system(f'unzip {self.datapath}/{zipfile} -d {self.datapath}')
            os.system(f'rm {self.datapath}/{zipfile}')
            print('')

        # Load file with harmonics
        filenames = glob.glob(f'{self.datapath}/RRL-CEP/*.fou')
        starfile  = random.choice(filenames)
        fourier   = np.loadtxt(starfile)    
        if self.verbose > 0:
            print(f'Using file {starfile} with frequencies:')
            print(fourier)
            
        # Convert units of input parameters
        time = self.time.value
        flux = np.zeros_like(time)
        lc = pd.DataFrame(data = {'time':time, 'flux':flux})

        # Generate the lightcurve in the dataframe
        components = len(np.array(fourier))
        for i in range(components):
            lc.flux = (lc.flux + fourier[i,1] *
                       np.sin((2*np.pi*fourier[i,0] * time) + fourier[i,2]))

        # Save magnitude to list
        self.lc['mag'] = - 2.5 * np.log10(lc.flux + 1)

        # plot light curve
        if args.plot:
            plt.figure(figsize=(10, 5))
            plt.plot(time, self.lc.mag*1e3, 'm-')
            plt.xlabel('Time [d]')
            plt.ylabel(r'$\delta m$ [mmag]')
            plt.xlim(np.min(time), np.max(time))
            plt.tight_layout()
            plt.show()




        
    def smbh_binary(self):

        """Function to generate a SMBH binary light curve.

        A Super Massive Black Hole (SMBH) binary system consist of several components,
        for which this model includes two of the effects:
        - The doppler boosting
        - The gravitational lensing effect

        """

        # Start script
        if self.verbose > 0:
            errorcode('module', '\nSMBH binary\n')

        # Fetch time array
        time = self.time.to('d').value

        # Fetch model parameters
        self.distribution = 'random'
        if self.distribution == 'random':
            from platosim.distribution import SMBHB
            smbhb = SMBHB()
            P, phi, Abeam, Aflare, tscale = smbhb.randomToyModel() 
        
        # Model doppler beaming effect
        #A   = 0.04          # [mag]
        #P   = 2 * 365.25  # [days]
        #phi = 0.1 * np.pi         # [rad] #np.random.uniform(low=0, high=2*np.pi) 
        flux_beam = Abeam * np.sin(2*np.pi * (time/P) + phi) + 1

        # Model the flare event
        #tscale    = 10     # [days] #np.ones(len(max)) * 10
        #amplitude = 0.08   # [mag]

        tmax = P * (1 - phi/(2*np.pi))
        flux_flare = np.ones_like(flux_beam)

        for tm in [tmax-P, tmax, tmax+P]:
            
            flux_flare += self.stellar_flare(tscale, tm, Aflare, asymmetry=1) - np.ones_like(flux_beam)

        if self.verbose:
            print(f'Model parameters of toy model:')
            print(f'Orbital period     : {P} day')
            print(f'Beaming amplitude  : {Abeam} mag')
            print(f'Flare time scale   : {tscale} day')
            print(f'Flare time maximum : {tmax} day')
            print(f'Flare amplitude    : {Aflare} mag')        
        # Combine model
        flux = flux_beam + flux_flare - 1
        mag  = -2.5 * np.log10(flux)
        
        # plot light curve
        if args.plot:
            fig, ax = plt.subplots(2, 1, figsize=(10,8))
            ax[0].plot(time, flux_beam,  '-', c='green')
            ax[0].plot(time, flux_flare, '-', c='orange')
            ax[0].plot(time, flux,       '-', c='royalblue')
            ax[0].set_xlabel('Time [d]')
            ax[0].set_ylabel(r'Relative flux')
            ax[0].set_xlim(time.min(), time.max())
            ax[1].plot(time, mag*1e3, '-', c='royalblue')
            ax[1].set_xlabel('Time [d]')
            ax[1].set_ylabel(r'$\delta m$ [mmag]')
            ax[1].set_xlim(time.min(), time.max())            
            plt.tight_layout()
            plt.show()

        # Return 
        self.lc['mag'] = mag.tolist()
        #self.std_params = [P, dphi*180/np.pi, a1, a2, scale]



        
        
    #--------------------------------------------------------------#
    #                       EXOPLAENTS TRANSITS                    #
    #--------------------------------------------------------------#
        
        
    def ldc(self):
        """
        To compute our custom limb-darkening transit duration coefficients that meet
        the PLATO transmission response function. We used the angle-dependent ("MU")
        Specific Intensity Spectra (SIS) from PHOENIX (Goettingen 2018), which exactly
        as above is a library of the stellar effective temperature, surface gravity,
        and metallicity. The limb darkening are naturally calclated for the exact same
        stellar parameter as used for the granulation and oscialltions. See links:
        
        RESOURCES
        ---------
        Webpage : https://phoenix.astro.physik.uni-goettingen.de/
        Download: http://phoenix.astro.physik.uni-goettingen.de/data/
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
            errorcode('warning', f"LD coefficients failed for (Teff, logg, Z) = ({self.Teff}, {self.logg}, {self.Z}")
            self.ldc = [0.430, 0.170]
        else:
            self.ldc = u[0]
            

        

            
    def exoplanet_model(self):
        """
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

        # Convert units of input parameters
        time = self.time.to('d')
        Ms   = self.M.to('kg')
        Rs   = self.R.to('m')
        Teff = self.Teff.to('K')
        
        # Extract bandpass
        wvl_tele = self.wvl_tele
        tra_tele = self.tra_tele

        # LOAD PLANET MODEL
        
        if args.planet == 'random':

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
            Mp, _, _ = mr_forecast.Rstat2M(mean=Rp.to('R_jup').value, std=0.01, unit='Jupiter',
                                           sample_size=1000, grid_size=1000, classify=classifier)
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
                params = load_exoplanet(args.planet)
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

        # DYNAMICS

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

        # First to fourth contact [d]
        # t1_tra = t0_tra_cen - t_tra_tot/2.
        # t2_tra = t0_tra_cen - t_tra_ful/2.
        # t3_tra = t0_tra_cen + t_tra_ful/2.
        # t4_tra = t0_tra_cen + t_tra_tot/2.

        if self.verbose > 0:
            errorcode('module', '\nPlanet Eclipses')
            print('')
            print("Planet name                  : {}".format(args.planet))
            print("Planet mass                  : {:.3f}".format(Mp.to('M_earth')))
            print("Planet radius                : {:.3f}".format(Rp.to('R_earth')))
            print('')
            print("Semimajor axis               : {:.3f} starRad".format(a.to('m')/Rs.to('m')))
            print("Eccentricity                 : {:.3f}".format(e))
            print("Inclination                  : {:.3f}".format(i.to('deg')))
            print("Argument of Periastron       : {:.1f}".format(w.to('deg')))
            print('')
            print("Orbital Period               : {:.3f}".format(P.to('d')))
            print("Transit-to-Occultation time  : {:.3f}".format(dt_c.to('d')))
            print("Time of emphemeris           : {:.3f}".format(t0.to('d')))
            print('')
            print("Total transit duration       : {:.3f}".format(t_tra_tot.to('h')))
            print("Full  transit duration       : {:.3f}".format(t_tra_ful.to('h')))
            print("In/Egress tra duration       : {:.3f}".format(tau_tra.to('min')))
            print('')
            print("Total occultation duration   : {:.3f}".format(t_occ_tot.to('h')))
            print("Full  occultation duration   : {:.3f}".format(t_occ_ful.to('h')))
            print("In/Egress occult. duration   : {:.3f}".format(tau_occ.to('min')))
            print('')
            print("Impact parameter Transit     : {:.3f}".format(b_tra))
            print("Impact parameter Occultation : {:.3f}".format(b_occ))
            print('')
            print("Limb darkening coefficients  : {:.3f}, {:.3f}".format(self.ldc[0], self.ldc[1]))

        # Check for model applicability
        if b_tra > 1 - Rp/Rs and b_tra <= 1 + Rp/Rs:
            errorcode('warning', 'Planetary model consist of grazing eclipses!')
        elif b_tra > 1 + Rp/Rs:
            errorcode('error', 'Planet model do not have any physical eclipses!')

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

            



    def exoplanet_transit(self):
        """
        Function to model exoplanet transits.

        In the following the exoplanet transits are being modelled with Batman:
        https://lweb.cfa.harvard.edu/~lkreidberg/batman/quickstart.html

        NOTE All times t0 and P can principly be anything as long as they are
        consistant. Here we make sure to use consistent reference time unit.
        """

        # Limb darkening model options:
        if args.ldm: limbDarkModel = args.lmd
        else: limbDarkModel = 'quadratic'

        # Initialize batman model
        batman_params = batman.TransitParams()
        batman_params.t0  = self.t0.to('d').value
        batman_params.per = self.P.to('d').value
        batman_params.a   = (self.a.to('m')/self.R.to('m')).value
        batman_params.ecc = self.e
        batman_params.inc = self.i.to('deg').value
        batman_params.w   = self.w.to('deg').value
        batman_params.rp  = (self.Rp.to('m')/self.R.to('m')).value
        batman_params.u   = self.ldc
        batman_params.limb_dark = limbDarkModel

        # Initializes transit model and extract light curve [ppm]
        model_tra = batman.TransitModel(batman_params, self.time.value)
        lc_tra    = (model_tra.light_curve(batman_params) - 1) * 1e6

        # True anomaly at each time: This will be used in our custom models later
        self.nu = model_tra.get_true_anomaly()

        # Time of periastron passage (calculated from t0)
        self.tau = model_tra.get_t_periastron(batman_params)

        # Return transit model
        self.lc['tran'] = lc_tra.tolist()
        self.exo_params = [self.Mp.to('M_earth').value, self.Rp.to('R_earth').value,
                           (self.a.to('R_sun')/self.R).value, self.P.to('d').value,
                           self.t0.to('d').value, self.e, self.i.to('deg').value,
                           self.w.to('deg').value, self.ldc[0], self.ldc[1]]






        

    def exoplanet_occultation(self):
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





        
    def exoplanet_beaming(self):
        """
        Doppler beaming model.
        """

        # Central wavelength of PLATO bandpass [m] TODO can this be done better?
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







    def exoplanet_ellipsoidal(self):
        """
        ELLIPSOIDAL DISTORTION
        """

        # Initialize and prepare model input:
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
                                        self.t0.to('d').value, self.P.to('d').value,
                                        self.dt_c.to('d').value,
                                        self.t0_tra_cen.to('d').value, self.t_tra_tot.to('d').value,
                                        self.t0_occ_cen.to('d').value, self.t_occ_tot.to('d').value,
                                        self.A_beam, self.A_elli)
        elif (self.time[-1] < self.P.to('d') + self.t0.to('d')):
            errorcode('warning', 'No phase plot, time series is shorter than the orbital period!')






    #--------------------------------------------------------------#
    #                      PROLOGUE AND SAVING                     #
    #--------------------------------------------------------------#
        

    def run_prolog(self):
        
        if self.verbose > 0:
            errorcode('module', '\nPrologue')
            
        # Compute delta magnitude of signal
        variable = ('std', 'dSct', 'gDor', 'Cep', 'SMBHB')
        if args.star in variable:
            dm = self.lc['mag']
        else:
            # Combine all signals
            
            # Granulation and pulsation are additive
            self.lc['comb'] = np.zeros(len(self.lc.time))
            if 'gran' in self.lc:
                self.lc['comb'] += self.lc.gran
            if 'puls' in self.lc:
                self.lc['comb'] += self.lc.puls
            if 'spot' in self.lc:
                self.lc['comb'] += self.lc.spot

            # Convert to relative flux to multiply with transits
            self.lc['comb'] = self.lc['comb'] / 1e6 + 1 
                
            # Spots and transits are multiplicative
            if 'tran' in self.lc:
                self.lc['comb'] *= (self.lc.tran / 1e6 + 1)
                
            # Convert to delta magnitude
            dF = self.lc['comb'].to_numpy()
            dm = - 2.5 * np.log10(dF)

            # Convert back again for plot (normalized to 0 [ppm])
            self.lc.comb = (self.lc.comb - 1) * 1e6
            
        # Convert to seconds
        if args.quarter:
            self.lc['time'] += self.timeStart * 86400

        # PLOT FINAL LIGHT CURVE

        if args.plot and args.star not in variable:
            pt.plot_final_lc(self.lc)
            plt.show()

        # OUTPUT PANDAS TABLE
        
        # Collect NaN info for standard stars
        try: self.star_params
        except: self.star_params = np.zeros(8) * np.nan
        try: self.spot_params
        except: self.spot_params = np.zeros(10) * np.nan
        try: self.puls_params
        except: self.puls_params = np.zeros(4) * np.nan
        try: self.exo_params
        except: self.exo_params = np.zeros(11) * np.nan
        try: self.std_params
        except: self.std_params = np.zeros(6) * np.nan
        else: self.spot_params[2] = self.std_params[0]

        # Create output table
        d = {'Ms_Msun': [self.star_params[0]],
             'Rs_Rsun': [self.star_params[1]],
             'Teff_K': [self.star_params[2]],
             'logg': [self.star_params[3]],
             'Z': [self.star_params[4]],
             'alpha': [self.star_params[5]],
             
             'BV': [self.spot_params[0]],
             'logRHK': [self.spot_params[1]],
             'Prot_day': [self.spot_params[2]],
             'Pmin_day': [self.spot_params[3]],
             'Pmax_day': [self.spot_params[4]],
             'lmax_deg': [self.spot_params[5]],
             'Pcyc_year': [self.spot_params[6]],
             'Povl_year': [self.spot_params[7]],
             'Acyc': [self.spot_params[8]],
             'is_deg': [self.spot_params[9]],

             'numax_muHz': [self.puls_params[0]],
             'deltanu_muHz': [self.puls_params[1]],
             'bol_coeff': [self.puls_params[2]],
             
             'dphi_deg': [self.std_params[1]],
             'a1': [self.std_params[2]],
             'a2': [self.std_params[3]],
             'scale': [self.std_params[4]],
             
             'Mp_Mearth': [self.exo_params[0]],
             'Rp_Rearth': [self.exo_params[1]],
             'a_Rstar': [self.exo_params[2]],
             'P_day': [self.exo_params[3]],
             't0_day': [self.exo_params[4]],
             'e': [self.exo_params[5]],
             'ip_deg': [self.exo_params[6]],
             'w_deg': [self.exo_params[7]],
             'u1': [self.exo_params[8]],
             'u2': [self.exo_params[9]]
        }
        df = pd.DataFrame(d)

        # Print table
        if self.verbose > 0 and args.star not in variable:
            print('\nUsed parameters space:')
            print(df.T) 

        # SAVE DATA
        
        if args.outfile:

            print(f'\nSaving output files')
            out = args.outfile[-3:]

            # Save to ascii
            if out == 'txt':
                np.savetxt(args.outfile, np.transpose([self.lc['time'], dm]), fmt=['%.1f', '%.8f'])
                if not self.star_source in ['SMBHB']:
                    df.to_feather(f'{args.outfile[:-4]}_parameters.ftr')
                    self.lc.to_feather(f'{args.outfile[:-4]}_components.ftr')

            # Save to numpy binary
            elif out == 'npy':
                np.save(args.outfile, np.transpose([self.time.to('s').value, dm]))

            # Save feather binary
            elif out == 'ftr':
                df = pd.DataFrame({'time': self.time.to('s').value, 'dmag': dm},
                                  columns=['time', 'dmag'])
                df.to_feather(args.outfile)
                df.to_feather(f'{args.outfile[:-4]}_params.{out}')
                    
            else:
                errorcode('Error', 'Output format not supported!')

        # Print for the simulation statistics
        if self.verbose > 0:
            toc = datetime.datetime.now()
            print('\nComputation time : {0} [hh:mm:ss]\n'.format(toc-tic))



            
#--------------------------------------------------------------#
#                PARSING COMMAND-LINE ARGUMENTS                #
#--------------------------------------------------------------#

software = '\nVariable Source Simulator\n'
parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', software))

parser.add_argument('-p', '--plot',    action='store_true', help='Flag to plot the synthetic models')
parser.add_argument('-v', '--verbose', metavar='NUM',  type=int, help='Verbosity level [0, 1, 2] (Default: 1)')
parser.add_argument('-o', '--outfile', metavar='NAME', type=str, help='Filename of output file (Use extension: ".txt"')
parser.add_argument('-x', '--xsource', metavar='INT',  type=int, help='Specific flag for KUL-TN-20 simulations [0, 1, 2, 3]')

obs_group = parser.add_argument_group('OBSERVATION')
obs_group.add_argument('--time',    metavar='DAYS', type=int, help='Time duration of simulation [days] (Default: 30 days)')
obs_group.add_argument('--quarter', metavar='NUM',  type=str, help='Quarter number or range of quaters to simulate (Default: False)')
obs_group.add_argument('--samp',    metavar='SECS', type=int, help='Time cadence of observation [seconds] (Default: 25 sec)')
obs_group.add_argument('--inst',    metavar='NAME', type=str, help='Observational instrument (Default: "PLATO")')

star_group = parser.add_argument_group('STAR')
star_group.add_argument('--star', metavar='NAME', type=str, help='Stellar variability source [<Object>, roAp, gDor, dSct] (Default: None)')
star_group.add_argument('--star_params', action='append', type=float, nargs=5, metavar=('M', 'R', 'Teff', 'logg', 'Z'),
                        help='Stellar model parameters [M: Msun, R: Rsun, Teff: K] (Default: None)')
star_group.add_argument('--gran',     metavar='RELATION', type=str, help='Scaling relation of Granulation [Kallinger2014, None] (Default: Kallinger2014)')
star_group.add_argument('--puls',     metavar='RELATION', type=str, help='Scaling relation of Pulsations  [Corsaro2013,   None] (Default: Corsaro2013)')
star_group.add_argument('--corr',     metavar='METHOD',   type=str, help='Scaling Correction method of p-modes [None] (Default: None)')
star_group.add_argument('--spot',     metavar='BOOL',     type=str, help='Inclusion of stellar spots    [True, False] (Default: True)')
star_group.add_argument('--pulslist', metavar='FILE',     type=str, help='Use custum list of pulsations {periods, amplitudes, phases}')

planet_group = parser.add_argument_group('EXOPLANET')
planet_group.add_argument('--planet', metavar='NAME',  type=str, help='Exoplanet variability source (Default: None)')
planet_group.add_argument('--planet_params', action='append', type=float, nargs=7, metavar=('t0', 'P', 'e', 'i', 'w', 'Rp', 'Mp'),
                          help='Planet model parameters [t0: days, P: days, i: deg, w: deg, Rp: Rearth, Mp: Mearth] (Default: None)')
planet_group.add_argument('--ldm',    metavar='MODEL', type=str, help='Limb Darkening model (Default: quadratic)')
planet_group.add_argument('--phase_curve', action='store_true',  help='Flag to include orbital phase curve model {Occultation, Beaming, Ellipsoidal} (Default: False)')

args = parser.parse_args()

#--------------------------------------------------------------#
#                            WORKFLOW                          #
#--------------------------------------------------------------#

# xsource: used for KUL20 stiching- and detrending
# 0 -> Std/constant (2 hamonics)
# 1 -> Gran, Puls
# 2 -> Gran, puls, Spot
# 3 -> Gran, Puls, Spots, Exo
# x -> Constant stars is any other number x
if args.xsource == 0:
    args.star = 'roAp'
if args.xsource in (1, 2, 3):
    args.star          = 'Sun'
    args.planet_params = False
if args.xsource == 2:
    args.spot          = True
    args.planet        = False
if args.xsource == 3:
    args.spot          = True
    args.planet        = 'random'
if not args.xsource in (0, 1, 2, 3):
    args.xsource = False

# Activate spot modulation by default
if args.spot is True or args.spot is None:
    args.spot = True

# Initialize instance of class
v = VarSim(args)

# Include stellar variability
if args.star == 'roAp':
    v.photometric_standard()
elif args.star == 'gDor':
    v.gamma_doradus()
elif args.star == 'dSct':
    v.delta_scuti()
elif args.star == 'Cep':
    v.classical_pulsator()
elif args.star == 'EB':
    v.eclipsing_binary()
elif args.star == 'SMBHB':
    v.smbh_binary()
else:
    # Select star
    v.stellar_source()
    v.stellar_spectrum()
    # Solar-like stars
    if args.star or args.star_params:
        if args.spot is True:
            v.stellar_activity()
        if not args.gran or not args.puls:
            v.stellar_gran_osc()    

    # Include exoplanet
    if args.planet or args.planet_params or args.planet == 'random':
        v.ldc()
        v.exoplanet_model()
        v.exoplanet_transit()
        v.exoplanet_occultation()
        v.exoplanet_beaming()
        v.exoplanet_ellipsoidal()
        if not args.xsource and args.plot:
            v.plot_phase_curve()

# Combine and save
v.run_prolog()


# if __name__ == "__main__":
#     main()

