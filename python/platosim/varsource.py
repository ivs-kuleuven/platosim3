#!/usr/bin/env python

"""
This script contains all relevant functions and classes that 
are used by the PLATOnium script "varsim.py".

NOTE This class needs the Poetry install: 
     >> poetry install --with platonium 
"""

# Built-in
import os
import glob
import math
import random
import urllib.request
from pathlib import Path

# PlatoSim standard
import h5py
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from scipy.stats import norm, truncnorm
from scipy.interpolate import interp1d, make_interp_spline
from astropy.io import fits
from astropy.table import Table
from astropy import units as u
from astropy import constants as c

# PLATOnium extra
from numba import njit
from PyAstronomy import funcFit, pyasl

# PlatoSim functions
import platosim.plot      as pt
import platosim.noise     as ns
import platosim.utilities as ut
from platosim.utilities import errorcode


#==============================================================#
#                        SOLAR-LIKE STARS                      #
#==============================================================#
    
    
class StellarFlares(object):

    """Model stellar flares.

    A simplistic analytical description of stellar flares described by
    an sudden flux increase followed by an exponential decay. Given the
    time of the time series the corresponding flux is returned including
    the wanted flares.
    """

    def __init__(self, time, seed=False):
        
        # Store array
        self.time = time
        self.rng  = ut.rng(seed)



    def initToyModelBeta0(self): # TODO under construction!

        """Uniform distribution of toy model.
        """

        # Time in [days] and ampl in [mmag]
        nflares = self.rng.integers(0, 10, 1)[0]
        self.tscale = self.rng.uniform(0.01, 0.1, nflares)
        self.tmax   = self.rng.uniform(0, self.time[-1], nflares)
        self.ampl   = self.rng.uniform(0, 0.2, nflares)

        return self.tscale, self.tmax, self.ampl

    
        
    def evaluate(self, tscale=False, tmax=False, ampl=False, asym=1, plot=False):

        """Analytic model of stellar flares.
        
        Parameters
        ----------
        tscale : float, ndarray
            Time scale duration of the flare(s) [days]
        tmax : float, ndarray
            Full time-width at half-maximum-flux of the flare(s) maximum intensity [days]
        ampl : float, ndarray
            Amplitude of the flare(s) [mmag]
        asym : float, ndarray
            Asymmetry factor of the flare(s)
        """

        # Check parsing
        if not tscale: tscale = self.tscale
        if not tmax:   tmax   = self.tmax
        if not ampl:   ampl   = self.ampl
        
        # Placeholders
        flux     = np.zeros_like(self.time)
        self.mag = np.zeros_like(self.time)

        # Sampling
        dt = np.diff(self.time)[0]
        
        # Secure that single flare works
        try: len(tmax)
        except: tmax = [tmax]
        
        # Loop over each flare event
        
        for m in range(len(tmax)):

            # Start and end of flare event
            t0 = (self.time[0]  - tmax[m])
            t1 = (self.time[-1] - tmax[m])

            # Time array during flare event
            tn = np.arange(t0, t1, dt)
            t  = tn / tscale[m]

            # Model parameters of flare
            B = asym
            C = 1/B
            b = -1.941 - 0.175 + 2.246 + 1
            c = 1 - 0.689

            # Loop over every time-step in the flare time interval:
            # NOTE: this is defined relative to this flares maxima
            # and put in units of the time-scale here the analytic
            # expressions for the rise and decay are used to determine
            # the flux of this flare

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
                    flux[i] += (0.689 * np.exp(-1.6    * t[i]*C) +
                                c     * np.exp(-0.2783 * t[i]*C))

                # No flare
                else:
                    flux[i] += 0

            # Convert to magnitude [mmag]
            self.mag = flux * ampl[m]
                    
        # plot light curve
        if plot: self.plot()
            
        # Return relative flux
        return self.mag



    def plot(self):

        """Function to plot result.
        """
        plt.figure(figsize=(9, 5))
        plt.plot(self.time, self.mag, 'k-')
        plt.xlabel('Time [d]')
        plt.ylabel(r'$\delta m$ [mmag]')
        plt.xlim(np.min(self.time), np.max(self.time))
        plt.tight_layout()
        plt.show()



            
    
class StellarSpots(object):
    
    """Class to generate rotational star spot modulations.

    This function simulate a synthetic noise-less light curve of main-sequence
    stars that include stellar activity in the form of cyclic spot modulations.

    Resources
    ---------
    Noyes et al.            (1984) : https://adsabs.harvard.edu/pdf/1984ApJ...279..763N
    Pillet et al.           (1993) : https://www.cambridge.org/core/journals/international-
                                     astronomical-union-colloquium/article/distribution-of-
                                     sunspot-decay-rates/9D2174592A1C0CD0E9DD8B1DF347B7FF#
    Baumann and Solanki     (2005) : https://www.aanda.org/articles/aa/abs/2005/45/aa3415-
                                     05/aa3415-05.html
    Mamajek and Hillenbrand (2008) : https://iopscience.iop.org/article/10.1086/591785/meta
    Llama et al.            (2012) : https://academic.oup.com/mnrasl/article/422/1/L72/
                                     971190?login=true
    Aigrain et al.          (2012) : https://academic.oup.com/mnras/article/419/4/3147/
                                     2908053?login=true
    Meunier et al.          (2019) : https://arxiv.org/abs/1911.05319

    Assumptions
    -----------
    - Model of spots only (missing e.g. faculae)
    - Model only valid for main-sequence stars

    Code courtesy
    -------------
    Suzanne Aigrain : Aigrain et al. (2012)
    """

    def __init__(self, seed=False):
        
        # Random number generator
        self.rng = ut.rng(seed)

        # Constants
        self.BV_SUN    = 0.656
        self.LRHK_SUN  = -5.025 # from Lorenzo-Oliveira et al. (2018, A&A 619, A73)
        self.PROT_SUN  = 27.0
        self.OMEGA_SUN = 2 * np.pi / (self.PROT_SUN * 86400)

        # Load table
        idir = os.getenv('PLATO_PROJECT_HOME') + '/inputfiles/data_varsim'
        self.t1 = Table.read(f'{idir}/varsim_meunier19a_t1.txt', format = 'ascii')


        
    ####################################
    # FROM TEFF TO ACTIVITY PARAMETERS #
    ####################################
    # All relations used are based on Meunier et al. (2019, A&A, 627, A56)
    # except where specified otherwise

    def get_stpar_from_teff(self, teff):

        # Check if Teff is within grid
        if teff < min(self.t1['Teff']) or teff > max(self.t1['Teff']):
            errorcode('Warning', 'Teff is outside range of conversion table. ' +
                      f'Returning B-V for Teff={min(self.t1["Teff"])}')
            if teff < min(self.t1['Teff']): return max(self.t1['BV'])
            if teff > max(self.t1['Teff']): return min(self.t1['BV'])

        g = interp1d(self.t1['Teff'], self.t1['BV'])
        return g(teff)

    
    def get_lrhk_from_S_and_bv(self, S, bv):
        # Cf Noyes et al. (1984, ApJ 279 763, Appendix a)
        lCcf = 1.13 * bv**3 - 3.91 * bv**2 + 2.84 * bv - 0.47 
        if bv < 0.63:
            x = 0.63 - bv
            lCcf += 0.135 * x - 0.814 * x**2 + 6.03 * x**3
        return -4 + np.log10(1.34) + lCcf + np.log10(S)


    def get_lrhk_from_bv(self, bv):
        if bv < 0.94:
            Smin = 0.144
        else:
            Smin = 0.0269231 * bv + 0.118892
        lrhkmin = self.get_lrhk_from_S_and_bv(Smin, bv)
        lrhkmax = -0.375 * bv - 4.4
        return self.rng.random() * (lrhkmax - lrhkmin) + lrhkmin

    
    def get_ltauc_from_bv(self, bv):
        # cf Noyes et al. (1984, ApJ 279 763, Eqn 4)
        x = 1.0 - bv
        if x > 0:
            return 1.362 - 0.166 * x + 0.025 * x**2 - 5.323 * x**3
        else:
            return 1.362 - 0.14 * x

        
    def get_prot_from_lrhk_and_bv(self, lrhk, bv):
        Ro = 0.808 - 2.966 * (lrhk + 4.52)
        delta = self.rng.random() * 0.4 - 0.2
        ltc = self.get_ltauc_from_bv(bv)
        return (Ro + delta) * 10**ltc 

    
    def get_prange_from_teff_and_prot(self, teff, prot):
        p0 = -3.485 + 2.47810e-4 * teff
        p1 = 1.597 - 1.3510e-4 * teff
        alpha = 10**(p0 + p1 * np.log10(prot))
        pmax = 2 * prot / (2 - alpha)
        pmin = pmax * (1 - alpha)
        return pmin, pmax

    
    def get_latrange(self): 
        lat_min = 0.0
        lat_max = 32.0 + 20.0 * self.rng.random()
        return lat_min, lat_max # in degrees

    
    def get_omega01_from_prange_and_latrange(self, pmin, pmax, lat_min, lat_max):
        omega_min = 2 * np.pi / pmax / 86400
        omega_max = 2 * np.pi / pmin / 86400
        s2min = np.sin(np.deg2rad(lat_min))**2
        s2max = np.sin(np.deg2rad(lat_max))**2
        omega_1 = (omega_min - omega_max) / (s2max - s2min)
        omega_0 = omega_max - omega_1 * s2min 
        return omega_0, omega_1 

    
    def get_omega_from_lat_and_omega01(self, lat, omega_0, omega_1):
        # In radians per second
        return omega_0 + omega_1 * np.sin(np.deg2rad(lat))**2


    def get_pcyc_from_prot(self, prot):
        delta = self.rng.random() * 0.6 - 0.3
        y = 0.84 * np.log10(1/prot) + 3.14 + delta
        return prot * 10**y

    
    def get_acyc_from_bv_and_lrhk(self, bv, lrhk):
        if bv < 0.851:
            Acyc_max = 0.727 * bv - 0.292
        else:
            Acyc_max = 0.727 * 0.851 - 0.292
        Acyc_min = max([0.28 * bv - 0.196, 0.342 * lrhk + 1.703, 0.005])
        return self.rng.random() * (Acyc_max - Acyc_min) + Acyc_min

    
    def get_arate_from_acyc(self, acyc):
        asun = self.get_acyc_from_bv_and_lrhk(self.BV_SUN, self.LRHK_SUN)
        return acyc/asun



    def regions(self, activity_rate=1, cycle_period=10, cycle_overlap=0, randspots=False,
                maxlat=70, minlat=0, tsim=1000, tstart=0, verbose=False):

        """ACTIVE REGION EMERGENCE

        According to Schrijver and Harvey (1994), the number of active regions
        emerging with areas in the range [A,A+dA] in a time dt is given by:

        n(A,t) dA dt = a(t) A^(-2) dA dt ,

        where A is the "initial" area of a bipole in square degrees, and t is
        the time in days; a(t) varies from 1.23 at cycle minimum to 10 at cycle
        maximum.

        The bipole area is the area within the 25-Gauss contour in the
        "initial" state, i.e. time of maximum development of the active region.
        The assumed peak flux density in the initial sate is 1100 G, and
        width = 0.2*bsiz (see disp_region). The parameters written onto the
        file are corrected for further diffusion and correspond to the time
        when width = 4 deg, the smallest width that can be resolved with lmax=63.

        In our simulation we use a lower value of a(t) to account for "correlated"
        regions.
        """
        
        nbin = 5                              # number of area bins
        delt = 0.5                            # delta ln(A)
        amax = 100.                           # orig. area of largest bipoles (deg^2)
        dcon = np.exp(0.5*delt)-np.exp(-0.5*delt)   # contant from integ. over bin
        latrmsd = 6
        atm     = 10 * activity_rate    
        # a(t) at cycle maximum (deg^2/day)
        # cycle period (days)
        # cycle duration (days)

        ncycle = int(cycle_period * 365)     # cycle length in days   
        nclen = int((cycle_period + cycle_overlap) * 365)
        fact = np.exp(delt*np.arange(nbin))  # array of area reduction factors
        ftot = fact.sum()                    # sum of reduction factors
        bsiz = np.sqrt(amax/fact)            # array of bipole separations (deg)
        tau1 = 5                             # first and last times (in days) for
        tau2 = 15                            #   emergence of "correlated" regions
        prob = 0.001                         # total probability for "correlation"
        nlon = 36                            # number of longitude bins
        nlat = 16                            # number of latitude bins       
        nday1 = 0                            # first day to be simulated
        ndays = int(tsim)                    # number of days to be simulated
        dt = 1

        # Initialize time since last emergence of a large region, as function
        # of longitude, latitude and hemisphere:
        tau = np.zeros((nlon,nlat,2),'int') + tau2
        dlon = 360. / nlon
        dlat = maxlat / nlat

        # Create arrays to store regions properties
        reg_tims = []
        reg_lats = []
        reg_lons = []
        reg_angs = []

        # Loop over time (in days):
        ncnt = 0
        ncur = 0
        start_day = 0

        for nd in range(ndays):
            nday = nd + nday1

            # Compute index of most recently started cycle:
            ncur_now = int(nday / ncycle)
            ncur_prev = int((nday-1) / ncycle)
            if ncur_now > ncur_prev:
                ncur = ncur + 1

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

                nstart = start_day        # start date of cycle
                if (nday-nstart) < nclen:  
                    ic = 1 - 2 * ((nc + 2) % 2) # +1 for even, -1 for odd cycle
                    phase = float(nday-nstart) / nclen # phase within the cycle

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
                        latrms = (maxlat/5.) - latrmsd * phase # rms latitude (degrees)
                        nlat1 = np.floor(max([maxlat * 0.9 - 1.2 * maxlat * phase, 0.0]) /
                                         dlat).astype(int) # first and last index
                        nlat2 = np.floor(min([maxlat + 15. - maxlat * phase, maxlat]) /
                                         dlat).astype(int)
                        nlat2 = min([nlat2, nlat - 1])

                    js = np.arange(nlat2 - nlat1).astype(int)

                    p = np.zeros(nlat)
                    for j in np.arange(nlat2-nlat1+1).astype(int) + nlat1:
                        p[j] = np.exp( - ((dlat * (0.5 + j) - latavg) / latrms)**2)
                    ru0 = ru0_tot * p / (p.sum() * nlon * 2)

                    # Loops over hemisphere and latitude:
                    for k in [0,1]:
                        for j in np.arange(nlat2-nlat1+1).astype(int) + nlat1:
                            # Emergence rates of largest regions per
                            # longitude/latitude bin (number per day):
                            r0 = ru0[j] + rc0[:,j,k]
                            rtot = r0.sum()
                            ssum = rtot * ftot
                            x = self.rng.random()
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
                                lon = dlon * (self.rng.random() + float(i))
                                lat = dlat * (self.rng.random() + float(j))
                                if (nday > tstart):
                                    reg_tims.append(self.rng.random() + nday)
                                    reg_lons.append(lon)
                                    if k == 0:                  # Insert on N hemisphere
                                        reg_lats.append(lat)
                                    else:
                                        reg_lats.append(-lat)
                                    x = self.rng.normal()
                                    while abs(x) > 1.6:
                                        x = self.rng.normal()
                                    y = self.rng.normal()
                                    while abs(y) >= 1.6:
                                        y = self.rng.normal()
                                    z = self.rng.random()
                                    if z > 0.14:
                                        # Tilt angle (degrees)
                                        ang = 0.5 * lat + 2.0 + 27. * x * y
                                    else:
                                        z = self.rng.normal()
                                        while z > 0.5:
                                            z = self.rng.normal()
                                        ang = np.deg2rad(z)
                                    reg_angs.append(ang)
                                    if verbose:
                                        print(reg_tims[-1], reg_lats[-1],
                                              reg_lons[-1], reg_angs[-1])
                                ncnt = ncnt + 1
                                if nb < 1:
                                    tau[i,j,k] = 0

        if verbose:
            print('Total number of regions: ', ncnt)

        reg_arr = np.zeros((4, len(reg_tims)))
        reg_arr[0] = np.array(reg_tims)
        reg_arr[1] = np.array(reg_lats)
        reg_arr[2] = np.array(reg_lons)
        reg_arr[3] = np.deg2rad(np.array(reg_angs))
        return reg_arr
    


    def spots(self, reg_arr, incl=None, omega_0=None, omega_1=0.0,
              dur=None, threshold=0.1):

        """Holds parameters for spots for a given star.
        
        Generate initial parameter set for spots (emergence times
        and initial locations.
        """

        # Set global parameters which are the same for all spots inclination [deg]
        if incl == None:
            self.incl = np.rad2deg(np.arcos(np.random.uniform()))
        else:
            self.incl = incl

        if omega_0 is None:
            omega_0 = self.OMEGA_SUN
            
        # Rotation and differential rotation [rad/s]
        self.omega_0 = omega_0
        self.omega_1 = omega_1
        
        # Regions parameters
        t0 = reg_arr[0,:]
        lat = reg_arr[1,:]
        lon = reg_arr[2,:]
        ang = reg_arr[3,:]

        # Keep only spots emerging within specified time-span with:
        # peak B-field > threshold
        if dur == None:
            self.dur = t0.max() 
        else:
            self.dur = dur
        l = (t0 < self.dur) * (ang > threshold)
        self.nspot = l.sum()
        self.t0 = t0[l]
        self.lat = lat[l]
        self.lon = lon[l]
        
        # The settings below are designed approximately match the distributions
        # used in Borgniet et al. (2015) and Meunier et al. (2019) spot sizes
        self.amax = ang[l]**2 * 300 * 1e-6
        
        # spot emergence and decay timescales
        mea = 15 * 1e-6
        med = 10 * 1e-6
        mu = np.log(med)
        sig = np.sqrt(2*np.log(mea/med))
        self.decay_rate = self.rng.lognormal(mean=mu, sigma=sig, size=self.nspot)


    def calci(self, time, i):

        """Single spot calculation.

        Evolve one spot and calculate its impact on the stellar flux.
        NOTE: Currently there is no spot drift or shear.
        """
        
        # Spot area (linear growth and decay)
        area        = np.zeros(len(time)) 
        decay_time  = self.amax[i] / self.decay_rate[i]
        emerge_time = decay_time / 10.0
        
        # exponential growth and decay
        l = time < self.t0[i]
        area[l] = self.amax[i] * np.exp(-(self.t0[i]-time[l]) / emerge_time)
        l = time >= self.t0[i]
        area[l] = self.amax[i] * np.exp(-(time[l]-self.t0[i]) / decay_time)
        
        # linear growth and decay
        # l = (time >= (self.t0[i]-emerge_time)) * (time < self.t0[i])
        # area[l] = self.amax[i] * (self.t0[i]-time[l]) / emerge_time
        # l = (time >= self.t0[i]) * (time < (self.t0[i]+decay_time))
        # area[l] = self.amax[i] * (1-(time[l]-self.t0[i]) / decay_time)

        # Rotation rate [rad/s]
        ome = self.get_omega_from_lat_and_omega01(self.lat[i], self.omega_0, self.omega_1)
        
        # Fore-shortening 
        phase = ome * time * 86400 + np.deg2rad(self.lon[i]) # [rad]
        beta  = (np.cos(np.deg2rad(self.incl)) * np.sin(np.deg2rad(self.lat[i])) +
                 np.sin(np.deg2rad(self.incl)) * np.cos(np.deg2rad(self.lat[i])) *
                 np.cos(phase))
        
        # Differential effect on stellar flux
        dF = - 2 * area * beta
        dF[beta < 0] = 0
        return area, ome, beta, dF


    def calc(self, time):

        """Calculations for all spots
        """
        
        N = len(time)
        M = self.nspot

        # Don't go beyond 5 Gb RAM memory!
        bytemax = int(5/8*1e9)
        X = N*M
        if X > bytemax:
            errorcode('warning', f'Too many spots -> Array of size larger than 5 Gb!')
            return None, None, None, None
        else:
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



    def evaluate(self, teff, time, dur, cadence_hours, incl=None, isim=0,
                 odir=None, verbose=False, save=False):

        """Generate spot modulated light curve.
        """
        
        if odir is None:
            odir = os.getcwd()
            
        # select parameters
        bv   = self.get_stpar_from_teff(teff)
        lrhk = self.get_lrhk_from_bv(bv)
        prot = self.get_prot_from_lrhk_and_bv(lrhk, bv)
        pmin, pmax       = self.get_prange_from_teff_and_prot(teff, prot)
        lmin, lmax       = self.get_latrange()
        omega_0, omega_1 = self.get_omega01_from_prange_and_latrange(pmin, pmax, lmin, lmax)
        pcyc     =  self.get_pcyc_from_prot(prot)
        clen     = pcyc / 365.
        coverlap = self.rng.random() * 0.1 * clen 
        acyc     = self.get_acyc_from_bv_and_lrhk(bv, lrhk)
        arate    = self.get_arate_from_acyc(acyc)
        if incl is None:
            incl = np.rad2deg(np.arccos(self.rng.random()))

        # simulate regions
        reg_arr = self.regions(activity_rate=arate,
                               cycle_period=clen, cycle_overlap=coverlap, verbose=False,
                               maxlat=lmax, minlat=lmin, tsim=dur+pcyc, tstart=0)

        # make the simulation start at a random point in the cycle
        reg_arr[0] -= self.rng.random() * pcyc

        # simulate LC
        self.spots(reg_arr, incl=incl, omega_0=omega_0, omega_1=omega_1,
                   threshold=0.1, dur=dur)
        
        # NOTE we decrease the sampling to 30 min for increased performance
        # NOTE this corresponds to every 72nd time point -> 30 * 60 / 25
        time0 = np.copy(time)
        time = time[::72]
        area, ome, beta, dF = self.calc(time)
        
        # Stop here if the RAM memory would have been overflown
        if dF is None:
            return None, None

        # Save data and figure
        if save:

            # save individual spot properties
            header = '{:6s} {:6s} {:6s} {:6s} {:8s} {:6s} {:6s}'.format('LAT','LON', 'PROT', 'T_MAX', 'A_MAX', 'TAU', 'TAU_R')
            flo.write('# {}\n'.format(header))
            header = '{:6s} {:6s} {:6s} {:6s} {:8s} {:6s} {:6s}'.format('deg','deg', 'days', 'days', 'muHem', 'days', 'periods')
            flo.write('# {}\n'.format(header))
            
            for i in range(self.nspot):
                # spot came too early or too late or was too short lived given cadence
                if area[i,:].max() == 0:
                    continue
                prot = 2*np.pi / ome[i] / 86400
                lifetime = self.amax[i] / self.decay_rate[i]
                str_ = '{:6.1f} {:6.2f} {:6.2f} {:6.2f} {:8.2e} {:6.2f} {:6.2f}'.format(self.lat[i], self.lon[i], prot, self.t0[i], self.amax[i] * 1e6, lifetime, lifetime/prot)
                flo.write('{}\n'.format(str_))

            # save LC
            X = np.zeros((2,len(time)))
            X[0,:] = time
            X[1,:] = dF.sum(0)
            lfile = os.path.join(odir, 'lightcurve_{:04d}.txt'.format(isim)) # save LC
            np.savetxt(lfile,X.T)

        # We interpolate back to original time grid of 25s cadence
        # NOTE Interpolate (piecewise cubic) into higher resolution grid
        #time_int = np.linspace(time0[0], time0[-1], len(time0)
        spline = make_interp_spline(time, dF.sum(0), k=3)
        flux   = spline(time0)

        # Finito!
        self.dF, self.dur, self.area, self.time = dF, dur, area, time
        self.params = [bv.tolist(), lrhk, arate, prot, pmin, pmax, clen, coverlap, lmax, incl]
        return flux, self.params

    

    def plot(self):

        """Plot spot modulation model.
        """

        fig, axes = plt.subplots(3,1, figsize=(11,8), sharex=True)
        ttl= ('Model: ' +
              r'${\rm AR} = $'+f'{self.params[2]:5.3f}, ' +
              r'${\rm CL} = $'+f'{self.params[6]:6.3f} yr, ' +
              r'$P_{\rm min} = $'+f'{self.params[4]:6.2f} d, ' +
              r'$P_{\rm max} = $'+f'{self.params[5]:6.2f} d, ' +
              r'$L_{\rm max} = $'+f'{self.params[8]:5.2f}'+r'$^{\circ}$, ' +
              r'$i = $'+f'{self.params[9]:6.2f}'+r'$^{\circ}$')
        axes[0].set_title(ttl, fontsize='18')
        axes[0].set_facecolor('yellow')
        axes[0].axhline(y=0, color='k', linestyle='--')
        for j in range(self.nspot):
            if self.t0[j] < -10:
                continue
            if self.t0[j] > self.dur:
                continue
            axes[0].plot(self.t0[j], self.lat[j], 'ko', alpha=0.8,
                         markersize=self.amax[j]*(1./3e-4)*5)
        axes[0].set_ylim(-90,90)
        axes[0].set_ylabel('Spot latitude [deg]')
        axes[1].plot(self.time, self.area.sum(0)*100, 'k-')
        axes[1].set_ylabel(r'Spot coverage [\%]')        
        axes[2].plot(self.time, self.dF.sum(0)*1e6, 'k-')
        axes[2].set_ylabel('Relative flux [ppm]')
        axes[2].set_xlim(0, self.dur)
        axes[2].set_xlabel('Time [days]')
        plt.tight_layout(h_pad=0.1)
        plt.show()
        
    

    

class SolarLikeOscillator(object):

    """Class to generate gravity oscillation time series.

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
    ---------------------
    De Ridder et al.   (2006) : https://academic.oup.com/mnras/article/365/2/595/976827
    Broomhall et al.   (2009) : Tables in paper
    """

    def __init__(self, time, star_params, path, seed=False):
        
        # Convert units [Ms => microHz in frequency]
        self.time    = time.to('Ms').value
        self.cadence = np.diff(time)[0].to('Ms').value
        self.Teff_sun = 5777.
        
        # Load stellar parameters [solar]
        self.Teff = star_params[0].to('K').value
        self.R    = star_params[1].to('R_sun').value
        self.M    = star_params[2].to('M_sun').value
        self.L    = star_params[3].to('L_sun').value
        self.T    = self.Teff / self.Teff_sun
        
        # Random number generator
        self.rng = ut.rng(seed)
                    
        # Solar frequency spectrum: mode line-widths and frequencies:
        # Frequencies of 96 modes from BiSON
        data = np.loadtxt(f'{path}/varsim_mainFitsBiSON.txt')
        freq = data[:,2]

        # Frequency of maximum power and the primary frequency splitting
        # Keldsen & Bedding (1995) eq. 10 and 9:
        numax_sun    = 3140.                                             # [muHz]
        deltanu_sun  = 134.9                                             # [muHz]
        self.numax   = self.M / self.R**2 / np.sqrt(self.T) * numax_sun  # [Sun] 
        self.deltanu = np.sqrt(self.M) * self.R**-1.5 * deltanu_sun      # [Sun]
            
        # Scaling solar distinct pulsation frequencies of the 96 modes of the Sun [microHz]
        self.freq = (freq-numax_sun)/deltanu_sun * self.deltanu + self.numax
        
        # From Chaplin et al. (2009) Eq. 1
        self.eta = np.exp(data[:,6])*np.pi


        
    def init_granulation(self, scaling='Kallinger2014'):

        """Initialize granulation model.
        """
        
        tau_gran_sun = 375.   # [s]

        # SELECT SCALING RELATION
        
        if scaling == 'KjeldsenBedding2011':  # TODO incomplete model?

            # Timescale from Keldsen & Bedding (2011) Eq. 9 [Sun]
            tau_gran = self.L * self.M**-1 * self.T**-3.5

            # Amplitude from Kallinger et al. 2014 table 3 TODO where is this equation from?
            A_gran_bol = self.L * tau_gran**0.5 * self.M**-1.5 * self.T**-2.75 * 41.

        elif scaling == 'Kallinger2014':
            # The idea to construct the PSD from Kallinger et al. (2014) Eq. 2 using only the
            # granulation component modelled by the contribution of 1-3 super-Lorentzian
            # functions. Here we only use 2 super-Lorentzian functions.

            #A_gran_bol_sun = 41.    # [ppm] Kallinger et al. (2014) 
            
            # Bolometric correction that scales with the Kepler bandpass: Sec. 5.4
            # Originally from Ballot et al. (2011); Michel et al. (2005)
            T0    = 5934.  # [K]
            C_bol = (self.Teff/T0)**0.8

            # Parameters of power law fits to granulation: Sec. 5.2, from Tab. 2:
            # a**2 correspond to the area under the super-Loretzian in PSD, hence,
            # the variance in the time series, and bolmetric correction account for the bandpass.
            # {b1, b2} are the characteristic frequencies defined by the chracteristic timescale tau
            self.a  = C_bol * 3710. * self.numax**-0.613 * self.M**-0.26  # [ppm]
            self.b1 = 0.317 * self.numax**0.970                           # [microHz]
            self.b2 = 0.948 * self.numax**0.992                           # [microHz]

        else:
            errorcode('error', 'Invalid scaling relation!')

        # Return model parameters
        return self.a, self.b1, self.b2
            


    def eval_granulation(self, a=False, b1=False, b2=False):

        """Model stochastic oscillations.
        """

        # Fetch input parameters
        if not a:  a  = self.a
        if not b1: b1 = self.b1
        if not b2: b2 = self.b2
        
        # Prepare PSD model
        Nfreq = int(len(self.time)/2 + 1)
        frequencies = np.arange(float(Nfreq)) / (Nfreq-1) * 0.5 / self.cadence

        # Add Gaussian noise TODO add new numpy rng!
        realPart = np.random.normal(0., 0.5, Nfreq)
        imagPart = np.random.normal(0., 0.5, Nfreq)

        # Construct PSD [ppm^2/muHz] and the full fourier spectrum
        psd1 = ut.superLorentzian(frequencies, b1, a)
        psd2 = ut.superLorentzian(frequencies, b2, a)
        psd  = psd1 + psd2

        fourier = np.sqrt(2. * psd * (Nfreq-1) / self.cadence) * (realPart + imagPart * 1j)
        self.signal_gran = np.real(np.fft.irfft(fourier))

        # Return generated signal
        return self.signal_gran



    def init_oscillations(self, scaling='Corsaro2013'):

        """Model stochastic oscillations.
        """
        
        # We use the solar frequency spectrum as a template for the pulsations but
        # scale the frequencies and amplitudes according to the scaling relations
        # (we do not scale the mode lifetimes) -> Michel et al. (2009)
        numax_sun      = 3140.                         # [muHz]
        deltanu_sun    = 134.9                         # [muHz]
        tau_puls_sun   = (2.88 * u.d).to('Ms').value   # [Ms]
        A_puls_bol_sun = 3.6                           # [ppm]

        # Handy definitions
        delta_Teff   = self.Teff_sun - self.Teff     # [K]
        numax_frac   = self.numax / numax_sun        
        deltanu_frac = self.deltanu / deltanu_sun
        
        # SELECT SCALING RELATION
        
        if scaling == 'KB1995Brown1991':
            # According to Corsaro et al. (2013) Eq. 6 [ppm]
            A_puls_bol = numax_frac**-1 * self.T**1.5 * A_puls_bol_sun

        elif scaling == 'KjeldsenBedding1995':
            # According to Kjeldsen and Bedding (1995), Eq. 4 using Eq. 8 [ppm]
            A_puls_bol = self.L * self.T**-1 * self.M**-1 * 4.7 * 550/623

        elif scaling == 'Mosser2010': # TODO
            r = 1.5  # free parameter
            tau0 = convert('d','Ms', 2.65)
            # According to Corsaro et al. 2013, equ. (24)
            tau_puls = convert('d', 'Ms', np.exp(delta_Teff / 601.) * 2.65)
            # According to Kjeldsen and Bedding 2011 equ. (6), [ppm]
            A_puls_bol = self.L * (tau_puls/tau0)**0.5 * self.M**-1.5 * (self.T)**-(1.25+r) *3.6
            
        elif scaling == 'Huber2011':
            # According to the relation by Huber et al. 2011b [ppm]
            s = 0.886
            t = 1.89
            r = 2.0
            A_puls_bol = self.L**s * self.M**-t * self.T**(1-r) * A_puls_bol_sun

        elif scaling == 'KjeldsenBedding2011':  # TODO
            # According to Corsaro et al. 2013, equ. (24)
            tau_puls = convert('d', 'Ms', np.exp(delta_Teff / 601.) * 2.65)
            r = 2.0
            # According to Kjeldsen and Bedding 2011 equ. (6), [ppm]
            tau0 = convert('d','Ms',2.65)
            A_puls_bol = self.L * (tau_puls/tau)**0.5 * self.M**-1.5 * self.T**-(1.25+r) * 3.6

        elif scaling == 'Corsaro2013huber2011':
            # Reestimated of Huber et al. (2011b)'s values by
            # Corsaro et al. (2013) using Bayesian inference
            s = 0.984
            t = 1.66
            r = 2.79
            # According to the relation by Huber et al. 2011b,
            # and the fit by Corsaro et al. (2013), Eq. (19) [ppm]
            A_puls_bol = (numax_frac**(2*s - 3*t) *
                          deltanu_frac**(-4*s + 4*t) *
                          self.T**(5*s - 1.5*t - r + 0.2) * 3.6)
            
        elif scaling == 'Corsaro2013':
            # According to Corsaro et al. 2013, Eq. 24:
            # NOTE The value of T0 was calibrated using Kepler RGs in the open clusters
            # NGC 6791 and NGC 6819, and a sample of MS and subgiant Kepler field stars.
            T0       = 601.                            # [K]
            tau0     = (2.65 * u.d).to('Ms').value     # [Ms]
            tau_puls = np.exp(delta_Teff / T0) * tau0  # [Ms]
            # According to Corsaro et al. (2013) Eq. 26 (Model 6) [ppm]
            r = -2.8
            t = 1.56
            A_puls_bol = (numax_frac**(2 - 3*t) *
                          deltanu_frac**(4*t - 4) *
                          (tau_puls/tau_puls_sun)**0.5 *
                          self.T**(4.55 - r - 1.5*t) * A_puls_bol_sun)
            
        # Amplitude [ppm]
        self.ampl = A_puls_bol * np.exp(-(self.freq-self.numax)**2/(2*(1.5*self.deltanu)**2))
        
        # Return model parameters
        return self.numax, self.deltanu, self.freq, self.ampl

            
            
    def eval_oscillations(self, time=False, freq=False, ampl=False, eta=False):

        """Compute time series of stochastically excited damped modes.
        
        Parameters
        ----------
        time : ndarray
            Time points [0..Ntime-1] (unit: e.g. Ms)
        freq : ndarray 
            Oscillation freqs [0..Nmodes-1] (unit: e.g. microHz)
        ampl : ndarray 
            Amplitude of each oscillation mode rms amplitude = ampl / sqrt(2.)
        eta : ndarray
            Damping rates (unit: e.g. (Ms)^{-1})

        Returns
        -------
        signal : ndarray
            Pulsation signal[0..Ntime-1]

        Resources
        ---------
        De Ridder et al., 2006, MNRAS 365, pp. 595-605.
        """

        # Parsing of arguments
        if not time: time = self.time
        if not freq: freq = self.freq
        if not ampl: ampl = self.ampl
        if not eta:  eta  = self.eta

        # Constants
        Ntime = len(time)
        Nmode = len(freq)

        # Set the kick (= reexcitation) timestep to be one 100th of the
        # shortest damping time. (i.e. kick often enough).
        kicktimestep = (1.0 / max(eta)) / 100.0

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
        for i in range(Nwarmup):
            amplsin = damp * amplsin + np.random.normal(np.zeros(Nmode), kick_amplitude)
            amplcos = damp * amplcos + np.random.normal(np.zeros(Nmode), kick_amplitude)

        # Initialize the last kick times for each mode to be randomly chosen
        # a little before the first user time point. This is to avoid that
        # the kicking time is always exactly the same for all of the modes.

        last_kicktime = np.random.uniform(time[0] - kicktimestep, time[0], Nmode)
        next_kicktime = last_kicktime + kicktimestep

        # Generate and return model [ppm -> normalised]
        self.signal_puls = pulsations(time, freq, eta, Ntime, Nmode, amplsin, amplcos,
                                      kicktimestep, kick_amplitude,
                                      last_kicktime, next_kicktime)
        return self.signal_puls

            
# Numba cannot be integrated into a class currently
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





#==============================================================#
#                         MASSIVE STARS                        #
#==============================================================#

    
class SurfaceModulations(object):

    """Class to generate variability of roAp stars.

    Rotationally modulated chemcically perculiar A-type (roAp) stars
    are common and readily available objects that efficiently can be 
    used as photometric calibrators in larger surveys. Rotational 
    periods can be drawn randomly from a uniform distribution between
    1 and 3 days. The shape of light curves of these stars is typically
    sinusoidal with frequent contribution of the second harmonic. 
    The amplitude is typically 10-30 mmag in the Kepler passband. 
    We here assume that Kepler passband is representative for the 
    PLATO passband.
    """

    def __init__(self, time, seed=False):
        
        self.time = time
        self.rng  = ut.rng(seed)
    

        
    def initToyModel(self, period_range=[1,3], amplitude_range=[10,30]):

        """Draw pulsations from uniform distribution.
        Author: Oleg Kochukhov (oleg.kochukhov@physics.uu.se)

        Typical ratational periods [1, 3] days
        Typical amplitude ranges [10, 30] mmag
        """
        
        # Random value of rotational period
        P = self.rng.uniform(period_range[0], period_range[1])

        # Relative amplitudes between main frequency and harmonic
        A = self.rng.uniform(0, 1)

        # Random phase offset [0, 2*pi] between main frequency and harmonic
        phi = self.rng.uniform(0, 2*np.pi)
        
        # Create light curve
        mag = np.cos(2*np.pi * (self.time/P)) + A * np.cos(4*np.pi * (self.time/P) + phi)

        # Scale amplitude within the range
        scale = self.rng.uniform(amplitude_range[0], amplitude_range[1]) / np.abs(mag.max())
        self.mag = mag * scale

        # Return model parameters
        return [P, phi, A, scale]


        
    def evaluate(self, plot=False):

        """Evaluate and return generated model.
        """        
        
        if plot:
            plt.figure(figsize=(10, 4))
            plt.plot(self.time, self.mag, 'k-')
            plt.xlabel('Time [d]')
            plt.ylabel(r'$\mathcal{P}$ [mmag]')
            plt.xlim(self.time.min(), self.time.max())
            plt.tight_layout()
            plt.show()
        
        return self.mag/1e3




#==============================================================#
#                        PULSATING STARS                       #
#==============================================================#


class GravityOscillator(object):

    """Class to generate gravity oscillation time series.
    """

    def __init__(self, time, power, seed=False):

        self.time  = time
        self.power = power

        # Random number generator
        self.rng = ut.rng(seed)
            

        
    def initToyModel(self, period_range, amplitude_range, nmodes=False):

        """Draw pulsations from uniform distribution.

        Parameters
        ----------
        period_range : list
            Range of pulsation periods [Pmin, Pmax] in unit of [days]
        amplitude_range : list
            Range of pulsation amplitudes [Amin, Amax] in units of [mmag]
        nmodes : int
            Number of pulsation modes to include

        TODO Model gives wrong results! High value at 0 c/d!
        """
        
        # Number of pulsation modes
        if not nmodes:
            nmodes = int(self.rng.normal(50, 5))
        
        # Generate pulsations using uniform distributions
        self.freq  = self.rng.uniform(0, 2*np.pi, nmodes)
        self.ampl  = self.rng.uniform(period_range[0], period_range[1], nmodes)
        self.phase = self.rng.uniform(amplitude_range[0], amplitude_range[1], nmodes)
        self.starname = 'g-Dor'

        
        
    def initGang2020(self, odir, starID=None):

        """Draw frequencies from Kepler g-Dor legacy.
        """

        # Name of folder on FTP server
        filename = 'varsource_gdor_gang2020'
        dataDir  = Path(f'{odir}/{filename}')
        
        # Select a random object from the list and load Fourier data
        if dataDir.is_dir():
            filenames = glob.glob(f'{odir}/{filename}/*.dat')
        else:
            zipfile = f'{filename}.zip'
            print(f'Downloading {zipfile} files..')
            ut.downloadFromFTP(filename=zipfile, outputDir=odir, server='plato')
            os.system(f'unzip {odir}/{zipfile} -d {odir}')
            os.system(f'rm {odir}/{zipfile}')

        # If requested, select specific star or else do a random draw
        if starID is None:
            starfile = Path(self.rng.choice(filenames))
        else:
            starfile = Path(filenames[starID-1])
            
        # Load file containing columns
        df = pd.read_csv(starfile, sep=' ', comment='#',
                         names=['freq', 'ampl', 'phase', 'snr'])

        # Else load the data
        self.freq  = 1/df.freq      # [days]
        self.ampl  = df.ampl * 1e3  # [mag]
        self.phase = df.phase       # [rad]
        self.snr   = df.snr
        self.starname = starfile.name

        

    def evaluate(self, plot=False):

        """Evaluate and return generated model.
        """

        return ns.timeSeriesFromFourier(self.time, self.freq, self.ampl, self.phase,
                                        plot=plot, title=self.starname)
   




class ClassicalPulsator(object):

    """Class to generate variability classical pulsators
    """

    def __init__(self, time, seed=False):
        
        self.time = time
        self.rng  = ut.rng(seed)


        
    def initFromFile(self, odir, starID=None):

        """Draw frequencies from Kepler g-Dor legacy.
        """

        # Name of folder on FTP server
        filenames = glob.glob(f'{self.idir}/RRL-CEP/*.fou')
        starfile = random.choice(filenames)
        
        # Select a random object from the list and load Fourier data
        if dataDir.is_dir():
            filenames = glob.glob(f'{odir}/{filename}/*.dat')
        else:
            zipfile = f'{filename}.zip'
            print(f'Downloading {zipfile} files..')
            ut.downloadFromFTP(filename=zipfile, outputDir=odir, server='plato')
            os.system(f'unzip {odir}/{zipfile} -d {odir}')
            os.system(f'rm {odir}/{zipfile}')

        # If requested, select specific star or else do a random draw
        if starID is None:
            starfile = Path(self.rng.choice(filenames))
        else:
            starfile = Path(filenames[starID-1])
            
        # Load file containing columns
        df = pd.read_csv(starfile, sep=' ', comment='#',
                         names=['freq', 'ampl', 'phase', 'snr'])

        # Else load the data
        self.starname  = starfile.name
        self.period    = df.freq        # [days]
        self.amplitude = df.ampl * 1e3  # [mag]
        self.phase     = df.phase       # [rad]
        self.snr       = df.snr


        
    def evaluate(self, plot=False):

        """Evaluate and return generated model.
        """

        return ns.timeSeriesFromFourier(self.time, self.freq, self.ampl, self.phase,
                                        plot=plot, title=self.starname)






#==============================================================#
#                         OTHER OBJECTS                        #
#==============================================================#



class EclipsingBinary(object):

    """Models Eclipsing Binaries (EBs).
    """


    def __init__(self, time, seed=None):

        """Open the HDF5 output file
        """

        self.time = time
        self.rng  = ut.rng(seed)



    def read_parameters_hdf5(self, file_name, verbose=False):

        """Read the full model parameters of the linear, sinusoid
        and eclipse models to an hdf5 file.

        Parameters
        ----------
        file_name: str
            File name (including path) for loading the results.
        verbose: bool
            If set to True, this function will print some information.

        Returns
        -------
        results: dict
            Contains:
            sin_mean: None, list[numpy.ndarray[float]]
                Parameter mean values for the linear and sinusoid model in the order they appear below.
                linear parameters: const, slope,
                sinusoid parameters: f_n, a_n, ph_n
            sin_err: None, list[numpy.ndarray[float]]
                Parameter error values for the linear and sinusoid model in the order they appear below.
                linear parameters: c_err, sl_err,
                sinusoid parameters: f_n_err, a_n_err, ph_n_err
            sin_hdi: None, list[numpy.ndarray[float]]
                Parameter hdi values for the linear and sinusoid model in the order they appear below.
                linear parameters: c_hdi, sl_hdi,
                sinusoid parameters: f_n_hdi, a_n_hdi, ph_n_hdi
            sin_select: None, list[numpy.ndarray[bool]]
                Sinusoids that pass certain selection criteria
                passed_sigma, passed_snr, passed_h
            ephem: None, numpy.ndarray[float]
                Ephemerides of the EB, p_orb and t_zero
            ephem_err: None, numpy.ndarray[float]
                Error values for the ephemerides, p_err and t_zero_err
            ephem_err: None, numpy.ndarray[float]
                Hdi values for the ephemerides, p_hdi and t_zero_hdi
            phys_mean: None, numpy.ndarray[float]
                Parameter mean values for the physical eclipse model in the order they appear below.
                ecosw, esinw, cosi, phi_0, log_rr, log_sb,
                extra parametrisations: e, w, i, r_sum, r_rat, sb_rat
            phys_err: None, numpy.ndarray[float]
                Parameter error values for the physical eclipse model in the order they appear below.
                ecosw_err, esinw_err, cosi_err, phi_0_err, log_rr_err, log_sb_err,
                Extra parametrisation: e_err, w_err, i_err, r_sum_err, r_rat_err, sb_rat_err
            phys_hdi: None, numpy.ndarray[float]
                Parameter hdi values for the physical eclipse model in the order they appear below.
                ecosw_hdi, esinw_hdi, cosi_hdi, phi_0_hdi, log_rr_hdi, log_sb_hdi,
                Extra parametrisations: e_hdi, w_hdi, i_hdi, r_sum_hdi, r_rat_hdi, sb_rat_hdi
            timings: None, numpy.ndarray[float]
                Eclipse timings of minima and first and last contact points, internal tangency
                and eclipse depth of the primary and secondary:
                t_1, t_2, t_1_1, t_1_2, t_2_1, t_2_2, t_b_1_1, t_b_1_2, t_b_2_1, t_b_2_2, depth_1, depth_2
            timings_err: None, numpy.ndarray[float]
                Parameter error values for the eclipse timings and depths:
                t_1_err, t_2_err, t_1_1_err, t_1_2_err, t_2_1_err, t_2_2_err,
                t_b_1_1_err, t_b_1_2_err, t_b_2_1_err, t_b_2_2_err, depth_1_err, depth_2_err
            timings_hdi: None, numpy.ndarray[float]
                Parameter hdi values for the eclipse timings and depths:
                t_1_hdi, t_2_hdi, t_1_1_hdi, t_1_2_hdi, t_2_1_hdi, t_2_2_hdi,
                t_b_1_1_hdi, t_b_1_2_hdi, t_b_2_1_hdi, t_b_2_2_hdi, depth_1_hdi, depth_2_hdi
            timings_indiv_err: None, numpy.ndarray[float]
                Parameter error values for the individual eclipse timings and depths:
                t_1_err, t_2_err, t_1_1_err, t_1_2_err, t_2_1_err, t_2_2_err,
                t_b_1_1_err, t_b_1_2_err, t_b_2_1_err, t_b_2_2_err, depth_1_err, depth_2_err
            var_stats: None, list[union(float, numpy.ndarray[float])]
                Variability level diagnostic statistics
                std_1, std_2, std_3, std_4, ratios_1, ratios_2, ratios_3, ratios_4
            stats: None, list[float]
                Some statistics: t_tot, t_mean, t_mean_s, t_int, n_param, bic, noise_level
            i_sectors: numpy.ndarray[int]
                Pair(s) of indices indicating the separately handled timespans
                in the piecewise-linear curve.
            text: list[str]
                Some information about the file and data:
                identifier, data_id, description and date_time
        """
        # check some input
        ext = os.path.splitext(os.path.basename(file_name))[1]
        if (ext != '.hdf5'):
            file_name = file_name.replace(ext, '.hdf5')
        # create the file
        with h5py.File(file_name, 'r') as file:
            identifier = file.attrs['identifier']
            description = file.attrs['description']
            data_id = file.attrs['data_id']
            date_time = file.attrs['date_time']
            t_tot = file.attrs['t_tot']
            t_mean = file.attrs['t_mean']
            t_mean_s = file.attrs['t_mean_s']
            t_int = file.attrs['t_int']
            n_param = file.attrs['n_param']
            bic = file.attrs['bic']
            noise_level = file.attrs['noise_level']
            # orbital period and time of deepest eclipse
            p_orb = np.copy(file['p_orb'])
            t_zero = np.copy(file['t_zero'])
            # the linear model
            # y-intercepts
            const = np.copy(file['const'])
            c_err = np.copy(file['c_err'])
            c_hdi = np.copy(file['c_hdi'])
            # slopes
            slope = np.copy(file['slope'])
            sl_err = np.copy(file['sl_err'])
            sl_hdi = np.copy(file['sl_hdi'])
            # sector indices
            i_sectors = np.copy(file['i_sectors'])
            # the sinusoid model
            # frequencies
            f_n = np.copy(file['f_n'])
            f_n_err = np.copy(file['f_n_err'])
            f_n_hdi = np.copy(file['f_n_hdi'])
            # amplitudes
            a_n = np.copy(file['a_n'])
            a_n_err = np.copy(file['a_n_err'])
            a_n_hdi = np.copy(file['a_n_hdi'])
            # phases
            ph_n = np.copy(file['ph_n'])
            ph_n_err = np.copy(file['ph_n_err'])
            ph_n_hdi = np.copy(file['ph_n_hdi'])
            # passing criteria
            passed_sigma = np.copy(file['passed_sigma'])
            passed_snr = np.copy(file['passed_snr'])
            passed_b = np.copy(file['passed_b'])
            passed_h = np.copy(file['passed_h'])
            # the physical eclipse model parameters
            ecosw = np.copy(file['ecosw'])
            esinw = np.copy(file['esinw'])
            cosi = np.copy(file['cosi'])
            phi_0 = np.copy(file['phi_0'])
            log_rr = np.copy(file['log_rr'])
            log_sb = np.copy(file['log_sb'])
            # some alternate parametrisations
            e = np.copy(file['e'])
            w = np.copy(file['w'])
            i = np.copy(file['i'])
            r_sum = np.copy(file['r_sum'])
            r_rat = np.copy(file['r_rat'])
            sb_rat = np.copy(file['sb_rat'])
            # eclipse timings
            t_1 = np.copy(file['t_1'])
            t_2 = np.copy(file['t_2'])
            t_1_1 = np.copy(file['t_1_1'])
            t_1_2 = np.copy(file['t_1_2'])
            t_2_1 = np.copy(file['t_2_1'])
            t_2_2 = np.copy(file['t_2_2'])
            t_b_1_1 = np.copy(file['t_b_1_1'])
            t_b_1_2 = np.copy(file['t_b_1_2'])
            t_b_2_1 = np.copy(file['t_b_2_1'])
            t_b_2_2 = np.copy(file['t_b_2_2'])
            depth_1 = np.copy(file['depth_1'])
            depth_2 = np.copy(file['depth_2'])
            # variability to eclipse depth ratios
            ratios_1 = np.copy(file['ratios_1'])
            ratios_2 = np.copy(file['ratios_2'])
            ratios_3 = np.copy(file['ratios_3'])
            ratios_4 = np.copy(file['ratios_4'])

        sin_mean = [const, slope, f_n, a_n, ph_n]
        sin_err = [c_err, sl_err, f_n_err, a_n_err, ph_n_err]
        sin_hdi = [c_hdi, sl_hdi, f_n_hdi, a_n_hdi, ph_n_hdi]
        sin_select = [passed_sigma, passed_snr, passed_b, passed_h]
        ephem = [p_orb[0], t_zero[0]]
        ephem_err = [p_orb[1], t_zero[1]]
        ephem_hdi = [p_orb[2:4], t_zero[2:4]]
        phys_mean = np.array([ecosw[0], esinw[0], cosi[0], phi_0[0], log_rr[0], log_sb[0],
                              e[0], w[0], i[0], r_sum[0], r_rat[0], sb_rat[0]])
        phys_err = np.array([ecosw[1], esinw[1], cosi[1], phi_0[1], log_rr[1], log_sb[1],
                             e[1], w[1], i[1], r_sum[1], r_rat[1], sb_rat[1]])
        phys_hdi = np.array([ecosw[2:4], esinw[2:4], cosi[2:4], phi_0[2:4], log_rr[2:4], log_sb[2:4],
                             e[2:4], w[2:4], i[2:4], r_sum[2:4], r_rat[2:4], sb_rat[2:4]])
        timings = np.array([t_1[0], t_2[0], t_1_1[0], t_1_2[0], t_2_1[0], t_2_2[0],
                            t_b_1_1[0], t_b_1_2[0], t_b_2_1[0], t_b_2_2[0], depth_1[0], depth_2[0]])
        timings_err = np.array([t_1[1], t_2[1], t_1_1[1], t_1_2[1], t_2_1[1], t_2_2[1],
                                t_b_1_1[1], t_b_1_2[1], t_b_2_1[1], t_b_2_2[1], depth_1[1], depth_2[1]])
        timings_hdi = np.array([t_1[2:4], t_2[2:4], t_1_1[2:4], t_1_2[2:4], t_2_1[2:4], t_2_2[2:4],
                                t_b_1_1[2:4], t_b_1_2[2:4], t_b_2_1[2:4], t_b_2_2[2:4], depth_1[2:4], depth_2[2:4]])
        timings_indiv_err = np.array([t_1[4], t_2[4], t_1_1[4], t_1_2[4], t_2_1[4], t_2_2[4],
                                      t_b_1_1[4], t_b_1_2[4], t_b_2_1[4], t_b_2_2[4], depth_1[4], depth_2[4]])
        var_stats = [ratios_1[0], ratios_2[0], ratios_3[0], ratios_4[0],
                     ratios_1[1:], ratios_2[1:], ratios_3[1:], ratios_4[1:]]
        stats = [t_tot, t_mean, t_mean_s, t_int, n_param, bic, noise_level]
        text = [identifier, data_id, description, date_time]
        # put everything in a dict
        results = {'sin_mean': sin_mean, 'sin_err': sin_err, 'sin_hdi': sin_hdi, 'sin_select': sin_select,
                   'ephem': ephem, 'ephem_err': ephem_err, 'ephem_hdi': ephem_hdi,
                   'phys_mean': phys_mean, 'phys_err': phys_err, 'phys_hdi': phys_hdi,
                   'timings': timings, 'timings_err': timings_err, 'timings_hdi': timings_hdi,
                   'timings_indiv_err': timings_indiv_err,
                   'var_stats': var_stats, 'stats': stats, 'i_sectors': i_sectors, 'text': text}
        if verbose:
            print(f'Loaded analysis file with identifier: {identifier}, created on {date_time}. \n'
                  f'data_id: {data_id}. Description: {description} \n')
        return results


        
    def initIJspeert2023(self, odir, starID=None):

        """Draw frequencies from Kepler g-Dor legacy.
        """

        # Name of folder on FTP server
        filename = 'varsource_EBs_kepler_ijspeert2021'
        dataDir  = Path(f'{odir}/{filename}')
        
        # Select a random object from the list and load Fourier data
        if dataDir.is_dir():
            folders = glob.glob(f'{odir}/{filename}/*')
        else:
            zipfile = f'{filename}.zip'
            print(f'Downloading {zipfile} files..')
            ut.downloadFromFTP(filename=zipfile, outputDir=odir, server='plato')
            os.system(f'unzip {odir}/{zipfile} -d {odir}')
            os.system(f'rm {odir}/{zipfile}')

        # If requested, select specific star or else do a random draw
        if starID is None:
            starDir = Path(self.rng.choice(folders))
        else:
            starDir = Path(filenames[starID-1])
            
        # Load file containing columns
        starfile = starDir / f'{starDir.name}_2.hdf5' 
        result   = self.read_parameters_hdf5(starfile)
        data     = result['sin_mean']

        # Else load the data
        self.freq  = data[2]      # [days]
        self.ampl  = data[3] * 1e3 #df.ampl * 1e3  # [mag]
        self.phase = data[4] #df.phase       # [rad]
        self.starname = str(starDir.stem)



    def evaluate(self, plot=False):

        """Evaluate and return generated model.
        """

        return ns.timeSeriesFromFourier(self.time, self.freq, self.ampl, self.phase,
                                        plot=plot, title=self.starname)
        

        
        


        
class SMBHB(object):

    """Models for Super Massive Black Hole Binaries (SMBHBs).
    """


    def __init__(self, time, seed=None):

        """Open the HDF5 output file
        """

        self.time = time
        self.rng  = ut.rng(seed)
        


        
    def __del__(self):

        """Destructor
        """

        pass



    def initToyModelSpikey(self):
        
        A   = 0.04         # [mag]
        P   = 2 * 365.25   # [days]
        phi = 0.1 * np.pi  # [rad] #np.random.uniform(low=0, high=2*np.pi) 
        tscale    = 10     # [days] #np.ones(len(max)) * 10
        amplitude = 0.08   # [mag]

        return

    
    def initToyModel(self):

        """Simplified description of a SMBH binary system.

        The model consist of two components:
          1) Doppler beaming
          2) Gravitational lens
        
        Assumptions
        -----------
        - Circular orbits (e = 0)
        - The orbital period is drawn from the uniform distribution [1.5, 2.5] years
        - The beaming amplitude is drawn from the uniform distribution [0, 0.1] mag
        - The lensing amplitude is twice that of the beaming amplitude
        - The lensing duration scales linearly with the orbital period [5, 12] days
        """

        # Orbital period [days]
        P_min = 1.5
        P_max = 2.5
        year  = 365.25
        P = self.rng.uniform(P_min, P_max) * year

        # Beaming amplitude [mag]
        A_beam = self.rng.uniform(0, 0.1)
        
        # Lens paramters [mmag and days]
        A_lens = 2 * A_beam

        # Lens time maximum [days]
        phi_lens  = self.rng.uniform(0, 2*np.pi)
        tmax_lens = P * (1 - phi_lens/(2*np.pi))
        
        # Lens duration [days]
        tdur_lens = ut.evalLinReg(np.array([P_min, P_max]),
                                  np.array([5, 10]), P/year)

        # Return parameters
        return P, A_beam, A_lens, phi_lens, tmax_lens, tdur_lens



    def evalLensingEvent(self, tmax=False, tscale=False, ampl=30, asym=1):

        """Simple analytic model for lensing event.
        
        Parameters
        ----------
        tmax : float, ndarray
            Time of maximum intensity due to lensing [days]
        tscale : float, ndarray
            Time scale duration of the flare(s) [days]
        ampl : float, ndarray
            Amplitude of the flare(s) [mmag]
        asym : float, ndarray
            Asymmetry factor of the flare(s)
        """

        # Convert units of input parameters
        time = self.time
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
            b = -1.941 - 0.175 + 2.246 + 1
            c = 1 - 0.689

            # Loop over every time-step in the flare time interval

            for i in range(len(t)):

                # Rise of lensing
                if t[i] <= 0:
                    flux[i] += 0.689 * np.exp(+1.6 * t[i])

                # Decay of lensing
                elif t[i] > 0:
                    flux[i] += 0.689 * np.exp(-1.6 * t[i])

                # Outside lens
                else:
                    flux[i] += 0

        # Return relative flux
        return flux * ampl + 1

    

    def evalToyModel(self, P, A_beam, A_lens, phi, tmax, tdur):

        """Uniform distribution of toy model.
        """
        
        # Model doppler beaming effect
        flux_beam = A_beam * np.sin(2*np.pi * (self.time/P) + phi) + 1

        # Model all lesing evens with P < duration of observation
        flux_lens = np.ones_like(flux_beam)
        for tm in [tmax-P, tmax, tmax+P]:
            flux = self.evalLensingEvent(tm, tdur, A_lens, asym=1)
            flux_lens += flux - np.ones_like(flux_beam)

        # Signals
        self.flux = flux_beam + flux_lens - 1
        self.flux_beam = flux_beam
        self.flux_lens = flux_lens
            
        # Return models
        return self.flux, self.flux_beam, self.flux_lens



    def plot(self):

        fig, ax = plt.subplots(1, 1, figsize=(9,4))
        ax.plot(self.time, self.flux_beam, '-', c='green')
        ax.plot(self.time, self.flux_lens, '-', c='orange')
        ax.plot(self.time, self.flux,      '-', c='royalblue')
        ax.set_xlabel('Time [d]')
        ax.set_ylabel(r'Relative flux')
        ax.set_xlim(self.time.min(), self.time.max())
        plt.tight_layout()
        plt.show()
    
    
#==============================================================#
#                         EXOPLANETS                           #
#==============================================================#


class LimbDarkening(funcFit.OneDFit):

    """Class for fitting the Limb Darkening Coefficients.
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
            data_int_VTA[:,jj], vind = pyasl.specAirVacConvert(wvl_int, data_int[:,jj],
                                                               direction="vactoair")

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

    """Model Doppler Beaming Effect (also called Boosting).

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
        #Teff  = self['Teff'].to('K')

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

    """Class for the calculation of the ellipsoidal effect.

    The model presented here is a fusion between the simple analytical description
    presented by Sphorer (2019) and adding a parametization solving the elliptical
    model. Higher order terms from Morris & Nafilan (1993) may be desired to have
    a look at for future improvemenets.

    Resources:
    ----------
    Sphorer (2019)          : https://arxiv.org/abs/1703.00496
    Murray & Correia (2011) : https://arxiv.org/abs/1009.1738v2
    Morris & Nafilan (1993) : http://cdsads.u-strasbg.fr/pdf/1993ApJ...419..344M
    Claret & Bloemen (2011) : TODO this is for TESS but we don't have g fror PLATO
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


    


class PlanetMRforecast():

    """Class to forecast the mass from a planets radius.
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
            ut.downloadFromFTP(hyper_file.name, hyper_file.parents[0], server='plato')

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
