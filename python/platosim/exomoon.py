# Built-in
import os
import sys
import glob
import json
import datetime

# PlatoSim standard
import numpy as np
import scipy as sp
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy import constants as c

# PlatoSim functions
import platosim.plot      as pt
import platosim.utilities as ut 
#from platosim.varsource    import Exomoon
from platosim.lightcurve   import LightCurve
from platosim.utilities    import getFunctions

# Pandora-moon
import pandoramoon as pandora
from pandoramoon.helpers import ld_convert, ld_invert
import ultranest
import ultranest.stepsampler
from ultranest import ReactiveNestedSampler
from ultranest.plot import cornerplot
import corner




def semimajor_axis(M, P):
    """Semimajor axis of system [SI units]
    """
    return (c.G.value * M * P**2 / (4 * np.pi**2))**(1/3)


def mass_moon(R_moon):
    """Mass-radius relation for rocky Earth-sized planets.
    """
    return (R_moon/c.R_earth.value)**3 * c.M_earth.value


def get_n_epochs(P_bary, t0_bary, t_dur):
    """Select data around planet transits.

    Parameters
    ----------
    P : float
        Orbital period [d]
    t0 :float
        Time of first transit [d]
    t_dur : float
        Duration of time series [d]
    """
    return round(np.floor((t0_bary + t_dur) / P_bary))


def cut_out_transits(df, P_bary, t0_bary, T_epoch):
    """Select data around planet transits.

    Parameters
    ----------
    df : pd.DataFrame
        Light curve with columns: 
        time [d] and flux [pp1]
    P : float
        Orbital period [d]
    t0 :float
        Time of first transit [d]
    T_epoch : float
        Duration around transit to make cut [d]
    """
    t_dur = df.time.iloc[-1]
    epochs = round(np.floor( (t0_bary + t_dur) / P_bary))
    df1 = pd.DataFrame()
    for i in range(epochs):
        t_transit = t0_bary + (P_bary * i)
        df0 = df[(df.time > t_transit - T_epoch) & (df.time < t_transit + T_epoch)]
        df1 = pd.concat([df1, df0])
    return df1


def fromMagToFlux(mag):
    """Convert relative magnitude to relative flux.

    Parameters
    ----------
    mag : float
        Input magnitude

    Return
    ------
    flux : ndarray
        Relative flux
    """
    return 10**(-0.4*mag)


def plot_lc(df, dv, di=None, alpha=0.5, figsize=(9,4)):
    """Select data around planet transits.
    
    Parameters
    ----------
    df : light curve with columns: time [d], flux [pp1]
    dv : input model with columns: time [d], flux [pp1]
    """
    fig, ax = plt.subplots(1, 1, figsize=figsize)
    ax.plot(df.time, df.flux, 'k.', ms=10, alpha=alpha)
    ax.plot(dv.time, dv.flux, '-', c='royalblue')
    if di is not None:
        ax.plot(di.time, di.flux, '-', c='orange')
    ax.set_xlabel("Time [days]")
    ax.set_ylabel("Normalized flux")
    plt.tight_layout()
    return fig, ax


def plot_epochs(df, dm, P_bary, t0_bary, T_epoch, n_epoch, alpha=0.5, figsize=(9,3)):
    """Select data around planet transits.

    Parameters
    ----------
    df : pd.DataFrame  
         Simulated light curve with columns: 
         time [d], flux [pp1], flux_err [pp1]
    dm : pd.DataFrame
        Pandora model with columns: 
        time [d], flux [pp1], flux_planet [pp1], flux_moon [pp1]
    n_epoch : float
        Number of transit epochs [d]
    T_epoch : float
        Duration of each epoch [d]
    t0 : float
        Time of first transit [d]
    P : float
        Period of transits [d]
        """
    # Interpolate model to data timings
    if df.shape[0] != dm.shape[0]:
        interp = sp.interpolate.make_interp_spline(dm.time, dm.flux, k=1)
        df['flux_model'] = interp(df.time)
    else:
        df['flux_model'] = dm.flux
        
    # Plot light curve and OC comprarison
    fig, ax = plt.subplots(2, n_epoch, figsize=figsize)        
    for i in range(n_epoch):
        # Plot simulated data
        ax[0,i].errorbar(df.time, df.flux, yerr=df.flux_err, fmt='.k', alpha=alpha, zorder=1)
        # Plot Pandora planet and moon model
        if 'flux_planet' in dm:
            ax[0,i].plot(dm.time, dm.flux_planet, '-', c='deeppink', zorder=2)
        if 'flux_moon' in dm:
            ax[0,i].plot(dm.time, dm.flux_moon,   '-', c='royalblue', zorder=3)
        # Plot injected model
        ax[0,i].plot(dm.time, dm.flux, '-', c='orange', zorder=4)                
        # Plot OC diagram
        ax[1,i].errorbar(df.time, df.flux-df.flux_model, yerr=df.flux_err,
                        fmt='.k', alpha=alpha)
        # Settings
        t_transit = t0_bary + (P_bary * i)
        for j in range(2):
            ax[j,i].set_xlim(t_transit-T_epoch/2 , t_transit+T_epoch/2)
            ax[j,i].set_xlabel("Time [days]")
    ax[0,0].set_ylabel("Normalized flux")
    ax[1,0].set_ylabel("Residuals")
    plt.tight_layout()
    return fig, ax


#--------------------------------------------------------------#
#                     PUBLIC ULTRANET CLASS                    #
#--------------------------------------------------------------#

class model_priors(object):
    """Initialise model priors.
    """
    def __init__(self):
        # Stellar parameters
        self.R_star         = 6.957e8  # [m]
        self.q1             = 0.399
        self.q2             = 0.381
        # Planet parameters
        self.M_planet       = 1.9e27  # [kg]
        self.t0_bary        = 0       # [day]
        self.per_bary       = 365.25  # [day]
        self.a_bary         = 215.1   # [R_star]
        self.r_planet       = 0.1     # [R_star]
        self.b_bary         = 1       # [0, 1]
        self.t0_bary_offset = 0       # [day]
        self.ecc_bary       = 0       # [0, 1]
        self.w_bary         = 0       # [0, 180] deg
        # Moon parameters
        self.M_moon         = 1e-8    # [kg]
        self.r_moon         = 1e-8    # [R_star]
        self.per_moon       = 10      # [day]
        self.tau_moon       = 0       # [0, 1]  
        self.Omega_moon     = 0       # [0, 180] deg
        self.i_moon         = 90      # [0, 90] deg
        self.ecc_moon       = 0       # [0, 1]
        self.w_moon         = 0       # [0, 180] deg
        
        # Other parameters
        self.epoch_distance = 365.25
        

def loglikelihood(x_data, x_err, x_model):
    """Simple chi-square log-likehood.
    """
    return -0.5 * np.nansum(((x_model - x_data) / x_err)**2)


def run_ultranest(df, priors, nsteps=1000, live_points=400, path='results', mode='planet-only'):    

    # Generate arrays used throughout notebook
    time = df.time.to_numpy()
    flux = df.flux.to_numpy()
    flux_err = df.flux_err.to_numpy()

    R_star   = priors.R_star
    q1       = priors.q1
    q2       = priors.q2
    M_planet = priors.M_planet
    t0_bary  = priors.t0_bary
    epoch_distance = priors.epoch_distance

    Rs = priors.R_star
    Mp = priors.M_planet
    Pb = priors.per_bary
    ab = priors.a_bary
    bb = priors.b_bary
    Rp = priors.r_planet
    tb = priors.t0_bary_offset
    Mm = priors.M_moon
    Rm = priors.r_moon
    Pm = priors.per_moon
    tm = priors.tau_moon
    om = priors.Omega_moon
    im = priors.i_moon
    
    # SELECT MODE FOR PLANET-ONLY OR PLANET-MOON
    
    if mode == 'planet-only':
        parameters = [
            'per_bary', 
            'a_bary', 
            'r_planet',
            'b_bary',
            't0_bary_offset',
        ]
        wrapped_params = [
            False,
            False, 
            False, 
            False, 
            False, 
        ]
        def prior_transform(cube):
            p = cube.copy()
            #-------------
            p[0] = cube[0] * (Pb[1] - Pb[0]) + Pb[0]
            p[1] = cube[1] * (ab[1] - ab[0]) + ab[0]   
            p[2] = cube[2] * (Rp[1] - Rp[0]) + Rp[0]   
            p[3] = cube[3] * (bb[1] - bb[0]) + bb[0]   
            p[4] = cube[4] * (tb[1] - tb[0]) + tb[0]
            return p
        def log_likelihood(p):
            # Convert q priors to u LDs (Kipping 2013)
            u1, u2 = ld_convert(q1, q2)
            # Calculate pandora model with trial parameters
            _, _, flux_trial, _, _, _, _ = pandora.pandora(
                R_star = R_star,
                u1     = u1,
                u2     = u2,
                # Planet parameters
                M_planet = M_planet,
                per_bary = p[0],
                a_bary   = p[1],
                r_planet = p[2],
                b_bary   = p[3],
                ecc_bary = 0, #priors.ecc_bary,
                w_bary   = 0, #priors.w_bary,
                t0_bary  = t0_bary,
                t0_bary_offset = p[4],   
                # Moon parameters
                M_moon     = 1e-8,  # Set negligible moon mass
                r_moon     = 1e-8,  # Set negligible moon size
                per_moon   = 10,    # Other moon parameter do not matter
                tau_moon   = 0,
                Omega_moon = 0,
                i_moon     = 0,
                ecc_moon   = 0,
                w_moon     = 0,
                # Other model parameters
                epoch_distance         = epoch_distance,
                supersampling_factor   = 1,
                occult_small_threshold = 0.01,
                hill_sphere_threshold  = 1.0,
                numerical_grid         = 25,
                time                   = time,
            )
            return loglikelihood(flux, flux_err, flux_trial)

    elif mode == 'planet-moon':
        parameters = [
            'R_star', 
            'per_bary', 
            'a_bary', 
            'r_planet',
            'b_bary',
            't0_bary_offset',
            'M_planet',
            'r_moon',
            'per_moon',
            'tau',
            'Omega_moon',
            'i_moon',
            'M_moon',
        ]
        wrapped_params = [
            False,
            False, 
            False, 
            False, 
            False, 
            False, 
            False, 
            False, 
            False,
            True,  # tau -> To save computation time
            False,
            False,  
            False,
        ]
        def prior_transform(cube):
            p = cube.copy()
            #----------------------------------
            p[0]  = cube[0]  * (Rs[1] - Rs[0]) + Rs[0]
            p[1]  = cube[1]  * (Pb[1] - Pb[0]) + Pb[0]
            p[2]  = cube[2]  * (ab[1] - ab[0]) + ab[0]   
            p[3]  = cube[3]  * (Rp[1] - Rp[0]) + Rp[0]   
            p[4]  = cube[4]  * (bb[1] - bb[0]) + bb[0]   
            p[5]  = cube[5]  * (tb[1] - tb[0]) + tb[0]
            #----------------------------------
            p[6]  = cube[6]  * (Mp[1] - Mp[0]) + Mp[0]
            p[7]  = cube[7]  * (Rm[1] - Rm[0]) + Rm[0]
            p[8]  = cube[8]  * (Pm[1] - Pm[0]) + Pm[0]
            p[9]  = cube[9]  * (tm[1] - tm[0]) + tm[0] 
            p[10] = cube[10] * (om[1] - om[0]) + om[0]
            p[11] = cube[11] * (im[1] - im[0]) + im[0]
            p[12] = cube[12] * (Mm[1] - Mm[0]) + Mm[0]
            return p
        def log_likelihood(p):
            u1, u2 = ld_convert(q1, q2)
            # Calculate pandora model with trial parameters
            _, _, flux_trial, _, _, _, _ = pandora.pandora(
                R_star = p[0],
                u1     = u1,
                u2     = u2,
                # Planet parameters
                per_bary       = p[1],
                a_bary         = p[2],
                r_planet       = p[3],
                b_bary         = p[4],
                ecc_bary       = 0,
                w_bary         = 0,
                t0_bary        = t0_bary,
                t0_bary_offset = p[5],   
                M_planet       = p[6],
                # Moon parameters
                r_moon         = p[7],
                per_moon       = p[8],
                tau_moon       = p[9],
                Omega_moon     = p[10],
                i_moon         = p[11],
                ecc_moon       = 0,
                w_moon         = 0,
                M_moon         = p[12],
                # Other model parameters
                epoch_distance         = epoch_distance,
                supersampling_factor   = 1,
                occult_small_threshold = 0.01,
                hill_sphere_threshold  = 1.0,
                numerical_grid         = 25,
                time                   = time,
                #cache=cache  # Can't use cache because free LDs
            )
            return loglikelihood(flux, flux_err, flux_trial)

    # Initialise sampler    
    sampler = ReactiveNestedSampler(
                parameters,
                log_likelihood, 
                prior_transform,
                wrapped_params=wrapped_params,
                log_dir=path,
                resume='overwrite'
                )

    # Set number of step
    sampler.stepsampler = ultranest.stepsampler.RegionSliceSampler(
        nsteps=nsteps,
        max_nsteps=10000,
        adaptive_nsteps='move-distance',
        )  

    # Run nested sampling
    tic  = datetime.datetime.now()
    results = sampler.run(min_num_live_points=live_points)
    toc  = datetime.datetime.now()
    print(f'Execution time: {toc-tic} [h:mm:ss]')

    # Always show the distributions    
    sampler.print_results()
    
    # Return result and sampler
    return results, sampler


def plot_corner(result, bestfit=False, values_input=None):
    """Select data around planet transits.
    
    Parameters
    ----------
    """    
    figure = corner.corner(
        result['samples'],
        smooth=1.5,
        color='royalblue',
        labels=result['paramnames'],
        show_titles=True,
        title_kwargs={"fontsize": 18},
        quantiles=[0.16, 0.5, 0.84],
    )
    if bestfit in ['maximum', 'median', 'mean']:
        if bestfit == 'maximum':
            values_bestfit = np.array(result['maximum_likelihood']['point'])
        else:
            values_bestfit = np.array(result['posterior'][bestfit])
        corner.overplot_lines(figure,  values_bestfit, color='deeppink', lw=1)
        corner.overplot_points(figure, values_bestfit[None], marker="s", color='deeppink')
    if values_input is not None:
        corner.overplot_lines(figure,  values_input, color='orange', lw=1)    
        corner.overplot_points(figure, values_input[None], marker="s", color='orange')
    return figure
