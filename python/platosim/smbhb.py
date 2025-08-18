#!/usr/bin/env python3

"""
This python module contains plot utilities used in the minimal 
PlatoSim installation and in the extra PLATOnium installation.
"""

import time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from astropy.table import Table
from astropy.time import Time
from astropy import units as u
from astropy import constants as c
from astropy.coordinates import SkyCoord
from ipywidgets import *
from tqdm import tqdm
from mpl_toolkits.axes_grid1 import make_axes_locatable

# PlatoSim libraries
import platosim.smbhb     as bh
import platosim.slurm     as sm
import platosim.utilities as ut
import platosim.plot      as pt
import platosim.noise     as ns
import platosim.starquery as sq
from platosim.lightcurve   import LightCurve
from platosim.slurm        import workerOverview
from platosim.matplotlibrc import setup_notebook
setup_notebook()

# Random number generator
rng = ut.rng(12345)

#---------------------------------------------------------------
#  FUNCTIONS FOR NOTEBOOK: 1. Star catalogues
#---------------------------------------------------------------

def fetch_gaia_info(df, NED=False):
    """Fetch Gaia info for each source in data frame.
    Use NASA/IPAC Extragalactic Database (NED).
    """
    for i in range(df.shape[0]):
        di = df.reset_index(drop=True).loc[i]
        if NED:
            ra, dec = di.RA, di.Dec
        else:
            ra, dec = di.ra, di.dec
        dx = sq.gaiaQueryCone(ra, dec, radius=0.001, mag_max=21)
        if dx.shape[0] > 1:
            dx = dx[dx.dis == 0]
        if i == 0:
            dq = dx
        else:
            dq = pd.concat((dq, dx))
            
    # Alter data frame
    dq = dq.reset_index(drop=True)
    dq = dq.drop(columns='dis')
    if NED:
        dq.insert(0, 'source', df['Object Name'])
        dq['z'] = df['Redshift (z)']
        
    return dq


def plot_aitoff(df_agn, df_all=False, df_lop=False, df_best=False, NED=False):
    """Fetch Gaia info for each source in data frame.
    Use NASA/IPAC Extragalactic Database (NED).
    """    
    if df_best is not False:
        df = df_best
    elif df_lop is not False:
        df = df_lop
    elif df_all is not False:
        df = df_all
    else:
        df = df_agn
    title = (f'Total: {df.shape[0]}, ' +
             f'LOPN1: {df[df.b > 0].shape[0]}, ' + 
             f'LOPS2: {df[df.b < 0].shape[0]}')
        
    # Plot PLATO AGNs
    fig, ax = pt.drawStarsInSkyAitoff(
        df_agn.ra, df_agn.dec, column=df_agn.ncam, cbarMap='Blues',
        cbarLabel=r'N-CAM visibility, $n_{\rm NCAM}$',
        title=title, fs=13, figsize=(10,7))

    # Plot all candidates
    if df_all is not False:
        if NED:
            ra, dec = df_all.RA, df_all.Dec
        else:
            ra, dec = df_all.ra, df_all.dec
        gal = SkyCoord(ra, dec, frame='icrs', unit=u.deg).galactic
        ax.scatter(-gal.l.wrap_at('180d').radian, gal.b.radian,
                   c='orange', marker='o', s=10, ec='w', lw=0.8, zorder=4)

    # Plot candidates within LOPs
    if df_lop is not False:
        gal = SkyCoord(ra, dec, frame='icrs', unit=u.deg).galactic
        ax.scatter(-gal.l.wrap_at('180d').radian, gal.b.radian,
                   c='k', marker='o', s=20, ec='w', lw=0.8, zorder=5)

    # Plot best candidates within LOPs
    if df_best is not False:
        gal = SkyCoord(df_best.ra, df_best.dec, frame='icrs', unit=u.deg).galactic
        ax.scatter(-gal.l.wrap_at('180d').radian, gal.b.radian,
                   c='limegreen', marker='o', s=20, ec='w', lw=0.8, zorder=5);
    
    return fig, ax


