#!/usr/bin/env python3

"""
Python module to with astro query functions used by "picsim".
"""

import os
import sys
import glob
import h5py
import math
import inspect
import pathlib
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable

from pylab import MaxNLocator
from colorama import Fore, Style
from prettytable import PrettyTable
from scipy.ndimage import median_filter

import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.simbad import Simbad
from astroquery.mast import Catalogs
from astroquery.gaia import Gaia




def ticQuery(star, radius=2, Vmax=18, outFile=None):

    """Query TIC catalog for stars around a given named source below a given V magnitude.

    Parameters
    ----------
    offs : float, optional; <0>
        Horizontal movement additional to default.
    dig : int, optional; <0>
        Number of decimals after the comma.
    side : string, optional; {<'left'>, 'right'}
        To choose the side of the y-axis notation.
    omit_last : bool, optional; <False>
        If True, the top y-axis-label is omitted.

    Returns
    -------
    locs : list
        List of y-tick locations.

    Note
    ----
    This is kind of a non-satisfying hack, which should be handled more
    properly. But it works. Functions to look at for a better implementation:
    ax.ticklabel_format
    ax.yaxis.major.formatter.set_offset_string
    """

    # Get the ticks
    locs, _ = plt.yticks()

    # Put the last entry into a string, ensuring it is in scientific notation
    # E.g: 123456789 => '1.235e+08'
    llocs = '%.3e' % locs[-1]

    # Get the magnitude, hence the number after the 'e'
    # E.g: '1.235e+08' => 8
    yoff = int(str(llocs).split('e')[1])

    # If omit_last, remove last entry
    if omit_last:
        slocs = locs[:-1]
    else:
        slocs = locs

    # Set ticks to the requested precision
    form = r'$%.' + str(dig) + 'f$'
    plt.yticks(locs, list(map(lambda x: form % x, slocs / (10 ** yoff))))

    # Define offset depending on the side
    if side == 'left':
        x_offs = -.18 - x_offs  # Default left: -0.18
    elif side == 'right':
        x_offs = 1 + x_offs  # Default right: 1.0

    # Plot the exponent
    plt.text(x_offs, y_offs, r'$\times10^{%i}$' % yoff, transform=
    plt.gca().transAxes, verticalalignment='top')

    # Return the locs
    return locs


def ticQuery(star, radius=2, Vmax=18, outFile=None):
    """
    Query TIC catalog for stars around a given named source below a given V magnitude.

    Parameters
    ----------
    star : str
        Name of the star to query around.
    radius : float
        Radius in arcmin to query around the star.
    Vmax : float
        Maximum V magnitude to query for.
    outFile : str
        Path of the output file to write to. If None, no file is written.

    Returns
    -------
    results : pandas.DataFrame
        DataFrame containing the results of the query. The named star will appear first if
        it is not removed by the Vmax cut.
    """


    # Get the coordinates of the star from Simbad
    
    result_table = Simbad.query_object(star)
    if result_table is None:
        raise ValueError(f"Could not find {star} in Simbad.")
    ra = result_table["RA"][0]
    dec = result_table["DEC"][0]
    coords = SkyCoord(ra, dec, unit=(u.hourangle, u.deg))

    # Query TIC for stars around the star, within the given radius
    
    results = Catalogs.query_region(coords, radius=radius * u.arcmin, catalog="TIC")
    if results is None:
        raise ValueError(f"Could not find any stars in TIC around {star}.")

    # Convert the results to a Pandas DataFrame
    
    results = results.to_pandas()
    results = results[results["Vmag"] < Vmax][["ra", "dec", "Vmag"]]

    # Optionally write the results to a txt file
    
    if outFile is not None:
        with open(outFile, "w") as f:
            f.write("# RA DEC Vmag\n")
            for i, row in results.iterrows():
                f.write(f"{row['ra']:.6f} {row['dec']:.6f} {row['Vmag']:.3f}\n")

    # That's it!
                
    return results





def gaiaQuery(star):
    """
    Query Gaia for a named star and return the Gaia DR2 ID.

    Parameters
    ----------
    star : str
        Name of the star to query around.

    Returns
    -------
    gaia_id : int
        The Gaia DR2 ID of the star.
    """

    # Query for target star
    result_table = Simbad.query_objectids(star)
    if result_table is None:
        raise LookupError(f"No Simbad results for {star} (probably not a star)")

    # If requested find all stars with 'radius'
    for row in result_table:
        if 'Gaia DR2' in row['ID']:
            gaia_id = row['ID']
            return int(gaia_id[9:])
    raise LookupError(f"No Gaia DR2 ID for {star} (probably a multiple star)")






def starQuery(star, radius=45):
    """
    Query Gaia for a named star and return the Gaia DR2 ID.

    Parameters
    ----------
    star : str
        Name of the star to query around.

    Returns
    -------
    gaia_id : int
        The Gaia DR2 ID of the star.
    """

    # Qucik check that target star exist
    result_table = Simbad.query_objectids(star)
    if result_table is None:
        raise LookupError(f"No Simbad results for {star} (probably not a star)")

    # Fetch the equatorial coordinates
    Simbad.reset_votable_fields()
    Simbad.remove_votable_fields('coordinates')
    Simbad.add_votable_fields('ra(:;A;ICRS;J2000)', 'dec(:;D;ICRS;2000)')
    table = Simbad.query_object(star, wildcard=False)
    coord = SkyCoord(ra=['{}h{}m{}s'.format(*ra.split(':')) for ra in table['RA___A_ICRS_J2000']], 
                     dec=['{}d{}m{}s'.format(*dec.split(':')) for dec in table['DEC___D_ICRS_2000']],
                     frame='icrs', equinox='J2000')
    raStar  = coord.ra.degree[0]
    decStar = coord.dec.degree[0]

    # Convert radius to from arcsec to deg
    radius /= 3600.

    # Gaia radius querymetric
    query_cone = f"""SELECT 
    TOP 10 
    source_id, ra, dec, phot_g_mean_mag
    FROM gaiadr2.gaia_source
    WHERE 1=CONTAINS(POINT(ra, dec), CIRCLE({raStar}, {decStar}, {radius}))
    """

    # Launch Gaia query 
    job     = Gaia.launch_job(query_cone)
    results = job.get_results()

    # Convert astropy results table into pandas df
    df = results.to_pandas()

    # Make sure that target is the first entry
    for row in result_table:
        if 'Gaia DR2' in row['ID']:
            gaia_id = int(row['ID'][9:])            
    row = df.index[df['source_id'] == gaia_id].tolist()
    dex = row + [i for i in range(len(df)) if i != row[0]]
    df = df.iloc[dex].reset_index(drop=True)
    
    # Calculate radial distance
    dist = np.sqrt( (df.ra - df.ra.iloc[0])**2 + (df.dec - df.dec.iloc[0])**2 ) * 3600
    df['dis'] = dist
    df = df.sort_values(by=['dis'])

    # Return data frame
    return df
