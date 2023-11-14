#!/usr/bin/env python3

"""
Python module to with astro query functions used by picsim.py.

NOTE This class needs the Poetry install: 
     >> poetry install --with platonium 
"""

# Built-in
import os
import sys
import glob
import math
import time
import inspect
import pathlib
import http.client as httplib
import urllib.parse as urllib
from xml.dom.minidom import parseString

# PlatoSim standard
import h5py
import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from astropy.io.votable  import parse
from astropy.coordinates import SkyCoord
from mpl_toolkits.axes_grid1 import make_axes_locatable
from pylab import MaxNLocator
from colorama import Fore, Style
from prettytable import PrettyTable
from scipy.ndimage import median_filter

# PLATOnium extra
from astroquery.simbad   import Simbad
from astroquery.mast     import Catalogs
from astroquery.gaia     import Gaia

import platosim.utilities as ut


#--------------------------------------------------------------#
#                          FUNCTIONS                           #
#--------------------------------------------------------------#


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

    """Query TIC catalog for stars around a given named source below a givenV.

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

    """Query Gaia for a named star and return the Gaia DR2 ID.

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






def simbadQuery(star, radius=45, maglim=21):

    """Query Gaia for a named star and return the Gaia DR2 ID.

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

    query_cone = f"""SELECT 
    DISTANCE( POINT({raStar},{decStar}), POINT(ra,dec) )
    AS dist, source_id, ra, dec,
    phot_g_mean_mag, bp_rp,
    parallax, parallax_error,
    pmra, pmdec, ruwe,
    teff_gspphot, logg_gspphot
    FROM gaiadr3.gaia_source AS cat
    WHERE 1=CONTAINS(POINT({raStar}, {decStar}),
    CIRCLE(cat.ra, cat.dec, {radius}))
    AND cat.phot_g_mean_mag < {maglim} 
    ORDER BY dist ASC
    """

    # Launch Gaia query 

    job     = Gaia.launch_job(query_cone)
    results = job.get_results()

    # Convert astropy results table into pandas df

    df = results.to_pandas()

    # Rename columns
    
    df = df.rename(columns={'source_id': 'gaiaDR3',
                            'phot_g_mean_mag': 'Gmag',
                            'bp_rp': 'BP_RP',
                            'parallax': 'plx',
                            'parallax_error': 'plxe',
                            'teff_gspphot': 'teff',
                            'logg_gspphot': 'logg'})

    
    # Make sure that target is the first entry    
    # for row in result_table:
    #     if 'source_id' in row['ID']:
    #         gaia_id = int(row['ID'][9:])            
    # row = df.index[df['source_id'] == gaia_id].tolist()
    # dex = row + [i for i in range(len(df)) if i != row[0]]
    # df = df.iloc[dex].reset_index(drop=True)

    # Convert Gmag to Pmag

    df['Pmag'] = ut.passbandConversionG2P(df.Gmag, df.BP_RP)
    
    # Relocate distance column [arcsec]
    
    df.dist = (df.dist-df.dist.iloc[0]) * 3600.
    #col = df.dist.values.tolist()
    #df  = df.drop(columns=['dist'])
    #df.insert(5, 'dist', col)

    # Sort and return
    
    return df.sort_values(by=['dist'])





def gaiaRegionQuerySmall(alpha, delta, radius=1, maglim=21):

    """Query sky region using Gaia DR3.

    Parameters
    ----------
    star : str
        Name of the star to query around.

    Returns
    -------
    gaia_id : int
        The Gaia DR2 ID of the star.
    """

    # Gaia query job cone

    query_cone = f"""SELECT 
    DISTANCE( POINT({alpha},{delta}), POINT(ra,dec) )
    AS dist, source_id, ra, dec,
    phot_g_mean_mag, bp_rp,
    parallax, parallax_error,
    pmra, pmdec, ruwe,
    teff_gspphot, logg_gspphot
    FROM gaiadr3.gaia_source AS cat
    WHERE 1=CONTAINS(POINT({alpha}, {delta}),
    CIRCLE(cat.ra, cat.dec, {radius}))
    AND cat.phot_g_mean_mag < {maglim} 
    ORDER BY dist ASC
    """

    # Launch Gaia query
    
    job     = Gaia.launch_job(query_cone)
    results = job.get_results()
    
    # Convert astropy results table into pandas df

    df = results.to_pandas()

    # Rename columns
    
    df = df.rename(columns={'source_id': 'gaiaDR3',
                            'phot_g_mean_mag': 'mag',
                            'bp_rp': 'BP_RP',
                            'parallax': 'plx',
                            'parallax_error': 'plxe',
                            'teff_gspphot': 'teff',
                            'logg_gspphot': 'logg'})
    
    # Relocate distance column [arcsec]
    
    df.dist = (df.dist-df.dist.iloc[0]) * 3600.
    col = df.dist.values.tolist()
    df  = df.drop(columns=['dist'])
    df.insert(5, 'dist', col)

    #
    df = df.fillna(0)
    df = df.astype({'teff':int})
    
    # Sort and return
    
    return df.sort_values(by=['dist'])






def gaiaRegionQuery(ra, dec, radius=1, maglim_min=0, maglim_max=17,
                    flag_astro=False, ofile=False):

    """Function to query a circular sky region from Gaia DR3.
    
    Parameters
    ----------
    ra : float, ndarray
        Right ascension of central point for query [deg]
    dec : float, ndarray
        Declination of central point for query [deg]
    radius : int, float
        Angular radius to search for stars within [deg]
    maglim : int, float
        Magitude limit to search for stars below
    ofile : str
        File name (without file extension) to be saved

    Return
    ------
    Feather file given by the file name destination

    Notes
    -----
    Columns can be found at: https://gea.esac.esa.int/archive/
    See: "Search" -> "Advanced (ADQL)" -> "Gaia Data Release 3"

    NOTE Query only valid for Gmag < 17, else gaps are intoduced!
    Use additional queries e.g. in bins of 1 mag for Gmah > 17.
    """

    # Fetch columns from catalogues

    if flag_astro:
        colname = ['gaia.source_id',
                   'gaia.ra',
                   'gaia.dec',
                   'gaia.phot_g_mean_mag',
                   'gaia.bp_rp',
                   'gaia.ag_gspphot',
                   'gaia.parallax', 'gaia.parallax_error',
                   'gaia.pmra',
                   'gaia.pmdec',
                   'gaia.ruwe',
                   'gaia.mh_gspphot',    'gaia.mh_gspphot_lower',    'gaia.mh_gspphot_upper',
                   'gaia.logg_gspphot',  'gaia.logg_gspphot_lower',  'gaia.logg_gspphot_upper',                   
                   'gaia.teff_gspphot',  'gaia.teff_gspphot_lower',  'gaia.teff_gspphot_upper',
                   'astro.radius_flame', 'astro.radius_flame_lower', 'astro.radius_flame_upper',
                   'astro.mass_flame',   'astro.mass_flame_lower',   'astro.mass_flame_upper',
                   'astro.lum_flame',    'astro.lum_flame_lower',    'astro.lum_flame_upper',
                   'astro.spectraltype_esphs']
        columns = ', '.join(colname)    
        query_base = f"""SELECT
        {columns}
        FROM gaiadr3.gaia_source AS gaia
        JOIN gaiadr3.astrophysical_parameters AS astro
          ON gaia.source_id = astro.source_id
        WHERE 1=CONTAINS(
          POINT(gaia.ra, gaia.dec),
          CIRCLE({ra}, {dec}, {radius}))
          AND gaia.phot_g_mean_mag > {maglim_min}
          AND gaia.phot_g_mean_mag < {maglim_max}        
        """

    else:
        colname = ['gaia.source_id',
                   'gaia.ra',
                   'gaia.dec',
                   'gaia.phot_g_mean_mag',
                   'gaia.bp_rp']
        columns = ', '.join(colname)    
        query_base = f"""SELECT
        {columns}
        FROM gaiadr3.gaia_source AS gaia
        WHERE 1=CONTAINS(
          POINT(gaia.ra, gaia.dec),
          CIRCLE({ra}, {dec}, {radius}))
          AND gaia.phot_g_mean_mag > {maglim_min}
          AND gaia.phot_g_mean_mag < {maglim_max}        
        """
        
    # We use the urllib to keep the connection open because
    # the sky regions are huge which fails with Gaia.lunch_job
    
    params = urllib.urlencode({"REQUEST"        : "doQuery",
                               "LANG"           : "ADQL",
                               "FORMAT"         : "votable_plain",
                               "PHASE"          : "RUN",
                               "JOBNAME"        : "PLATO catalog",
                               "JOBDESCRIPTION" : "None", 
                               "QUERY"          : query_base})
    
    headers = {"Content-type": "application/x-www-form-urlencoded",
               "Accept"      : "text/plain"}

    # Information about server
    
    host = "gea.esac.esa.int"
    port = 443
    pathinfo = "/tap-server/tap/async"
    connection = httplib.HTTPSConnection(host, port)
    connection.request("POST", pathinfo, params, headers)

    # Get status
    
    response = connection.getresponse()
    location = response.getheader("location")
    jobid    = location[location.rfind('/')+1:]
    connection.close()

    # Check job status, wait until finished

    while True:
        connection = httplib.HTTPSConnection(host, port)
        connection.request("GET",pathinfo+"/"+jobid)
        response = connection.getresponse()
        data = response.read()

        # XML response: parse it to obtain the current status
        # (you may use pathinfo/jobid/phase entry point to avoid XML parsing)
        dom = parseString(data)
        phaseElement = dom.getElementsByTagName('uws:phase')[0]
        phaseValueElement = phaseElement.firstChild
        phase = phaseValueElement.toxml()

        # Check finished
        if phase == 'COMPLETED':
            break
        
        # Wait and repeat
        time.sleep(0.2)

    connection.close()

    # Get results
    
    connection = httplib.HTTPSConnection(host, port)
    connection.request("GET",pathinfo+"/"+jobid+"/results/result")
    response = connection.getresponse()
    data = response.read().decode('iso-8859-1')

    # Write output file
    
    if not ofile:
        ofile = 'starcatGaiaDR3.vot'
    outputFile = open(ofile, 'w')
    outputFile.write(data)
    outputFile.close()
    connection.close()

    # Load output file into a pandas df
    
    votable = parse(ofile)
    df = ut.votable2pandas(votable)
    os.remove(ofile)
    
    # Rename columns

    if flag_astro:
        df = df.rename(columns={'source_id': 'gaiaDR3',
                                'phot_g_mean_mag': 'Gmag',
                                'bp_rp': 'BP_RP',
                                'ag_gspphot': 'Ag',
                                'parallax': 'plx',
                                'parallax_error': 'plx_err',
                                'mh_gspphot': 'Z',
                                'mh_gspphot_lower': 'Z_low',
                                'mh_gspphot_upper': 'Z_upp',
                                'logg_gspphot': 'logg',
                                'logg_gspphot_lower': 'logg_low',
                                'logg_gspphot_upper': 'logg_upp',
                                'teff_gspphot': 'Teff',
                                'teff_gspphot_lower': 'Teff_low',
                                'teff_gspphot_upper': 'Teff_upp',
                                'radius_flame': 'R',
                                'radius_flame_lower': 'R_low',
                                'radius_flame_upper': 'R_upp',
                                'mass_flame': 'M',
                                'mass_flame_lower': 'M_low',
                                'mass_flame_upper': 'M_upp',
                                'lum_flame': 'L',
                                'lum_flame_lower': 'L_low',
                                'lum_flame_upper': 'L_upp',                                
                                'spectraltype_esphs': 'spec'})
        df = df.sort_values(by=['BP_RP'])
        # Round Teff column
        df = df.fillna(-1)
        df = df.astype({'Teff':int,
                        'Teff_low':int,
                        'Teff_upp':int})
        df = df.replace({-1:np.nan})
    else:
        df = df.rename(columns={'source_id': 'gaiaDR3',
                                'phot_g_mean_mag': 'Gmag',
                                'bp_rp': 'BP_RP'})

    # Convert names to string
    
    df.gaiaDR3 = df.gaiaDR3.astype(str)
        
    # Reset index and return

    return df.reset_index(drop=True)
