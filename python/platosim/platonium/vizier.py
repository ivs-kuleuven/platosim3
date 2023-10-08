#!/usr/bin/env python3

"""
This script creates a coherent PLATO catalogue that covers the
FoV of each camera group. It uses the Gaia DR3 catalogue and
makes query of 9 grid points distribution around the requested
PLATO pointing field. The output catalogues (one for each group)
is to be used directly by "platonium", either for the generation
of full-frame CCD images or for exploitation of the PLATO-CS.
"""

# Python standard
import os
import argparse

# PlatoSim standard
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from astropy.coordinates import SkyCoord
from astropy import units as u

# PLATOnium extra
from pathlib import Path
from tqdm import tqdm
import ligo.skymap.plot

# PlatoSim functions
import platosim.utilities       as ut
import platosim.starquery       as sq
import platosim.referenceFrames as rf
from platosim.utilities    import errorcode
from platosim.matplotlibrc import latex
latex()


#==============================================================#
#                         BEGIN CLASS                          #
#==============================================================#


class Vizier(object):

    """Class to generate a Gaia DR3 catalogue.
    
    This class makes it straight forward to generate star catalogue
    directly from the Gaia DR3.
    """

    def __init__(self, args):

        # Global parameters

        self.field   = args.field
        self.maglim  = args.maglim
        self.plot    = args.plot
        self.verbose = args.verbose
        self.bright  = args.bright

        if not self.field:
            self.field = 'LOPS2'
            
        if not self.maglim:
            self.maglim = 21

        # Output directory
        if args.outdir:
            self.odir = Path(args.outdir).resolve()
        elif args.project:
            self.project = args.project
            self.odir = Path(os.getenv('PLATO_WORKDIR')) / self.project / "input"
        else:
            errorcode('error', 'An output destination is required! Use either -o or --project')

        # VERBOSITY (a.k.a log level) -> Identical to PlatoSim usage
        # verbose = 0: Cluster mode: Disabling print and warnings, and no log files are saved
        # verbose = 1: Default mode: Print details to bash but do not save log files
        # verbose = 3: Debug mode  : Print details to bash and saves all log files        
        self.verbose = args.verbose
        if self.verbose in [0, 1]:
            import warnings
            warnings.filterwarnings("ignore")        
        else:
            self.verbose = 3
            
        # Constants [deg]

        self.rGroup = 19.0                       # Max radius to query stars within a group
        self.cGroup = 9.2                        # Diagonal opening angle of groups
        self.aGroup = np.sqrt(self.cGroup**2/2)  # Horizontal/vertical opening angle of group

        self.r = 12.2                            # Radius of grid point search
        self.a = 17.2                            # Distance between equidistant grid points
        self.c = np.sqrt(2*self.a**2)            # Diagonal of grid
        
        # Get pointing of platform [rad] 
                
        self.alpha, self.delta, self.kappa = ut.getPointingField(self.field, unit='rad')

        # Find the camera group pointings [rad]
        
        raGroups, decGroups = rf.getCameraGroupCoordinates(self.alpha, self.delta, self.kappa)
        self.raGroups  = np.append(raGroups,  self.alpha)
        self.decGroups = np.append(decGroups, self.delta)

            


        
    def queryGaiaDR3(self):

        """Query function for Gaia DR3.

        This function query a circular area from the Gaia DR3 given a
        equatorial pointing of the PLATO spacecraft. It uses a grid of 
        9 circular regions to search for stars (in order not to exceed
        the RAM memory) and concatenates these grids into a final Gaia
        custom catalogue spanning each of the camera groups.

        Created on Mon Oct 17, 2022
        Authors: Juan Cabrera & Nicholas Jannsen
        Adapted from:
        https://www.cosmos.esa.int/web/gaia-users/archive/programmatic-access
        """
        
        # Make a grid in azimuth and tilt angles and find the sky coordinates [deg]

        azim = [-45,   0,  45,
                -90,   0,  90,
                225, 180, 135]
        tilt = [self.c, self.a, self.c,
                self.a,      0, self.a,
                self.c, self.a, self.c]

        # Get grid points [deg]
        
        x = np.zeros(len(azim))
        y = np.zeros(len(azim))
        for i in range(len(azim)):
            x[i], y[i] = rf.platformToTelescopePointingCoordinates(self.alpha,
                                                                   self.delta,
                                                                   self.kappa,
                                                                   np.deg2rad(azim[i]),
                                                                   np.deg2rad(tilt[i]))
        self.raGrid  = np.rad2deg(x)
        self.decGrid = np.rad2deg(y)
            
        # Plot grid used to fetch Gaia DR3 stars

        if self.plot:
            self.plotGrid()

        # Query stars within each grid point FOV

        if self.verbose > 0:
            print('Start Gaia DR3 query')

        filename = self.odir / f'starcat_GaiaDR3_{self.field}'
        N = len(self.raGrid)
        
        for i in tqdm(range(N), bar_format=ut.tqdmBar()):
            sq.gaiaRegionQuery(self.raGrid[i], self.decGrid[i], radius=self.r,
                               maglim=self.maglim, ofile=f'{filename}_grid{i+1}')

        # Combine all the data sets

        if self.verbose > 0:
            print('Load and combine all star found from grid search')

        dfs = [pd.read_feather(f'{filename}_grid{i+1}.ftr') for i in range(N)] 
        df  = pd.concat(dfs, axis=0, ignore_index=True)

        # Remove grid files

        for i in range(N):
            Path(f'{filename}_grid{i+1}.ftr').unlink()

        # Remove duplicates

        df = df.drop_duplicates(subset=['designation'])

        # Remove distance column

        df = df.drop(columns='dist')

        # Rename columns

        df = df.rename(columns={'designation': 'gaiaDR3',
                                'phot_g_mean_mag': 'mag',
                                'parallax': 'plx',
                                'parallax_error': 'plxe',
                                'teff_gspphot': 'teff',
                                'logg_gspphot': 'logg'})

        # Keep digit ID only

        df.gaiaDR3 = df.gaiaDR3.str[9:]

        # Add brightest stars not available in Gaia

        if self.bright:
            sirius  = {'gaiaDR3':'Sirius',  'ra':101.2871667, 'dec':-16.7161167, 'mag':-1.46}
            canopus = {'gaiaDR3':'Canopus', 'ra': 95.9879167, 'dec':-52.6956611, 'mag':-0.72}
            epscma  = {'gaiaDR3':'epsCMa',  'ra':104.6564583, 'dec':-28.9720861, 'mag': 1.50}
            df = df.append([sirius, canopus, epscma], ignore_index=True)

        # Keep only stars within the camera group FOV
        
        for i in range(5):

            if self.verbose > 0:
                print(f'Creating catalogue for camera group {i+1}')

            raStars  = np.deg2rad(df.ra.to_numpy())
            decStars = np.deg2rad(df.dec.to_numpy())

            # Calculate angular distance
            dOA = np.arccos(np.sin(self.decGroups[i]) * np.sin(decStars) +
                            np.cos(self.decGroups[i]) * np.cos(decStars) *
                            np.cos(self.raGroups[i] - raStars))
            df0 = df[np.rad2deg(dOA) < 19]
            
            # Select output filename

            ofile = f'{filename}_group{i+1}.ftr'

            # Save new catalogue
            
            df0.reset_index(inplace=True)
            df0.to_feather(ofile)

            
            




    def plotGrid(self):

        # Shorten names
        r, a, c = self.r, self.a, self.c
        rg, ag  = self.rGroup, self.aGroup

        # Plot the grid in cartesian coordinates
            
        fig, ax = plt.subplots(figsize=(7,7))

        ax.add_artist(plt.Circle(( 0, a), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle((-a, a), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle(( a, a), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle(( 0, 0), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle((-a, 0), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle(( a, 0), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle(( 0,-a), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle((-a,-a), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle(( a,-a), r, color='k', alpha=.2))
        ax.add_artist(plt.Circle(( ag, ag), rg, color='b', alpha=.2))
        ax.add_artist(plt.Circle(( ag,-ag), rg, color='b', alpha=.2))
        ax.add_artist(plt.Circle((-ag,-ag), rg, color='b', alpha=.2))
        ax.add_artist(plt.Circle((-ag, ag), rg, color='b', alpha=.2))

        # Settings
        ax.set_xlim(-30, 30)
        ax.set_ylim(-30, 30)
        ax.set_aspect('equal')
        plt.show()

        # Plot the grid points on the sky with the stars

        fig = plt.figure(figsize=(7,7))

        platform = SkyCoord(np.rad2deg(self.alpha), np.rad2deg(self.delta), frame='icrs', unit=u.deg)
        ax = plt.axes(projection='astro degrees zoom', center=platform,
                      radius='35 deg', rotate='180 deg')

        # Plot pointing of platform
        ax.plot(platform.ra.deg, platform.dec.deg, '*', c='k', mfc='magenta', ms=25,
                    transform=ax.get_transform('world'))

        # Plot the pointing of each camera group
        #if not self.fcam:
        colors = ['b', 'limegreen', 'yellow', 'r']
        for i, c in zip(range(4), colors):
            ax.plot(np.rad2deg(self.raGroups[i]), np.rad2deg(self.decGroups[i]),
                    'o', ms=13, c=c, mec='k', transform=ax.get_transform('world'))

        # Plot grid points
        grid = SkyCoord(self.raGrid, self.decGrid, frame='icrs', unit=u.deg)
        ax.plot(grid.ra.deg, grid.dec.deg, 'o', ms=10, c='k',mec='k',
                transform=ax.get_transform('world'))      

        # Settings
        ax.scalebar((0.05, 0.05), 10 * u.deg).label()
        ax.compass(0.95, 0.05, 0.1)
        ax.grid(color='gray')
        ax.set_xlabel('RA [deg]')
        ax.set_ylabel('Dec [deg]')
        plt.show()

            
#--------------------------------------------------------------#
#                PARSING COMMAND-LINE ARGUMENTS                #
#--------------------------------------------------------------#

parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', '\nPLATO FoV Catalogue Query'))

out_group = parser.add_argument_group('I/O PARAMETERS')
out_group.add_argument('-p', '--plot',    action='store_true',      help='Flag for plotting')
out_group.add_argument('-v', '--verbose', action='store_true',      help='Flag for verbosity')
out_group.add_argument('-o', '--outdir',  metavar='PATH', type=str, help='Output directory')
out_group.add_argument('--project',       metavar='NAME', type=str, help='PLATOnium project (overwrites --odir)')

que_group = parser.add_argument_group('QUERY PARAMETERS')
que_group.add_argument('--field',  metavar='NAME', type=str,   help='LOP (SPF, NPF, LOPS2, LOPN)')
que_group.add_argument('--maglim', metavar='MAG',  type=float, help='Maximum magnitude to query (default: 17)')
que_group.add_argument('--bright', action='store_true',        help='Flag add Yale bright stars catalogue')


# Initialize instance of class
args = parser.parse_args()
x = Vizier(args)
x.queryGaiaDR3()
print('')
