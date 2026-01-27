#!/usr/bin/env python3

# Built-in
import os
import glob
import warnings
import argparse
from pathlib import Path

# PlatoSim standard
import natsort
import numpy as np
import pandas as pd

# PlatoSim functions
import platosim.utilities as ut
from platosim.utilities  import errorcode
from platosim.lightcurve import LightCurve

# Ignore warnings
warnings.filterwarnings("ignore")

#--------------------------------------------------------------#
#                        START OF SCRIPT                       #
#--------------------------------------------------------------#

# Parse arguments
parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)

man_group = parser.add_argument_group('MANDATORY PARAMETERS')
man_group.add_argument('starID', type=int, help='Star ID')
man_group.add_argument('sdir',   type=str, help='Simulation directory with feather simulations')
man_group.add_argument('idir',   type=str, help='Input directory containing star folders')
man_group.add_argument('odir',   type=str, help='Output directory to save analysis')

opt_group = parser.add_argument_group('OPTIONAL PARAMETERS')
opt_group.add_argument('--detrend',    action='store_true', help='Perform auto-poly detrending')
opt_group.add_argument('--bin_size',   metavar='FLOAT', type=float, help='Time bin size [sec]')
opt_group.add_argument('--clip_sigma', metavar='FLOAT', type=float, help='Sigma-clip threshold')
opt_group.add_argument('--flux_err',   action='store_true', help='Calculate flux errors')
opt_group.add_argument('--gaps',       action='store_true', help='Use instrumentGAP.tab file')
opt_group.add_argument('-v', '--verbose', action='store_true', help='Flag print to bash')
opt_group.add_argument('-c', '--clean',   action='store_true', help='Flag to remove camera data')

args = parser.parse_args()

# User defined parameters
verbose    = args.verbose
bin_size   = args.bin_size
clip_sigma = args.clip_sigma
gaps       = args.gaps
flux_err   = args.flux_err

# Activate detrending
if args.detrend:
    detrend = 'poly'
else:
    detrend = False

# File paths
star = f'{args.starID}'.zfill(9)
sdir = Path(args.sdir).resolve() / star
idir = Path(args.idir).resolve()
odir = Path(args.odir).resolve()
odir_final = odir / 'lightcurve'
odir_table = odir / 'table'

errorcode('software', f'\nReducing star {star}')

# Create output directory
odir_final.mkdir(parents=True, exist_ok=True)
odir_table.mkdir(parents=True, exist_ok=True)
os.system(f'chmod 755 {odir_final}')
os.system(f'chmod 755 {odir_table}')

# Define filenames
filename_final = f'lc_{star}'
filename_table = f'table_{star}'
filename_ftr = odir_final / f'{filename_final}.ftr'
filename_tab = odir_table / f'{filename_table}.ftr'
filename_gap = idir / 'instrumentGAP.tab'

# EXTRACT FINAL LIGHT CURVE

# Construct light curve object
lcs = LightCurve(sdir, 'multi')

fig, ax = lcs.plot_multi(suffix='ftr', group=False, camera=False, quarter=False, 
                         flux_median=144, alpha=0.5, figsize=(9,5))

# Save simulation table
if args.verbose:
    print('Saving simulation table')
ds = lcs.stat_sim_table(filename_tab)

# Merge ligth curves
lc = lcs.merge(suffix='ftr',
               verbose=verbose,
               detrend=detrend,
               flux_group_mean=True,
               binsize=bin_size,
               clip_sigma=clip_sigma,
               flux_offset=True,
               flux_error=flux_err)

# Introducing data gaps
if gaps and filename_gap.is_file():
    if args.verbose:
        print('Introducing data gaps')
    df = lc.gaps(filename_gap, replace=True)
    df = df.dropna().reset_index(drop=True)
else:
    df = lc.data()

# Save feather with modes
df.to_feather(filename_ftr)
os.system(f'chmod 755 {filename_ftr}')

# Remove output and starshadow files
if args.clean:
    os.system(f'rm -r {str(idir)}')
