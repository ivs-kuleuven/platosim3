#!/usr/bin/env python3

# Built-in
import os
import warnings
import argparse
from pathlib import Path

# PlatoSim standard
import numpy as np
import pandas as pd

# PLATOnium standard
import star_shadow as ss

# PlatoSim functions
import platosim.utilities as ut
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
man_group.add_argument('idir',   type=str, help='Input directory containing star folders')
man_group.add_argument('odir',   type=str, help='Output directory to save analysis')
man_group.add_argument('--clean', action='store_true', help='Flag to remove camera data')
args = parser.parse_args()

# File paths
star = f'{args.starID}'.zfill(9)
idir = Path(args.idir).resolve() / star
odir = Path(args.odir).resolve()
odir_final = odir / 'final'
odir_table = odir / 'table'
odir_modes = odir / 'modes'
odir_final.mkdir(parents=True, exist_ok=True)
odir_table.mkdir(parents=True, exist_ok=True)
odir_modes.mkdir(parents=True, exist_ok=True)
os.system(f'chmod 755 {odir_final}')
os.system(f'chmod 755 {odir_table}')
os.system(f'chmod 755 {odir_modes}')
filename_final = f'lc_{star}'
filename_table = f'table_{star}'
filename_modes = f'modes_{star}'
filename_dat = odir_final / f'{filename_final}.dat'
filename_ftr = odir_final / f'{filename_final}.ftr'
filename_tab = odir_table / f'{filename_table}.ftr'
filename_mod = odir_modes / f'{filename_modes}.ftr'


# EXTRACT FINAL LIGHT CURVE

# Merge ligth curves
lcs = LightCurve(idir, 'multi')
lc = lcs.merge(flux_group_mean=True, binsize=600)

# Extract light curve
df = lc.data()
df = df.dropna(subset=['flux']).reset_index(drop=True)

# Save final light curve
ds = df
ds.to_feather(filename_ftr)
os.system(f'chmod 755 {filename_ftr}')

# Save simulation table
lc.stat_sim_table(filename_tab)


# STAR SHADOW ANALYSIS

# Prepare light curve for starshadow
dt = df
dt.time /= 86400. 
dt['flux_err'] = np.ones_like(dt.flux)
dt.to_csv(filename_dat, sep=' ', index=False, header=False)

# Perform prewhitening using STARSHADOW
ss.analyse_lc_from_file(str(filename_dat), save_dir=str(odir_final), stage='freq',
                        overwrite=True, verbose=False)

# Load file containing columns
folder_hdf5   = odir_final  / f'{filename_final}_analysis'
filename_hdf5 = folder_hdf5 / f'{filename_final}_analysis_2.hdf5' 
result = ss.utility.read_parameters_hdf5(filename_hdf5, verbose=False)

# Save data into feather file
mean = result['sin_mean']
err  = result['sin_err']
df = pd.DataFrame()
df['freq']      = mean[2]      # [c/d]
df['freq_err']  = err[2]       # [c/d]
df['ampl']      = mean[3]*1e6  # [ppm]
df['ampl_err']  = err[3]*1e6   # [ppm]
df['phase']     = mean[4]      # [rad]
df['phase_err'] = err[4]       # [rad]
df = df.sort_values('freq').reset_index(drop=True)
df.to_feather(filename_mod)
os.system(f'chmod 755 {filename_mod}')

# Remove starshadow files
filename_dat.unlink()
os.system(f'rm -r {folder_hdf5}')
if args.clean:
    os.system(f'rm -r {str(idir)}')
