#!/usr/bin/env python3

# Built-in
import os
import warnings
import argparse
from pathlib import Path

# PlatoSim standard
import numpy as np
import pandas as pd

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
man_group.add_argument('path',   type=str, help='Input directory containing star folders')

opt_group = parser.add_argument_group('OPTIONAL PARAMETERS')
opt_group.add_argument('-v', '--verbose', action='store_true', help='Flag print to bash')

args = parser.parse_args()

# File paths
if args.verbose: print(f'Processing star {args.starID}')
star = f'{args.starID}'.zfill(9)
path = Path(args.path).resolve()
idir = path / 'final'
odir = path / 'lightcurve'
odir.mkdir(parents=True, exist_ok=True)
filename_idir = idir / f'lc_{star}.ftr'
filename_odir = odir / f'lc_{star}.ftr'

# Construct light curve object
lc = LightCurve(filename_idir, mode='final')

# Sigma clipping
if args.verbose: print(f'Sigma-clipping')
df = lc.clip(model='wotan', flux_unit='ppt', sigma_lower=4, sigma_upper=4, replace=True)
df = df.dropna()

# Flux offset correction
if args.verbose: print(f'Corrrecting flux offset')
flux_offset = df.flux.median() - 1
df.flux    -= flux_offset        

# If requested save output file
if args.verbose: print('Saving new light curve')
df.reset_index(drop=True, inplace=True)
df.to_feather(filename_odir)
os.system(f'chmod 755 {filename_odir}')
