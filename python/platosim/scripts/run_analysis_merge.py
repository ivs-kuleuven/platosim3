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
man_group.add_argument('idir',   type=str, help='Input directory containing star folders')
man_group.add_argument('odir',   type=str, help='Output directory to save analysis')
args = parser.parse_args()

# File paths
star = f'{args.starID}'.zfill(9)
idir = Path(args.idir).resolve() / star
odir = Path(args.odir).resolve()

#folders = natsort.natsorted(glob.glob(f'{idir}/*'))
#folder = Path(f)
#starID = folder.stem

errorcode('software', f'\nReducing star {star}')
    
odir_table = odir / 'table'
odir_final = odir / 'final'

odir_table.mkdir(parents=True, exist_ok=True)
odir_final.mkdir(parents=True, exist_ok=True)

os.system(f'chmod 755 {odir_table}')
os.system(f'chmod 755 {odir_final}')

filename_tab = odir_table / f'lc_{star}.tab'
filename_ftr = odir_final / f'lc_{star}.ftr'

# Merge ligth curves
lcs = LightCurve(idir, 'multi')

print('Unpacking data')
lcs.unpack()

lcs.stat_sim_table(ofile=filename_tab, clean=True)

lcs.merge(ofile=filename_ftr, suffix='ftr', binsize=1/6,
          flux_group_mean=True, flux_offset=True)

lcs.remove()
