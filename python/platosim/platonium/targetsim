#!/usr/bin/env python3

"""
Usage examples:
  $ ./targetsim -o </path/to/outdir> <target name>

This script generates target and contaminant catalogues for a named target
and saves them to the specified output directory as .ftr files.

Notes on PLATO fields:
  Number of targets in PIC: 320,743
  Number of targets in SPF: 163,772
  Number of targets in NPF: 156,971
  Uses the de-reddened Gaia V magnitude obtained from Gaia
  photometry using a calibration technique being valid for:
  0.5 < Bp-Rp < 5.17

"""

import argparse
import datetime
import os
import urllib.request
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

# PlatoSim imports
import platosim.utilities as ut
from platosim.utilities import errorcode 
from platosim.starquery import gaiaQuery

# Start time tracking
tic = datetime.datetime.now()
warnings.filterwarnings('ignore')


# ==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                  #
# ==============================================================#


def load_targets(pic_tar):
    """
    Function to load target catalogue.
    Also needed for more efficient script speed
    """
    ID = pic_tar[:, 0].astype(float).astype(int)  # PICidDR1: PIC-ID-DR1 from Gaia DR2
    ra = pic_tar[:, 1].astype(np.float64)  # ra: ICRS RA
    dec = pic_tar[:, 2].astype(np.float64)  # decl: ICRS Dec
    mag = pic_tar[:, 3].astype(np.float64)  # gaiaV: De-reddened V mag from Gaia colour photometry
    sample = pic_tar[:, 4].astype(float).astype(int)  # sampleFlag: Bitmaskdefining PIC samples
    Teff = pic_tar[:, 5].astype(np.float64)  # teff: Stellar effective temperature [K]
    R = pic_tar[:, 6].astype(np.float64)  # radius: Stellar radius [R_sun]
    M = pic_tar[:, 7].astype(np.float64)  # mass: Stellar mass [M_sun]
    ncams = pic_tar[:, 8].astype(float).astype(int)  # nCameraObs: EOL number of cameras seeing the star
    field = pic_tar[:, 9].astype(str)
    return ID, ra, dec, mag, sample, Teff, R, M, ncams, field


def load_contaminants(pic_con):
    """
    Function to load contaminant catalogue.
    Also needed for more efficient script speed
    """
    ID_con = pic_con[:, 0].astype(float).astype(int)
    ra_con = pic_con[:, 1].astype(float)
    dec_con = pic_con[:, 2].astype(float)
    distance = pic_con[:, 3].astype(float)
    mag_con = pic_con[:, 4].astype(float)
    return ID_con, ra_con, dec_con, distance, mag_con


parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', '\nPIC targetsim'))
pic_group = parser.add_argument_group('PIC INFO (required)')

pic_group.add_argument('outdir', type=str, help='Output directory to save catalogue and plots')
pic_group.add_argument('target', type=str, help='Name of the target star to query around')

con_group = parser.add_argument_group('CONTAMINANTS (optional)')
con_group.add_argument('--dmag', type=int, metavar='MAG',
                       help='Delta magnitude inclusion threshold of target-to-contaminant (Default: 5 mag)')
con_group.add_argument('--dist', type=int, metavar='ARCSEC',
                       help='Radial distance threshold to include a contaminant [30, 45, 60] (Default: 60 arcsec)')

args = parser.parse_args()
outputDir = Path(args.outdir).resolve()
if not outputDir.exists():
    errorcode('error', f'Output directory does not exist: {outputDir}')
targetName = args.target

dmagConLimit = args.dmag
disConLimit = args.dist
inputDir = Path(os.getenv("PLATO_PROJECT_HOME") + "/inputfiles/data_picsim.")

inputFileTar = inputDir / 'PIC110target'
inputFileCon = inputDir / 'PIC110contaminant'
inputFileDict = inputDir / 'gaiaDR2-PIC110id'

# Prefix name of output files
outputPrefix = 'starcat'
outputPrefixTar = outputPrefix + '_targets.ftr'
outputPrefixCon = outputPrefix + '_contaminants.ftr'

outputFileTar = outputDir / outputPrefixTar
outputFileCon = outputDir / outputPrefixCon

# Visibility by N-cams: PIC gives EOL estimate, but we here assume that no camera are lost!
cams = {'6': 6, '12': 12, '18': 18, '24': 24}
cams_EOL = [[5, 11, 16, 17, 22], [6, 12, 18, 18, 24]]

# ==============================================================#
#                       LOAD PIC TARGETS                        #
# ==============================================================#

if dmagConLimit is None:
    dmagConLimit = 5  # Default
elif dmagConLimit > 22:
    errorcode('warning', 'The PIC (Gaia DR2) is only complete up to 22th magnitude!')
# Fetch the data targets and their contaminants from the PIC.


if disConLimit is None:
    disConLimit = 60  # Default equal to 4 pixels
elif disConLimit not in [30, 45, 60]:
    errorcode('error', 'Not a valid contaminant-to-target distance! Use either 30, 45, or 60 arcsec')

errorcode('module', '\nTarget star')
# Fetch PIC catalogue from FTP server

if not inputFileTar.with_suffix('.npy').exists():
    errorcode('message', 'Inauguration: Welcome to TargetSim!')
    print(f'Downloading PIC 1.1.0 catalogue...')
    url_tar = 'ftp://plato:miSotalP@ftp.ster.kuleuven.be/PIC110target.npy'
    url_con = 'ftp://plato:miSotalP@ftp.ster.kuleuven.be/PIC110contaminant.npy'
    with urllib.request.urlopen(url_tar) as response, open(inputFileTar.with_suffix('.npy'), "wb") as out_file:
        data = response.read()
        out_file.write(data)
    with urllib.request.urlopen(url_con) as response, open(inputFileCon.with_suffix('.npy'), "wb") as out_file:
        data = response.read()
        out_file.write(data)

if not inputFileDict.with_suffix('.npy').exists():
    errorcode('message', 'Inauguration: Welcome to TargetSim!')
    print(f'Downloading Gaia DR2 to PLATO PIC 1.1.0 ID Table...')
    url_dict = 'ftp://plato:miSotalP@ftp.ster.kuleuven.be/gaiaDR2-PIC110id.npy'
    with urllib.request.urlopen(url_dict) as response, open(inputFileDict.with_suffix('.npy'), "wb") as out_file:
        data = response.read()
        out_file.write(data)

pic_tar = np.load(inputFileTar.with_suffix('.npy'))

ID, ra, dec, mag, sample, Teff, R, M, ncams, field = load_targets(pic_tar)
# Create pandas data array
d = {'ID': ID, 'ra': ra, 'dec': dec, 'mag': mag, 'sample': sample,
     'Teff': Teff, 'R': R, 'M': M, 'ncams': ncams, 'field': field}

df = pd.DataFrame(d, columns=['ID', 'ra', 'dec', 'mag', 'sample',
                              'Teff', 'R', 'M', 'ncams', 'field'])
# Replace EOL N-CAM visibility with BOL
ncams = df['ncams'].replace(cams_EOL[0], cams_EOL[1])

# Save each N-CAM visibility into separate catalogues
df = df.assign(ncams=ncams)

# ==============================================================#
#                    Find Target in the PIC                     #
# ==============================================================#

# Remove any star with NaN for any parameter
df = df.dropna()
# Removing duplicate stars (if any)
df = df.drop_duplicates(subset=['ID'])
print(args.target)
# Get gaia ID of target star
try:
    gaiaID = gaiaQuery(targetName)
except LookupError:
    errorcode('error', 'Target star not found in Gaia DR2')
print(f'Gaia ID: {gaiaID}')
# Read gaia to pic mapping
gaia2pic = dict(np.load(inputFileDict.with_suffix('.npy')))
# Get pic ID of target star
if gaiaID not in gaia2pic:
    errorcode('error', "Target star found in Gaia DR2, but not PIC110. "
                       "\nThis could be because:"
                       "\n\t1. The star is not in PIC110 (see A&A 653, A98 (2021))"
                       "\n\t2. The star is in PIC110, but it was removed from TargetSim"
                       " due to its catalogue entry being incomplete")
picID = gaia2pic[gaiaID]
print(f'PIC ID:  {picID}')
pd.set_option('display.max_columns', None)
df = df[df['ID'] == picID]
print(f'RA, Dec: {df["ra"].iloc[0]:.6f}, {df["dec"].iloc[0]:.6f} ({"Northern" if df["field"].iloc[0] == "N" else "Southern"} Field)')
print(f'Vmag:    {df["mag"].iloc[0]:.2f}')
print(f'N-cams:  {df["ncams"].iloc[0]}')
print(f'Teff:    {df["Teff"].iloc[0]:.0f} K')
print(f'Radius:  {df["R"].iloc[0]:.2f} R_sun')
print(f'Mass:    {df["M"].iloc[0]:.2f} M_sun')

# FETCH STELLAR CONTAMINANTS FOR EACH TARGET STAR
errorcode('module', '\nPIC contaminants')

pic_con = np.load(inputFileCon.with_suffix('.npy'))
ID_con, ra_con, dec_con, distance, mag_con = load_contaminants(pic_con)

# Create pandas data array
dd = {'ID': ID_con, 'ra': ra_con, 'dec': dec_con, 'mag': mag_con, 'dis': distance}
dc = pd.DataFrame(dd, columns=['ID', 'ra', 'dec', 'mag', 'dis'])

print(f'Fetching contaminants within {disConLimit} arcsec from the target')

# Create empty data-frame to append to
dfc = pd.DataFrame(columns=['ID', 'ra', 'dec', 'mag', 'dis'])
numConPerTar = []

# Find corresponding contaminants to our target
dcc = dc[dc['ID'] == df['ID'].iloc[0]]

# We will only add stars within a user-defined distance of our target
dcc = dcc[dcc['dis'] < disConLimit]

# Convert to PLATO passband using host star Teff
Teff = np.full(len(dcc), df['Teff'].iloc[0])
teff = pd.DataFrame({'Teff': Teff}, columns=['Teff'])
dcc['mag'] = ut.passbandConversionV2P(dcc['mag'].values, teff['Teff'].values)

# Add only contaminants brigther than threshold
dcc = dcc[dcc['mag'] < df['mag'].iloc[0] + dmagConLimit]

# Append contaminants to contaminant list
dfc = pd.concat([dfc, dcc])
print("\nTable of contaminants:")
print(dfc)

# Store number of contaminants per target for statistics below
numConPerTar.append(len(dcc))

# Let's add the number of contaminants to the target catalog
df['ncon'] = np.array(numConPerTar).tolist()
print(f'\nNumber of contaminants: {df["ncon"].iloc[0]}')

# ==============================================================#
#                          SAVE OUTPUT                         #
# ==============================================================#

errorcode('module', '\nPrologue')

# We reset the index in order to save to feather
df = df.reset_index()
dfc = dfc.reset_index()
print(f'Saving file {outputFileTar}')
df.to_feather(outputFileTar)
print(f'Saving file {outputFileCon}')
dfc.to_feather(outputFileCon)

# Output file name of HPC data
outputVarFileName = f'{outputDir}/data_hpc_{outputPrefix}.txt'
print(f'Saving file {outputVarFileName}')
headerVar = 'PIC, M, R, Teff'
outputVarData = np.transpose([df['ID'].to_numpy(), df['M'].to_numpy(),
                              df['R'].to_numpy(), df['Teff'].to_numpy()])
np.savetxt(outputVarFileName, outputVarData, header=headerVar, comments='', fmt=['%i', '%0.3f', '%0.3f', '%i'],
           delimiter=',')

toc = datetime.datetime.now()
print('Total execution time: {0} [hh:mm:ss]'.format(toc - tic))
print('')
