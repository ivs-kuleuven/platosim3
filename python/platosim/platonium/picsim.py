#!/usr/bin/env python3

"""
This script is an integrated part of PlatoSim's toolkit PLATOnium.
This script take the PIC inputcatalog (csv format) and fetch the 
correct entries for target stars and their contaminants. Notice that
the catalog will only be saved if an ouput directory and prefix (-o)
are given as input. This script can also be used to re-plot a old 
output ascii catalog produced by this software. Note that the PIC
 contaminant catalog is complete up to 20 mag.

Notes on PLATO fields:
  Number of targets in PIC: 320,743
  Number of targets in SPF: 163,772
  Number of targets in NPF: 156,971
  Uses the de-reddened Gaia V magnitude obtained from Gaia
  photometry using a calibration technique being valid for:
  0.5 < Bp-Rp < 5.17

Note on Sample flag:
  Number of targets in P1:  15,094 (SPF:   6,817, NPF:   6,892)
  Number of targets in P2:   1,385 (SPF:     717, NPF:     668)
  Number of targets in P4:  33,032 (SPF:  16,866, NPF:  16,166)
  Number of targets in P5: 272,617 (SPF: 132,571, NPF: 140,046)
  Note that P2 is a bright sub-sample of P1.
  Note that P5 is on-board processed lightcures only.

Usage examples:
  $ picsim 100 P1 SPF -p
  $ picsim 1000 P5 NPF --camera 24
  $ picsim all all SPF --camera 6 --group 1 --mag 9.5-12.2 -o </path/to/outdir>
"""

# Python standard
import os
import shutil
import random
import warnings
import argparse
import datetime

# PlatoSim standard
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from prettytable import PrettyTable
from tqdm import tqdm

# PlatoSim functions
import platosim.plot       as pt
import platosim.utilities  as ut
import platosim.starquery as starQuery
from platosim.utilities  import errorcode, copyInputYAML
from platosim.simulation import Simulation

# Disable warnings
warnings.simplefilter("ignore")

# Start time tracking
tic = datetime.datetime.now()

                    
#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#


parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', '\nPIC of Destiny'))

pic_group = parser.add_argument_group('PIC INFO (required)')
pic_group.add_argument('stars',  type=str, help='Number of targets saved in catalog (Select "all" for all stars)')
pic_group.add_argument('sample', type=str, help='PIC samples [P1, P2, P4, P5, all]  (Select "all" for all samples)')
pic_group.add_argument('field',  type=str, help='PLATO fields [NPF, SPF]')

out_group = parser.add_argument_group('I/O PARAMETERS')
out_group.add_argument('-p', '--plot', action='store_true', help='Flag to plot target and contaminant catalog creation')
out_group.add_argument('-s', '--save', action='store_true', help='Flag to save target and contaminant catalog into a PlatoSim like starcatalog format')
out_group.add_argument('-t', '--targ', action='store_true', help='Flag to only show catalogue illustration (mainly testing)')
out_group.add_argument('-o', '--outdir', type=str, metavar='STR', help='Output directory to save catalogue and plots')
out_group.add_argument('--project',     type=str, metavar='NAME', help='Name project folder within $PLATO_WORKDIR')
out_group.add_argument('-i', '--incat',  type=str, nargs='*',     help='PIC input target- and contaminant files (Format: csv, txt, ftr)')
out_group.add_argument('-u', '--unique', type=str, metavar='STR', help='Parse old picsim target catalogue to select new unique stars (Format: ftr)')

que_group = parser.add_argument_group('QUERY OPTIONS (optional)')
que_group.add_argument('--query', type=str, metavar='NAME',   help='Simbad target for query using Gaia DR2')
que_group.add_argument('--dmag',  type=int, metavar='MAG',    help='Delta magnitude inclusion threhold of target-to-contaminatn (Default: 5 mag)')
que_group.add_argument('--dist',  type=int, metavar='ARCSEC', help='Radial distance to fecth contaminants [30, 45, 60] (Default: 45 arcsec)')

obs_group = parser.add_argument_group('OBSERVABLES (optional)')
obs_group.add_argument('--ncams', metavar='NUM',   type=str, help='N-CAM visibility [6, 12, 18, 24]')
obs_group.add_argument('--group', metavar='NUM',   type=int, help='Camera-group visibility [1, 2, 3, 4]')
obs_group.add_argument('--mag',   metavar='RANGE', type=str, help='P magnitude range of PIC targets (float or dash-range)')
obs_group.add_argument('--spec',  metavar='TYPE',  type=str, help='Spectral type [F, G, K]')
obs_group.add_argument('--lum',   metavar='CLASS', type=str, help='Luminosity class [V, IV]')
obs_group.add_argument('--ntime', metavar='NUM',   type=int, help='Number of individual timeseries to be generated')

args = parser.parse_args()

numTargets    = args.stars
starSample    = args.sample
platoField    = args.field

plot          = args.plot
subfield      = args.save
onlyTargets   = args.targ
outputDir     = args.outdir
project       = args.project
inputFiles    = args.incat
oldCatalogue  = args.unique

targetQuery   = args.query
dmagConLimit  = args.dmag
disConLimit   = args.dist

ncamVisible   = args.ncams
camGroup      = args.group
magRange      = args.mag
specType      = args.spec
lumClass      = args.lum
numTimeseries = args.ntime

# Add latex font if catalogue is saved
if outputDir is None:
    from platosim.matplotlibrc import setup
    setup()
else:
    from platosim.matplotlibrc import latex
    latex()

# Activate plot if `-t` is parsed
if onlyTargets: plot = True

# Handle size of samples

numTargetsTotalInPIC = 320743  # NOTE Hard-code value from PIC1.1.0

if numTargets[0].isdigit():
    numTargets = int(numTargets)
elif numTargets == 'all':
    if numTimeseries is not None:
        nTargets   = int(numTimeseries/6.) + 1
        numTargets = int(numTimeseries/6.) + 1
    else:
        nTargets = 'all'
        numTargets = numTargetsTotalInPIC

# Detect input file format
fileFormat = '.ftr'

# Default is to use the latest PIC catalog saved
# Else several old catalogs can be parsed for replotting
# E.g. use: -intar starcat**targets.txt --incon starcat**contaminants.txt
inputDir     = Path(os.getenv("PLATO_PROJECT_HOME")) / 'inputfiles/data_picsim'
inputFileTar = inputDir / 'PIC110target.npy'
inputFileCon = inputDir / 'PIC110contaminant.npy'

# Prefix name of output files
outputPrefix = 'starcat'
if args.sample: outputPrefix += '_' + args.sample
if args.field:  outputPrefix += '_' + args.field
if args.group:  outputPrefix += '_Group' + str(args.group)
if args.ncams:  outputPrefix += '_Ncam' + args.ncams
if inputFiles is not None: outputPrefix += '_NewCat'
outputPrefixTar = outputPrefix + '_targets'      + fileFormat
outputPrefixCon = outputPrefix + '_contaminants' + fileFormat

# Name space for output files

if outputDir is not None or project is not None:
    # Resolve absolute output directory
    if outputDir is not None:
        outputDir = Path(outputDir).resolve()
    elif project is not None:
        outputDir = Path(os.getenv('PLATO_WORKDIR')) / project / 'input'
    # Set file names
    if subfield:
        outputFileFOV = (outputDir / outputPrefix).joinwith_suffix(fileFormat)
    else:
        outputFileTar = outputDir / outputPrefixTar
        outputFileCon = outputDir / outputPrefixCon

        
#-----------------------------------------------------------------
# TODO implement properly such dmag and dist can be requested too
if targetQuery:
    df = starQuery(targetQuery)
    print(f'Catalogue for {targetQuery}:')
    print(df)
    
    if outputDir is not None:
        df.to_csv(outputDir / f"starcat_{targetQuery}.txt", sep=' ', header=False)    
    exit()
#-----------------------------------------------------------------

        
# Select PIC sample
bitmask = {'P1': 1,
           'P2': 3,  # Notice 2 is stated in documentation but 3 is correct..
           'P4': 4,
           'P5': 8}

if starSample == 'all':
    samplePIC = None
else:
    try: bitmask[starSample]
    except KeyError: errorcode('error', 'Not a valid PIC sample! Use either P1, P2, P4, or P5')
    else: samplePIC = bitmask[starSample]

# Select magnitude range or convert string to actual list range

if magRange is None: magRange = [0, 20]  # Select all stars
else: magRange = ut.convertMagnitudeRange(magRange)

# Check for PLATO pointing field: ICRS-Ra-Dec (deg) and Galactic-long.-lat. (deg)

alpha, delta, kappa = ut.getPointingField(platoField, unit='deg')

if platoField in ['NPF']:
    pointingField = 'N'
else:
    pointingField = 'S'
    
# Visibility by N-cams:
# NOTE PIC gives EOL estimate but we here assume that no camera are lost!

cams     = {'6':  6, '12': 12, '18': 18, '24': 24}
cams_EOL = [[5, 11, 16, 17, 22], [6, 12, 18, 18, 24]]

if ncamVisible is not None:
    try: cams[ncamVisible]
    except KeyError: errorcode('Not valid number of N-CAMs! Use either 6, 12, 18, or 24')
    else: pass

# Check cam-group entry

if camGroup is not None:
    if camGroup not in [1, 2, 3, 4]:
        errorcode('error', 'Not valid cam-group number! Use either 1, 2, 3, or 4')
    else: pass

# Contaminant limits

if dmagConLimit is None:
    dmagConLimit = 5  # Default
elif dmagConLimit > 22:
    errorcode('warning', 'The PIC (Gaia DR2) is only complete up to 22th magnitude!')

if disConLimit is None:
    disConLimit = 45  # Default equal to 2 pixels
elif disConLimit not in [30, 45, 60]:
    errorcode('error', 'Not a valid contaminant-to-target distance! ' +
              'Use either 30, 45, or 60 arcsec')

#==============================================================#
#                       LOAD PIC TARGETS                       #
#==============================================================#
# Fetch the datatargets and their contaminants from the PIC.

errorcode('module', '\nPIC targets')

def load_targets(pic_tar):
    """
    Function to load target catalogue. 
    Also needed for more efficient script speed
    """
    ID     = pic_tar[:,0].astype(float).astype(int)  # PICidDR1: PIC-ID-DR1 from Gaia DR2
    ra     = pic_tar[:,1].astype(np.float64)  # ra: ICRS RA
    dec    = pic_tar[:,2].astype(np.float64)  # decl: ICRS Dec
    mag    = pic_tar[:,3].astype(np.float64)  # gaiaV: De-reddened V mag from Gaia colour photometry
    sample = pic_tar[:,4].astype(float).astype(int)  # sampleFlag: Bitmaskdefining PIC samples
    Teff   = pic_tar[:,5].astype(np.float64)  # teff: Stellar effective temperature [K]
    R      = pic_tar[:,6].astype(np.float64)  # radius: Stellar radius [R_sun]
    M      = pic_tar[:,7].astype(np.float64)  # mass: Stellar mass [M_sun]
    ncams  = pic_tar[:,8].astype(float).astype(int) # nCameraObs: EOL number of cameras seeing the star
    field  = pic_tar[:,9].astype(str)
    return ID, ra, dec, mag, sample, Teff, R, M, ncams, field

# Fetch PIC catalogue from FTP server
try:
    pic_tar = np.load(inputFileTar.with_suffix('.npy'))
except:
    inuaguration = True
    errorcode('message', 'Inuaguration: Welcome to the PIC of Destiny!')
    print(f'Downloading PIC 1.1.0 catalogue..')
    ut.downloadFromFTP(filename=inputFileTar.name, outputDir=inputDir, server='plato')
    ut.downloadFromFTP(filename=inputFileCon.name, outputDir=inputDir, server='plato')
    
# Load old PIC catalog for re-plot an old catalogue
if inputFiles is not None:

    # Params
    inuaguration = False
    N = int(len(inputFiles)/2.)

    # Check if one or more catalogs are parsed
    # NOTE using *.ftr will order in alphabetic order
    if len(inputFiles) > 2:
        inputFileTarget = inputFiles[1::2]
        inputFileContam = inputFiles[0::2]
    else:
        inputFileTarget = [inputFiles[1]]
        inputFileContam = [inputFiles[0]]

    # Append each catalog if several are parsed
    for i in range(N):
        if i == 0:
            DF  = pd.read_feather(inputFileTarget[i])
            DFC = pd.read_feather(inputFileContam[i])
        else:
            DF  = DF.append(pd.read_feather(inputFileTarget[i]))
            DFC = DFC.append(pd.read_feather(inputFileContam[i]))

# If it is the first run of this software on a CSV format catalogue we save
# catalog to a binary file. This dramitically speed-up any future data reads
try:
    pic_tar = np.load(inputFileTar.with_suffix('.npy'))
except:
    # Load ascii catalogue
    data = np.genfromtxt(inputFileTar.with_suffix('.csv'), delimiter=',',
                        usecols=[0, 3, 5, 53, 55, 56, 58, 60, 71])
    ID     = data[:,0]
    ra     = data[:,1]
    dec    = data[:,2]
    mag    = data[:,3]
    sample = data[:,4]
    Teff   = data[:,5]
    R      = data[:,6]
    M      = data[:,7]
    ncams  = data[:,8]
    # String field needs to be loaded seperately: PLATO field: N=North, S=South
    field = np.loadtxt(inputFileTar.with_suffix('.csv'), delimiter=',',
                       usecols=[68], dtype=str)
    # Load catalogue
    pic_tar = np.transpose([ID, ra, dec, mag, sample, Teff, R, M, ncams, field])
    ID, ra, dec, mag, sample, Teff, R, M, ncams, field = load_targets(pic_tar)
    # TODO for PIC2.0.0
    #----------------------------------------
    # Convert V Jonhson-Cousin to P passband
    mag = ut.passbandConversionV2P(mag, Teff)
    #----------------------------------------
    # Save dataset to binary file:
    errorcode('message', 'Inuaguration! Welcome to the PIC of Destiny!')
    print('PIC catalogue is saved to binary format for optimised performance!')
    inuaguration = True
    output_file_tar = inputFileTar.with_suffix('.npy')
    output_data_tar = np.transpose([ID, ra, dec, mag, sample, Teff, R, M, ncams, field])
    np.save(output_file_tar, output_data_tar)
else:
    # If saved open binary table
    ID, ra, dec, mag, sample, Teff, R, M, ncams, field = load_targets(pic_tar)

# Create pandas data array
d = {'ID': ID,  'ra': ra, 'dec': dec, 'mag': mag, 'sample': sample,
     'Teff': Teff, 'R': R, 'M': M, 'ncams': ncams, 'field': field}
df = pd.DataFrame(d, columns=['ID', 'ra', 'dec', 'mag', 'sample',
                              'Teff', 'R', 'M', 'ncams', 'field'])

# Replace EOL N-CAM visibility with BOL
ncams = df['ncams'].replace(cams_EOL[0], cams_EOL[1])
df = df.assign(ncams = ncams)

# Save each N-CAM visibility into seperate catalogues
if inputFileTar.with_suffix('.npy').is_file():

    # Save SPF for each group
    df_SPF = df[df['field'] == 'S']
    df_SPF06 = df_SPF[df_SPF['ncams'] == 6]
    df_SPF12 = df_SPF[df_SPF['ncams'] == 12]
    df_SPF18 = df_SPF[df_SPF['ncams'] == 18]
    df_SPF24 = df_SPF[df_SPF['ncams'] == 24]

    df_SPF06 = df_SPF06[df_SPF06.columns[1:3]].to_numpy()
    df_SPF12 = df_SPF12[df_SPF12.columns[1:3]].to_numpy()
    df_SPF18 = df_SPF18[df_SPF18.columns[1:3]].to_numpy()
    df_SPF24 = df_SPF24[df_SPF24.columns[1:3]].to_numpy()

    np.save(inputDir / 'SPF-NCAM06.npy', df_SPF06)
    np.save(inputDir / 'SPF-NCAM12.npy', df_SPF12)
    np.save(inputDir / 'SPF-NCAM18.npy', df_SPF18)
    np.save(inputDir / 'SPF-NCAM24.npy', df_SPF24)

    # Save NPF for each group
    df_NPF = df[df['field'] == 'N']
    df_NPF06 = df_NPF[df_NPF['ncams'] == 6]
    df_NPF12 = df_NPF[df_NPF['ncams'] == 12]
    df_NPF18 = df_NPF[df_NPF['ncams'] == 18]
    df_NPF24 = df_NPF[df_NPF['ncams'] == 24]

    df_NPF06 = df_NPF06[df_NPF06.columns[1:3]].to_numpy()
    df_NPF12 = df_NPF12[df_NPF12.columns[1:3]].to_numpy()
    df_NPF18 = df_NPF18[df_NPF18.columns[1:3]].to_numpy()
    df_NPF24 = df_NPF24[df_NPF24.columns[1:3]].to_numpy()

    np.save(inputDir / 'NPF-NCAM06.npy', df_NPF06)
    np.save(inputDir / 'NPF-NCAM12.npy', df_NPF12)
    np.save(inputDir / 'NPF-NCAM18.npy', df_NPF18)
    np.save(inputDir / 'NPF-NCAM24.npy', df_NPF24)

#==============================================================#
#                  SELECT PIC TARGETS ON INPUT                 #
#==============================================================#

if fileFormat != '.txt':

    # CHECK PARAMETERS FROM THE PIC ITSELF

    # Remove any star with NaN for any parameter
    df = df.dropna()

    # Removing dublicate stars (if any)
    df = df.drop_duplicates(subset=['ID'])

    # Check pointing field
    if platoField is not None:
        df = df[df['field'] == pointingField]

    # Check stellar sample
    if samplePIC is not None:
        df = df[df['sample'] == samplePIC]

    # Check parsing of old catalogue to select unique new targets
    if oldCatalogue:
        try:
            df_old = pd.read_feather(oldCatalogue)
        except:
            errorcode('error', 'Old picsim target catalogue do not exist!')
        else:
            cond = df['ID'].isin(df_old['ID']) # Boolen
            df   = df.drop(df[cond].index)

    # Seperate dwarf (MS) and sub-gaint (post MS) stars
    def ms_limit(Teff):
        return 1 + 2 * 1e-7 * (Teff - 4000)**2
    ds = df[df['R'] < ms_limit(df['Teff'])]
    sg = df[df['R'] > ms_limit(df['Teff'])]
    # Seperate dwarf spectral types
    df_dK = ds[ ds['Teff'] < 5300]
    df_dG = ds[(ds['Teff'] > 5300) & (ds['Teff'] < 5900)]
    df_dF = ds[ ds['Teff'] > 5900]
    # Seperate sub-giants spectral types
    df_sgK = sg[ sg['Teff'] < 5300]
    df_sgG = sg[(sg['Teff'] > 5300) & (sg['Teff'] < 5900)]
    df_sgF = sg[ sg['Teff'] > 5900]

    # Check evolutionary stage 
    if lumClass is not None:
        if lumClass == 'V':
            df = ds
        elif lumClass == 'IV':
            df = sg
        else:
            errorcode('error', 'Not valid evolutionary stage for PIC stars!')

    # Check spectral type
    if specType is not None:
        if specType == 'K':
            df = df_dK if lumClass == 'V' else (df_sgK if lumClass == 'IV' else pd.concat([df_dK, df_sgK]))
        elif specType == 'G':
            df = df_dG if lumClass == 'V' else (df_sgG if lumClass == 'IV' else pd.concat([df_dG, df_sgG]))
        elif specType == 'F':
            df = df_dF if lumClass == 'V' else (df_sgF if lumClass == 'IV' else pd.concat([df_dF, df_sgF]))
        else:
            errorcode('error', 'Not valid spectral type for PIC stars!')

    # Select stars depending on camera visibility
    # NOTE this will be done properly later but here it is done to catch
    #      the failure when the stars requested are larger than catalogue

    if ncamVisible is not None:
        df = df[df['ncams'] == int(ncamVisible)]

    # Allow to select only one cam-group
    if camGroup:
        sim = Simulation('picsim')
        dex = sim.getStarsWithinCameraGroup(df['ra'].to_numpy(), df['dec'].to_numpy(),
                                            alpha, delta, kappa, camGroup)
        
        # Select max number of stars if too many is requested
        max_stars = np.sum(dex[0])
        if numTargets > max_stars:
            numTargets = max_stars
            errorcode('warning', f'Too many stars requested - using max of {max_stars} instead!')
        df = df[dex[0]]

    # Check P passband magnitude range
    if magRange is not None:
        df = df[(df['mag'] > magRange[0]) & (df['mag'] < magRange[1])]


    # Check if too many stars are selected
    if len(df['ID']) > numTargetsTotalInPIC:
        errorcode('error', 'More stars selected than available in the PIC! See -h for sample notes')
    elif numTargets > len(df['ID']) and samplePIC is None and args.stars != 'all':
        errorcode('error', 'More stars selected than available from chosen PIC sample and/or PLATO field! See -h for sample notes')

    # CHECK OBSERVABILITY BY EACH CAMERA

    # Since the PIC are bigger than actual FOV (18.89 deg radius) we need to
    # Loop over each camera-group and avoid stars that are close to the FOV edge
    # NOTE the following also checks if the subfield can be placed on a CCD
    camGroup = False
    if camGroup:

        # CONFIGURE SPACECRAFT

        # # Set the telescope group ID
        # sim["Telescope/GroupID"] = camGroup        

        # # Solar panel orientation: 0, 90, 180, and 270 degrees for Q1, Q2, Q3, and Q4
        # solarPanelOrientation = sim["Platform/SolarPanelOrientation"] = math.fmod(self.quarter * 90., 360.)
        # self.solarPanelOrientation = np.deg2rad(float(solarPanelOrientation))

        # # Select PLATO pointing field from inputfile (e.g. SPF, NPF, etc.)
        # pointingField  = sim["ObservingParameters/StarCatalogFile"]
        # raPlatformDeg  = sim["ObservingParameters/RApointing"]  = alpha
        # decPlatformDeg = sim["ObservingParameters/DecPointing"] = delta

        # # Camera-group Alt (tilt) and Az
        # tiltTelescope    = sim["CameraGroups/TiltAngle"][self.group-1]
        # azimuthTelescope = sim["CameraGroups/AzimuthAngle"][self.group-1]

        # # Include spacecraft Pointing Repeatability Error (PRE) between consecutive quarters
        # # NOTE: Will only be included if the file "PRE.txt" is available in the input folder 
        # if os.path.exists(self.inputDir + '/PRE.txt'):
        #     PRE = np.loadtxt(self.inputDir + '/PRE.txt')
        #     dex = np.where(PRE[:,0] == self.quarter)[0]
        #     raPlatformDeg         += PRE[dex, 1][0]
        #     decPlatformDeg        += PRE[dex, 2][0]
        #     solarPanelOrientation += np.deg2rad(PRE[dex, 3][0])

        # # Include Absolute Pointing Error (APE) due to camera misalignments
        # # NOTE: Will only be included if the file "APE.txt" is available in the input folder
        # if os.path.exists(self.inputDir + '/APE.txt'):
        #     APE = np.loadtxt(self.inputDir + '/APE.txt')
        #     dex = (self.group - 1) * 6 + self.camera - 1
        #     tiltTelescope    += APE[dex, 0]
        #     azimuthTelescope += APE[dex, 1]

        print('Checking star visibility for each camera group..')
        dexGroup1, d1 = getStarsWithinCamGroup(1, raPF, decPF, df['ra'].to_numpy(), df['dec'].to_numpy())
        dexGroup2, d2 = getStarsWithinCamGroup(2, raPF, decPF, df['ra'].to_numpy(), df['dec'].to_numpy())
        dexGroup3, d3 = getStarsWithinCamGroup(3, raPF, decPF, df['ra'].to_numpy(), df['dec'].to_numpy())
        dexGroup4, d4 = getStarsWithinCamGroup(4, raPF, decPF, df['ra'].to_numpy(), df['dec'].to_numpy())

        # NOTE check how many stars that are observable in each group
        #-----------------
        # starGroup = np.sum([dexGroup1*1, dexGroup2*1, dexGroup3*1, dexGroup4*1], axis=1)
        # print(starGroup)
        #-----------------

        # NOTE the following is a small table for checking radial distance to OA
        # A ghost within 8 deg from the OA will produce a point-like ghost
        #----------------
        # Pmag = ut.passbandConversionV2P(mag, Teff)
        # t2 = PrettyTable(['PIC', 'V', 'P', 'd_G1 (deg)', 'd_G2 (deg)', 'd_G3 (deg)', 'd_G4 (deg)'])
        # for i in range(len(ID)):
        #     t2.add_row([int(ID[i]), '{:.2f}'.format(mag[i]), '{:.2f}'.format(Pmag[i]), '{:.2f}'.format(d1[i]),
        #                 '{:.2f}'.format(d2[i]), '{:.2f}'.format(d3[i]), '{:.2f}'.format(d4[i])])
        # print('\nUser have requested data table of all targets')
        # print(t2); exit()
        #----------------

        # Find the actual N-Cam visibility
        ncams = (dexGroup1*1 + dexGroup2*1 + dexGroup3*1 + dexGroup4*1) * 6

        # Remove stars that are not observable at all
        dex = dexGroup1 + dexGroup2 + dexGroup3 + dexGroup4
        ID    = ID[dex].astype(int)
        ra    = ra[dex]
        dec   = dec[dex]
        mag   = mag[dex]
        Teff  = Teff[dex]
        R     = R[dex]
        M     = M[dex]
        ncams = ncams[dex]

        # Again choose stars from N-Cam visibility
        if ncamVisible:
            dex = ncams == int(ncamVisible)
            ID    = ID[dex].astype(int)
            ra    = ra[dex]
            dec   = dec[dex]
            mag   = mag[dex]
            Teff  = Teff[dex]
            R     = R[dex]
            M     = M[dex]
            ncams = ncams[dex]

        # Allow to select only one cam-group
        if camGroup:
            dex = starsWithinCamGroup(camGroup, raPF, decPF, ra, dec)
            ID    = ID[dex].astype(int)
            ra    = ra[dex]
            dec   = dec[dex]
            mag   = mag[dex]
            Teff  = Teff[dex]
            R     = R[dex]
            M     = M[dex]
            ncams = ncams[dex]

    # We need this here in order to run the above
    if inputFiles is not None:
        df = DF

    # Only after the above cuts can we select randomly a sub-sample (if applicable)
    # Use random.choices(data, weights=some weight, k=numTargets) if a weigthing is needed
    # NOTE args.stars being a string is needed to select "all" in one of the samples too
    # NOTE we need to shuffle the sample if a specific number of timeseries are needed
    if args.stars != 'all' or numTimeseries is not None:
        df = df.sample(n = numTargets)

    # Redefine the number fo targets
    numTargets = len(df['ID'])

    # Randomly choose sample until the number of timeseries are met
    if numTimeseries:
        numImagettes = 0
        for i in range(numTargets):
            # Add camera visibilities
            numImagettes += df['ncams'].to_numpy()[i]
            # The data it cut by the index that reached the limit
            if numImagettes >= numTimeseries:
                print(i)
                df = df.drop(df.index[:i])
                break

    # MAKE AN OVERVIEW TABLE

    cameras = [6, 12, 18, 24]
    countStars = [df[df['ncams'] == 6].count()[0],  df[df['ncams'] == 12].count()[0],
                  df[df['ncams'] == 18].count()[0], df[df['ncams'] == 24].count()[0]]
    countImagettes = [countStars[i] * cameras[i] for i in range(4)]

    # Last collection of stars
    if fileFormat == '.txt': numTargets = len(ncams)
    else: numTargets = np.sum(countStars)

    # Make table
    t = PrettyTable(['N-cams', 'Targets', 'Time series'])
    for i in range(4):
        t.add_row([str(cameras[i]), countStars[i], countImagettes[i]])
    tt = 'User catalog contains {0} stars of {1} time series drawn from:'.format(numTargets, np.sum(countImagettes))
    print(tt)
    print(t)

# Correct for missing masses using mass-luminosity relation (0.43 < M/Msun < 2)
# https://en.wikipedia.org/wiki/Mass%E2%80%93luminosity_relation
# NOTE Perhaps better in PIC2.0.0?
#---------------------------------
# for i in range(len(ID)):
#     if np.isnan(M[i]):
#         M[i] = R[i]**(1/2) * Teff[i]/5778.
        #print(M[i], R[i])
#---------------------------------

# PLOT STARS

# Plot targets in aitoff Galactic sky projection
fig0, ax = pt.drawStarsInSkyAitoff(df['ra'], df['dec'], df['mag'])
if plot: plt.show()

# Save if both pointing fields are given
if platoField is None:
    fig0[0].savefig(outputDir+f'plot_{outputPrefix}_allsky.png', bbox_inches='tight', dpi=200)
    exit()

# Plot Zoom-in on FOV
if len(df['ID']) > 200: plotmag = None
else: plotmag = df['mag'].to_numpy()
title = f'{platoField} {starSample} sample'
fig1, ax = pt.plotPlatoFOV(platoField, df['ra'].to_numpy(), df['dec'].to_numpy(), magStars=plotmag, showGroups=True, showLegend=True, title=title)
if plot: plt.show()

# Plot sample distribution in Teff vs. Radius
title = f'{platoField} {starSample} sample'
fig2, ax = pt.plotTeffvsRadius(ds, df_dK, df_dG, df_dF,
                               sg, df_sgK, df_sgG, df_sgF,
                               df, ms_limit, title)
if plot: plt.show()

# Stop here when flag -t
if onlyTargets: exit()

# NOTE use this to save target catalogue only!
#---------------------------------
# outputTarData = np.transpose([ra, dec, mag])
# np.savetxt(outputFileTar, outputTarData, fmt=['%f', '%f', '%f'])
# exit()
#---------------------------------

#==============================================================#
#                       PIC CONTAMINANTS                       #
#==============================================================#

# LOAD PIC CONTAMINANTS
# We here use a distance of 45 arcsec provided by the PIC catalog to
# include contaminants beloging to a parent target star. Notice that there
# are more than 8 million contaminants stars in the original catalog,
# hence, this step can take some time if the large catalog are to be used.
errorcode('module', '\nPIC contaminants')

def load_contaminants(pic_con):
    """
    Function to load target catalogue. 
    Also needed for more efficient script speed
    """
    ID_con   = pic_con[:,0].astype(float).astype(int)
    ra_con   = pic_con[:,1].astype(float)
    dec_con  = pic_con[:,2].astype(float)
    distance = pic_con[:,3].astype(float)
    mag_con  = pic_con[:,4].astype(float)
    return ID_con, ra_con, dec_con, distance, mag_con

# We here allow to load in an already existing txt file for the creation of new plots
if inputFiles is not None and fileFormat == '.ftr':
    dfc = DFC

else:
    # If it is the first run of this software we save catalog to a binary file
    # This dramitically speed-up any future data reads
    try:
        pic_con = np.load(inputFileCon.with_suffix('.npy'))
    except:
        errorcode('message', 'Inuaguration! Welcome to the PIC of Destiny!')
        print('PIC catalogue is saved to binary format for optimised performance!')
        print('Be patient, this will take around 1 min the first time only!')
        # Load ascii catalogue
        pic_con = np.loadtxt(inputFileCon.with_suffix('.csv'), delimiter=',', skiprows=1, usecols=[2, 6, 8, 5, 19])
        ID_con, ra_con, dec_con, distance, mag_con = load_contaminants(pic_con)
        # TODO for PIC2.0.0
        #----------------------------------------
        # Convert to PLATO passband
        #mag_con = ut.passbandConversionV2P(mag_con, Teff_con)
        #----------------------------------------
        # Save dataset to binary file:
        output_file_con = inputFileCon.with_suffix('.npy')
        output_data_con = np.transpose([ID_con, ra_con, dec_con, distance, mag_con])
        np.save(output_file_con, output_data_con)
    else:
        # If saved open binary table
        ID_con, ra_con, dec_con, distance, mag_con = load_contaminants(pic_con)

    # Create pandas data array
    dd = {'ID': ID_con,  'ra': ra_con, 'dec': dec_con, 'mag': mag_con, 'dis': distance}
    dc = pd.DataFrame(dd, columns=['ID', 'ra', 'dec', 'mag', 'dis'])

    # FETCH STELLAR CONTAMINANTS FOR EACH TARGET STAR

    print('Fetching contaminants within {0} arcsec from each target'.format(disConLimit))

    # Create empty data-frame to append to 
    dfc = pd.DataFrame(columns=['ID', 'ra', 'dec', 'mag', 'dis'])
    numConPerTar = []

    for starNo in tqdm(range(numTargets), bar_format=ut.tqdmBar()):

        # Find corresponding contaminants to our target
        dcc = dc[dc['ID'] == df['ID'].iloc[starNo]]

        # We will only add stars within a user-defined distance of our target
        dcc = dcc[dcc['dis'] < disConLimit]

        # Convert to PLATO passband using host star Teff
        # NOTE we should change this in the future!
        Teff = np.full(len(dcc), df['Teff'].iloc[starNo])
        teff = pd.DataFrame({'Teff': Teff}, columns=['Teff'])
        dcc['mag'] = ut.passbandConversionV2P(dcc['mag'].values, teff['Teff'].values)

        # Add only contaminants brigther than threshold
        dcc = dcc[dcc['mag'] < df['mag'].iloc[starNo] + dmagConLimit]

        # Append contaminats to contaminant list
        dfc = dfc.append(dcc)

        # Store number of contaminats per target for statistics below
        numConPerTar.append(len(dcc))

    # Let's add the number of contaminats to the target catalog
    df['ncon'] = np.array(numConPerTar).tolist()


# OVERVIEW TABLE OF CATALOGUE
#-----------------------------
# gaia = pic_tar[:,1]  #Gaia ID identifier col 2
# t2 = PrettyTable(['PIC',  'RA (deg)', 'Dec (deg)', 'P', 'N-Cams', 'Con60'])
# for i in range(len(ID)):
#     t2.add_row([int(ID[i]), ra[i], dec[i], '{:.2f}'.format(mag[i]), int(ncams[i]), int(numConPerTar[i])])
# print('\nUser have requested data table of all targets')
# print(t2)
#----------------------------

# DISTRIBUTION PLOTS

# Necessary to define again in case a PIC sample is defined
magRange = [np.min(df['mag']), np.max(df['mag'])]
fig3 = plt.subplots(2, 2, figsize=(14,8))
try: pt.plotStellarSampleDistributions(fig3, df['mag'].to_numpy(), dfc['mag'].to_numpy(),
                                       magRange, df['ncon'].to_numpy(), dfc['dis'].to_numpy())
except: fig3 = False
if plot: plt.show()

#==============================================================#
#                          SAVE OUTPUT                         #
#==============================================================#

errorcode('module', '\nPrologue')

if outputDir is not None:
    res = 200
    print('Saving all plots to {0}'.format(outputDir))
    fig0.savefig(outputDir / f'plot_{outputPrefix}_allsky.png',       bbox_inches='tight', dpi=res)
    fig1.savefig(outputDir / f'plot_{outputPrefix}_pointing.png',     bbox_inches='tight', dpi=res)
    fig2.savefig(outputDir / f'plot_{outputPrefix}_TeffvsRadius.png', bbox_inches='tight', dpi=res)
    if not fig3 is False: fig3[0].savefig(outputDir / f'plot_{outputPrefix}_distribution.png', bbox_inches='tight', dpi=res)

    # Save textfile with information about PIC stars
    log = ('PIC Catalogue include\n' +
           f'PLATO field           : {platoField}\n' +
           f'PLATO sample          : {starSample}\n' +
           f'PLATO camera-group    : {camGroup}\n' +
           f'Magnitude range       : {magRange[0]:.2f}-{magRange[1]:.2f} mag\n' +
           f'Luminosity class      : {lumClass}\n' +
           f'Spectral type         : {specType}\n' +
           f'Contaminants distance : {disConLimit} arcsec\n' +
           f'Contaminant magnitude : {dmagConLimit} mag\n' +
           tt + '\n' + t.get_string())
    with open(outputDir / f'{outputPrefix}.log','w') as file:
        file.write(log)

    # Only save the catalogs if they don't exist
    if fileFormat != '.txt' and subfield is False:

        # We reset the index in order to save to feather
        df  = df.reset_index()
        dfc = dfc.reset_index()
        print(f'Saving file {outputFileTar}')
        df.to_feather(outputFileTar)
        print(f'Saving file {outputFileCon}')
        dfc.to_feather(outputFileCon)

        # Output file name of HPC data
        outputVarFileName = f'{outputDir}/cluster_{outputPrefix}.txt'
        print(f'Saving file {outputVarFileName}')
        headerVar = 'PIC, M, R, Teff'
        outputVarData = np.transpose([df['ID'].to_numpy(), df['M'].to_numpy(),
                                      df['R'].to_numpy(),  df['Teff'].to_numpy()])
        np.savetxt(outputVarFileName, outputVarData, header=headerVar, comments='',
                   fmt=['%i', '%0.3f', '%0.3f', '%i'], delimiter=',')

    if fileFormat == '.csv' and subfield is True:
        print(f'Saving file {outputFileFOV}')
        outputConData = np.transpose([np.append(df['ra'].to_numpy(),  dfc['ra'].to_numpy()),
                                      np.append(df['dec'].to_numpy(), dfc['dec'].to_numpy()),
                                      np.append(df['mag'].to_numpy(), dfc['mag'].to_numpy())])
        np.savetxt(outputFileFOV, outputConData, fmt=['%.6f', '%.6f', '%.3f'])

# Finito!
toc = datetime.datetime.now()
print('Total execution time: {0} [hh:mm:ss]'.format(toc-tic))
print('')
