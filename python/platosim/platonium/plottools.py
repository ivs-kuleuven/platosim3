#!/usr/bin/env python3

"""
Given one or more HDF5 PlatoSim output files in original format, this script
can be used to produce all the necessary plots to understand the data output
from PlatoSim. Likewise, it's a tool to make the necessary plots for a PLATO
technical note for when we need to simulate data for the PLATO Consortium.
This code do not make.

User examples:
  $ python simTechNote output.hdf5
  $ python simTechNote *.hdf5 -j -s -v

The following diagnostic plots can be generated with this software:

* Subfield animation
    - Create a video of simulated imagettes/subfields
* Astrophysical variability
    - Raw and corrected time series showing mask-updates
    - PSD of raw time series
    - Star positions time series (if available)
* Spacecraft AOSC Jitter
    - Time series of yaw, pitch, roll
    - Power spectral density (PSD) of yaw, pitch, roll
* Camera Thermo-Elastic Distortion (TED)
    - Time series of yaw, pitch, and roll
    - Power Spectral Density (PSD) of yaw, pitch, and roll
* Time dependent instrumental characteristics
    - Transmission efficiency  : Time series
    - FEE temperature variation: Time series and PSD
    - CCD temperature variation: Time series and PSD

Author: Nicholas Jannsen
"""

import os
import h5py
import argparse
import datetime
import numpy as np

import matplotlib.pyplot as plt
import matplotlibrc

import platosim.referenceFrames as rf
from platosim.simulation import Simulation
from platosim.utilities import errorcode, normalize
from platosim.plot import (plotSubfieldAnimation, plotPhotometry,
                           plotMultiCameraAndQuarterPhotometry,
                           plotPSD)
\
# Start measure of execution time

tic = datetime.datetime.now()

#==============================================================#
#                   USER DEFINED PARAMETERS                    #
#==============================================================#

baseDir  = os.getenv("PLATO_PROJECT_HOME")
workDir  = os.getenv("PLATONIUM")
modelDir = workDir + '/models'
plot = True
save = True

# Constants

rad2arcsec = 648000 / np.pi
day2sec    = 60*60*24.

# Take time of coputations

tic = datetime.datetime.now()

#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

parser = argparse.ArgumentParser(epilog=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)

parser.add_argument('outputFiles', action='append', type=str, nargs='*',  help='HDF5 output file with simulated data')
parser.add_argument('-a', '--animation',  action='store_true', help='Produce animation of imagette time series')
parser.add_argument('-p', '--photometry', action='store_true', help='Produce plots of time series and PSD')
parser.add_argument('-s', '--starcoor',   action='store_true', help='Produce plots of injected variability signals')
parser.add_argument('-j', '--jitter',     action='store_true', help='Produce plots of spacecraft jitter')
parser.add_argument('-d', '--drift',      action='store_true', help='Produce plots of telescope drift')
parser.add_argument('-f', '--fee',        action='store_true', help='Produce plots of time dependent FEE temperature variations')
parser.add_argument('-c', '--ccd',        action='store_true', help='Produce plots of time dependent CCD temperature variations')
parser.add_argument('-t', '--te',         action='store_true', help='Produce plots of time dependent transmission efficiency')
parser.add_argument('-prime',             action='store_true', help='Uses the test data produced Prime/PPT')
parser.add_argument('-o', metavar='outputDir', type=str, help='Destination to save produced plots')

args = parser.parse_args()
filenames = args.outputFiles[0]
argA      = args.animation
argP      = args.photometry
argS      = args.starcoor
argJ      = args.jitter
argD      = args.drift
argF      = args.fee
argC      = args.ccd
argT      = args.te
fromPrime = args.prime
outputDir = args.o

# Select all plots if all are False
argsArray = np.array([argA, argP, argS, argJ, argD, argT, argF, argC])
if not any(argsArray):
    argA = True
    argP = True
    argS = True
    argJ = True
    argD = True
    argT = True
    argF = True
    argC = True

# Fetch information about observation
n = len(filenames)
group   = filenames[0][-11]
camera  = filenames[0][-9]
quarter = filenames[0][-6]

# Print code name to bash
errorcode('software', '\nTechnical Document Simulator\n')

#==============================================================#
#                      LOAD OUTPUT FILE(S)                     #
#==============================================================#

# First check file format

if not h5py.is_hdf5(filenames[0]):
    errorcode('error', 'File is not of the format HDF5!')
else:
    f = h5py.File(filenames[0], 'r')

    # TIME AND IMPORTANT INPUT

    #time = np.array(f['Telescope/Time'])

    raPointing      = f['InputParameters/ObservingParameters'].attrs['RApointing']
    decPointing     = f['InputParameters/ObservingParameters'].attrs['DecPointing']
    numImages       = f['InputParameters/ObservingParameters'].attrs['NumExposures']

    print('Observation by Camera {0}.{1} Q{2}'.format(group, camera, quarter))
    print('Telescope pointing (RA, Dec) : {}, {}'.format(raPointing, decPointing))
    print('Number of exposures          : {}'.format(numImages))

#--------------------------------------------------------------#
#                            ANIMATION                         #
#--------------------------------------------------------------#

if argA:

    errorcode('message', '\nAnimation')

    # ANIMATION OF IMAGETTES

    fig = plt.figure(figsize=(10,10))
    plotSubfieldAnimation(fig, filenames[0], showGrid=True, showStarPositions="PIC",
                          useTitle='Imagette of PIC 17157958 (10.83 mag)',
                          skipNimages=1000, colorMap='Spectral_r', outputFileName='animation')

#--------------------------------------------------------------#
#                           PHOTOMETRY                         # TODO
#--------------------------------------------------------------#

if argP is True:

    errorcode('message', '\nPhotometry')

    # TIMESERIES

    print('Plotting photometric timeseries of star')

    if len(filenames) == 1:
        fig = plt.figure(figsize=(10,10))
        plotPhotometry(fig, filenames[0], COB=False, title='Photometry')
    else:
        fig = plt.figure(figsize=(20,15))
        plotMultiCameraAndQuarterPhotometry(fig, filenames)
        exit()
    # POWER SPECTRAL DENSITY

    exit()
    print('Plotting PSD of photometric time series')
    titlePhotPSD   = 'PSD of photometric time series'
    legendsPhotPSD = ['Raw PSD']
    plotPSD(time, [flux*1e-6], titlePhotPSD, legendsPhotPSD, freqlim=1e-1)

    toc = datetime.datetime.now()
    print('Execution time for PlatoSim : {0} [hh:mm:ss]'.format(toc-tic))


    # TODO add to more light curves can be shown
    # if len(filenames) > 1:
    #     fig = plt.figure(figsize=(20,15))
    #     plotPhotometryComparison(fig, filenames)
    #     exit()


#--------------------------------------------------------------#
#                        SPACECRAFT JITTER                     #
#--------------------------------------------------------------#

if argJ is True:

    errorcode('message', '\nSpacecraft Jitter')

    # Fetch data using the PlatoSim script

    xJitter = np.array(f['ACS/Yaw'])
    yJitter = np.array(f['ACS/Pitch'])
    zJitter = np.array(f['ACS/Roll'])


    if len(data) == 3:
        print('Red Noise RMS yaw, pitch, and roll {0}, {1}, {2} [arcsec]'.format(data[0], data[1], data[2]))

    if len(data) == 1:

        # TIMESERIES

        print('Plotting jitter timeseries for yaw, pitch, and roll')
        rmsJitter    = data[0]
        timeJitter   = data[1][0]
        dataJitter   = data[1][1:]
        titleJitter  = 'Timeseries of spacecraft jitter'
        plotYawPitchRollJitter(time/(60*60), dataJitter, titleJitter, rms=rmsJitter, save=[outputDir, 'Jitter'])

        # POWER SPECTRAL DENSITY

        print('Plotting jitter PSD of raw, pitch, and roll jitter')
        titleJitterPSD   = 'PSD of spacecraft jitter'
        legendsJitterPSD = ['Yaw', 'Pitch', 'Roll', 'Jitter']
        plotYawPitchRollPSD(timeJitter, dataJitter, titleJitterPSD, legendsJitterPSD, freqlim=1e-4, save=[outputDir, 'Jitter'])

        # PARAMETER CORRELATION

        print('Plotting Jitter correlations between yaw, pitch, and roll')
        plotYawPitchRoll(timeJitter, dataJitter, rmsJitter, clabel='Time [seconds]', save=[outputDir, 'Jitter'])

#--------------------------------------------------------------#
#                    THERMAL ELASTIC DRIFT                     #
#--------------------------------------------------------------#

if argD is True:

    errorcode('message', '\nTelescope Thermo-Elastic Drift')

    # Fetch data using the PlatoSim script

    photometryClass = PhotometricFile(filenames[0])
    data = photometryClass.getThemoElasticDrift()

    if len(data) == 3:
        print('Red Noise RMS yaw, pitch, and roll {0}, {1}, {2} arcsec'.format(data[0], data[1], data[2]))

    if len(data) == 1:

        # TIMESERIES

        print('Plotting timeseries of drfit for raw, pitch, and roll ')
        rmsDrift    = data[0]
        timeDrift   = data[1][0]
        dataDrift   = data[1][1:]
        titleDrift  = 'Timeseries of camera thermo-elastic drift'
        labelsDrift = ['Time [days]', 'Yaw', 'Pitch', 'Roll', 'arcsec']
        plotTimeserie(time/day2sec, dataDrift, titleDrift, labelsDrift, rms=rmsDrift, save='Drift')

        # POWER SPECTRAL DENSITY

        print('Plotting PSD of drfit for raw, pitch, and roll ')
        titleDriftPSD   = 'PSD of camera thermo-elastic drift'
        legendsDriftPSD = ['Yaw', "Pitch", "Roll"]
        plotPSD(time, dataDrift, titleDriftPSD, legendsDriftPSD, lowlim=1e-6, save='Drift')

        # PARAMETER CORRELATIONS

        print('Plotting drift correlations between yaw, pitch, and roll')
        plotYawPitchRoll(time, dataDrift, rmsDrift, 'Time [seconds]', save='Drift')

#--------------------------------------------------------------#
#                    TRANSMISSION EFFICIENCY                   # TODO
#--------------------------------------------------------------#


if argT is True:

    # FEE TEMPERATURE VARIATION TODO

    feeTemperature = f['InputParameters/FEE'].attrs['Temperature']
    if feeTemperature == 'Nominal':
        feeTemperatureNominal = f['InputParameters/FEE'].attrs['NominalOperatingTemperature']
    else:
        feeTemperatureFileName = f['InputParameters/FEE'].attrs['TemperatureFileName'].decode('utf-8') # TODO file is missing!

# if useTE is True and argT is True:

#     errorcode('message', '\nTransmission efficiency')

#     print('Plotting timeseries of Transmission Efficiency')
#     plt.figure(figsize=(12,13))
#     plt.title('Timeseres of Transmission Efficiency', fontweight='bold')
#     plt.plot(time/day2sec, transmissionEfficiency, 'k-')
#     plt.axhline(teBOL, color='b', linestyle=':')
#     #plt.axhline(teEOL, color='r', linestyle=':')
#     plt.xlabel('Time [days]')
#     plt.ylabel('Efficiency')
#     plt.show()

#--------------------------------------------------------------#
#                    FEE TEMPERATURE VARIATION                 # TODO
#--------------------------------------------------------------#

#if useTemperatureVariations is True:


#--------------------------------------------------------------#
#                              PSF                             # TODO
#--------------------------------------------------------------#

    # psfsub = simFile.getPsf("rebinnedPSFsubPixel")
    # #psfsub = array(outputFile["/PSF/rebinnedPSFsubPixel"])
    # number_subpix = outputFile["/InputParameters/SubField"].attrs["SubPixels"]
    # sigmaPSF = phot.computePSFsigma(psfsub, number_subpix)
    # print("Every pixel has {0}x{0} subpixels".format(number_subpix))
    # print("sigma of the PSF = ", sigmaPSF)
    # plt.pcolormesh(psfsub, cmap=cm.nipy_spectral, vmax=0.01)
    # plt.grid(True, which='major', axis='both', linestyle='-', color='w')
    # plt.xticks(np.arange(0, psfsub.shape[0], number_subpix))
    # plt.yticks(np.arange(0, psfsub.shape[1], number_subpix))
    # plt.xlim(0, psfsub.shape[0])
    # plt.ylim(0, psfsub.shape[1])
    # print("showing the PSF at subpixel level...")
    # plt.show()


# Measure time of script

toc = datetime.datetime.now()
print('Execution time : {0} [hh:mm:ss]\n'.format(toc-tic))



# NOTE
# This block can be used to calculate the PSD for data with uneven time sampling
# Each of PSD computation (yaw, pitch, roll) takes 3.5 hours using 3 threads!
# The jitter files are to be found the extcloud/Platoman/models folder
#----------------------------------------------------------
#jitterFileName = workDir + '/models/jitterFiles_2020-01-10/01_PLATO_PDR_FPM_02_longrun_APE.csv'
# timeStart = datetime.datetime.now()
# fx, sx = pergrams.scargle(time, xJitter, threads='max', norm='density')
# fy, sy = pergrams.scargle(time, yJitter, threads='max', norm='density')
# fz, sz = pergrams.scargle(time, zJitter, threads='max', norm='density')
# timeEnd = datetime.datetime.now() - timeStart
# print('Execution time : {0} [hh:mm:ss]'.format(timeEnd))
# print('Loading PSD for yaw, pitch, and roll..')
# power_model_x = np.loadtxt(modelDir + '/jitterPSDyaw.txt')
# power_model_y = np.loadtxt(modelDir + '/jitterPSDpitch.txt')
# power_model_z = np.loadtxt(modelDir + '/jitterPSDroll.txt')
#----------------------------------------------------------

# NOTE
# with open(inputfile) as file:
#     content = yaml.load(file, Loader=yaml.FullLoader)
#     raPointing  = content.get('ObservingParameters').get('RApointing')
