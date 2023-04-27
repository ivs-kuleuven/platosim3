#!/usr/bin/env python3
"""
This is an example script on how to compute long time series by splitting them up
and run them in parallel. The script specifies an array of jobs, for which we can
derive the total number of jobs, and the sequential number of the current job.
For full description run usage function: python runSimParallel.py
"""

import os
import sys
import time
from math import ceil
from colorama import Fore, Style, Back
from platosim.simfile import SimFile
from platosim.simulation import Simulation
timeStart = time.time()

#--------------------------------------------------------------#
#                PARSING COMMAND-LINE ARGUMENTS                #
#--------------------------------------------------------------#

usage = \
"""
Usage : {0} <arrayTaskID>

arrayTaskID : the task ID in the array of current task

This is an example script on how to compute long time series by splitting them up
and run them in parallel. The script specifies an array of jobs, for which we can
derive the total number of jobs, and the sequential number of the current job.

Don't run this script directly with Python, but launch the script runSimParallel.sh
which will in turn launch this python script. For example if you run on the prompt:

    $ module load worker
    $ wsub -t 1-10 -batch jobSequential.pbs

the time series will be into 10 different parts that are computed simultaneously,
and written to 10 different HDF5 output files.
""".format(sys.argv[0])

# Define error message and usage function
def error(): print(Fore.RED + Style.BRIGHT + '[ERROR]: Wrong input!' + Style.RESET_ALL)
def help(): print(Fore.BLUE + Style.BRIGHT + usage + Style.RESET_ALL); exit()

# Print usage if no arguments are given
if len(sys.argv) == 1: help()

# Check for opvious mistakes
if len(sys.argv) > 2:
    error(); help(); exit()

# The only input argument parsed is the ongoing job ID
jobNr = int(sys.argv[1])

#--------------------------------------------------------------#
#                       SETUP SIMULATION                       #
#--------------------------------------------------------------#

# Specify the absolute paths and filenames of input and output
inputDir   = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"
inputFile  = inputDir + "/inputfile.yaml"
outputDir  = os.getcwd()
outputFile = "outputSequential"

# Specify the total number of exposures of the entire unpartitioned time series
Nexps = 105

# Spefify the number of smaller simulations (jobs) that will run in parallel
# If selectSequence is True:  Njobs is the image range of a single job
# If selectSequence is False: Njobs is the total number of parallel simulations
selectSequence = True
N = 10

# Load simulation object (DO NOT CHANGE!)
sim = Simulation(outputFile+"_{0:04d}".format(jobNr), inputFile, outputDir=outputDir)

# Set the simulation parameters after need below
sim["PSF/Model"]                       = "MappedGaussian"
sim["ObservingParameters/RApointing"]  = 180.0
sim["ObservingParameters/DecPointing"] = -70.0
sim["Platform/SolarPanelOrientation"]  = 0.0
sim["Telescope/AzimuthAngle"]          = 0.0
sim["Telescope/TiltAngle"]             = 0.0

# Set the random seeds so that every segment of the partitioned time series uses different seeds
sim["RandomSeeds/ReadOutNoiseSeed"]  = 1424949740 + jobNr * 11111
sim["RandomSeeds/PhotonNoiseSeed"]   = 1433320336 + jobNr * 11111
sim["RandomSeeds/JitterSeed"]        = 1433320381 + jobNr * 11111
sim["RandomSeeds/FlatFieldSeed"]     = 1425284070 + jobNr * 11111
sim["RandomSeeds/DriftSeed"]         = 1433429158 + jobNr * 11111
sim["RandomSeeds/CosmicSeed"]        = 1494750830 + jobNr * 11111
sim["RandomSeeds/DarkSignalSeed"]    = 1468838669 + jobNr * 11111

#--------------------------------------------------------------#
#                       RUN SIMULATION                         #
#--------------------------------------------------------------#

# If Njobs should be interperted as an array sequence
if selectSequence is True:
    Njobs = int(ceil(float(Nexps) / N))
else:
    Njobs = N

# Calculate the nr of the begin exposure for this job
numberExposuresPerJob    = int(ceil(float(Nexps) / Njobs))
beginExposureNrOfThisJob = (jobNr-1) * numberExposuresPerJob + 1

# Set sequential input parameter for current job
sim["ObservingParameters/BeginExposureNr"] = beginExposureNrOfThisJob

# For the very last job the number of exposures needs to be adjusted because what's left
# after (Njobs-1) jobs may not be the full 'numExposuresPerJob'
if jobNr == Njobs:
    numberExposuresPerJob = Nexps - (Njobs-1) * numberExposuresPerJob
    sim["ObservingParameters/NumExposures"] = numberExposuresPerJob
else:
    sim["ObservingParameters/NumExposures"] = numberExposuresPerJob

# Print ongoing job for trouble shooting
inputs = {"a": Njobs, "b": Nexps, "c": jobNr, "d": beginExposureNrOfThisJob,
          "e": beginExposureNrOfThisJob+numberExposuresPerJob-1}
print("Monitoring (Njobs, Nexp)=({a}, {b}): jobNr={c}, jobExp={d}-{e}".format(**inputs))

# Run the simulation
#simFile = sim.run()

# Print ellapsed time to output
timeEnd = time.time()
print("Simulation No. {0:04d}".format(jobNr) + " took {} seconds".format(ceil(timeEnd-timeStart)))
