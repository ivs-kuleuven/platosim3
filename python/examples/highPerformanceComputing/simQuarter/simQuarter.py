#!/usr/bin/env python3

"""
Written by Joris De Ridder and updated Nicholas Jannsen.
In this example we show how you can simulate an arbitary number of mission quarters. Since a nescessary rotation (along the roll axis) of the spacecraft platform is required in order to repoint the solar panels every 90 days, simulations cannot realistically extend beyond a quarter of a year. This together with the fact that the PLATO spacecraft will be equipped with 4 camera groups consisting of 6 cameras each, constrains the efficient use of nodes and cores. In order to make the simulation as realistic as possible, random seats of various intrinsic- and instrumental effects needs to be included, meaning that each camera within each camera group needs to be simulated independently (which indeed are the raw imaging output of PLATO). The MapReduce approach here boils down to defining a function that simulates a timeseries of 88 days (where 2 days are used for rotation and stabilazation), given an `inputfile.yaml`, the camera group (1-4), the camera within that group (1-6), and a quarter (1-N) out of N quarters simulated.
"""

import os
import sys
import math
import numpy as np
from colorama import Fore, Style, Back
from platosim.simulation import Simulation
from platosim.referenceFrames import getCCDandPixelCoordinates
from platosim.referenceFrames import platformToTelescopePointingCoordinates
from platosim.referenceFrames import sunSkyCoordinatesAwayfromPlatformPointing
from platosim.referenceFrames import CCD

#--------------------------------------------------------------#
#                     CONFIGURE PARAMETERS                     #
#--------------------------------------------------------------#

inputDir     = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"
outputDir    = os.getcwd() + "/"
simPrefix = "run1"

#--------------------------------------------------------------#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#--------------------------------------------------------------#

# Help usage function:
usage = \
"""
Usage: python {0} <inputfile> <cameraGroupNr> <cameraNr> <quarterNr> [<logLevel>]

inputfile    : PlatoSim yaml inputfile
cameraGroupNr: either 1, 2, 3, 4
cameraNr     : either 1, 2, 3, 4, 5, 6
quarterNr    : either 1, 2, 3, 4, 5, 6, 7, 8 (2 years of mission)
logLevel     : either 1, 2, 3. Least verbose: 1 (default), most verbose: 3.

Before running this script, check the user defined parameterss within and modify
them for your needs. Then run the script by

Example: $ python simQuarter.py inputfile.yaml 2 5 6
Example: $ python simQuarter.py inputfile.yaml 2 4 5 3
""".format(sys.argv[0])

def error(): print(Fore.RED + Style.BRIGHT + '[ERROR]: Wrong input!' + Style.RESET_ALL)
def help(): print(Fore.BLUE + Style.BRIGHT + usage + Style.RESET_ALL); exit()

# Print usage if no arguments are given:
if (len(sys.argv) == 1): help()

# Check for opvious mistakes
if (len(sys.argv) < 5) or (len(sys.argv) > 6): error(); help()

# Handle log level
if len(sys.argv) == 6:
    logLevel = int(sys.argv[5])
else:
    logLevel = 1

#--------------------------------------------------------------#
#                   START OF PLATOSIM SETTINGS                 #
#--------------------------------------------------------------#

# Plato has 4 groups, of each 6 telescopes. Each quarter of the year, the platform
# (not the telescopes!) is rotated along its roll axis to repoint the solar panels
# towards the Sun.

# Select which camera from the arguments with which the script is called

inputFile = sys.argv[1]
group     = int(sys.argv[2])
telescope = int(sys.argv[3])
quarter   = int(sys.argv[4])

print("Using " + inputFile + " as inputfile")
print("Configuring PlatoSim for telescope {0}.{1} Q{2}".format(group, telescope, quarter))
#print("Writing output to " + outputDir + simPrefix + "_group{0}_camera{1}_Q{2}.hdf5".format(group, telescope, quarter))


# Output will be stored in the output file

outputFilePrefix = simPrefix + "_{0:1d}_{1:1d}_{2:1d}".format(group, telescope, quarter)
sim = Simulation(outputFilePrefix, inputFile)
sim.outputDir = outputDir

# Set the simulation parameters that are the same for any quarter and for any telescope
# The subfield will be selected so that it's right on the platform pointing axis.
# This ensures that it's visible by all cameras.

raPlatform  = sim["ObservingParameters/RApointing"]    # [deg]
decPlatform = sim["ObservingParameters/DecPointing"]   # [deg]
raCenter    = np.deg2rad(raPlatform)                   # [rad]
decCenter   = np.deg2rad(decPlatform)                  # [rad]

numColumnsSubField = sim["SubField/NumColumns"]
numRowsSubField    = sim["SubField/NumRows"]

# Set the telescope group ID, this is needed for the subfield calculations later on.

sim["Telescope/GroupID"] = group

# Set the quarter specific parameters

sim["RandomSeeds/JitterSeed"] = 2033429158 + 100 * quarter
sim["Platform/SolarPanelOrientation"] = math.fmod(quarter * 90., 360.)  # 0, 90, 180, and 270 degrees for Q1, Q2, Q3, and Q4

cycleTime = sim["ObservingParameters/CycleTime"]
numExposuresCoveringOneQuarter = 90. * 86400. / cycleTime               # One quarter is 90 days
numExposures = (90. - 2.) * 86400. / cycleTime                          # Two days lost because of platform roll + thermal stabilisation
#numExposures = 100                                                     # For testing only
sim["ObservingParameters/NumExposures"] = int(numExposures)
sim["ObservingParameters/BeginExposureNr"] = (quarter-1) * int(numExposuresCoveringOneQuarter)

# Attempt to set a subfield around the specified coordinates on one of the 4 CCDs of the telescope.
# This will fail (return value == False) if the subfield is not visible by any of the 4 CCDs or
# that the subfield is too large to entirely fit on a CCD.
# If successful, the function sets the CCD and subfield parameters correctly in the 'sim' object.

isSuccessful = sim.setSubfieldAroundCoordinates(raCenter, decCenter, numColumnsSubField, numRowsSubField, normal=True)

if isSuccessful:
    # Make sure that the following random seeds differ for each telescope and for each quarter
    # We assume a maximum of 8 quarter and 4 camera groups

    randomSeedOffset = 1000 * 8 * quarter + 10 * 4 * group + telescope

    sim["RandomSeeds/ReadOutNoiseSeed"] = 1424949740 + randomSeedOffset
    sim["RandomSeeds/PhotonNoiseSeed"]  = 1533320336 + randomSeedOffset
    sim["RandomSeeds/FlatFieldSeed"]    = 1633320381 + randomSeedOffset
    sim["RandomSeeds/DriftSeed"]        = 1733429158 + randomSeedOffset
    sim["RandomSeeds/CosmicSeed"]       = 1494750830 + randomSeedOffset

	# Run the PlatoSim simulator
    # logLevel can 1 (least verbose) to 3 (most verbose)

    print("Launching PlatoSim for {0} exposures".format(numExposures))
    simFile = sim.run(logLevel=logLevel)

else:
    print("Sub-field does not lay entirely on any of the CCDs of telescope {0}.{1} in Q{2}".format(group, telescope, quarter))
