
"""
Usage: simQuarter.py <camaraGroupNr> <cameraNr> <quarterNr>

cameraGroupNr: either 1,2,3, or 4
cameraNr: either 1,2,3,4,5 or 6
quarterNr: either 1,2,3,4,5,6,7 or 8

Example: $ python3 simQuarter.py 2 5 6
"""


import os
import sys
import math

import numpy as np
from simulation import Simulation
from referenceFrames import getCCDandPixelCoordinates
from referenceFrames import platformToTelescopePointingCoordinates
from referenceFrames import sunSkyCoordinatesAwayfromPlatformPointing
from referenceFrames import CCD



inputDir = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"


#--- Check the number of input arguments

if len(sys.argv) != 4:
    print("Usage:   $ python3 simQuarter.py <cameraGroupNr> <cameraNr> <quarterNr>")
    print("Example: $ python3 simQuarter.py 2 5 6")
    exit(1)



#--- Configuration parameters

inputFile = inputDir + "/inputmagali.yaml"
outputDir = os.getenv("PLATO_PROJECT_HOME") 
outputPrefix = "Run1"
print("Using " + inputFile + " as inputfile")
print("Writing output to " + outputDir + outputPrefix + "_Q*_group*_camera*.hdf5")

raPlatform  = np.deg2rad(171.675)       # Platform right ascension pointing coordinate 
decPlatform = np.deg2rad(3.005)         # Platform declination pointing coordinate

raCenter  = np.deg2rad(171.675)         # Right ascension on which to centre the subfield
decCenter = np.deg2rad(3.005)           # Declination on which to centre the subfield

numColumnsSubField = 20                 # Number of columns in the modelled sub-field [pixels]
numRowsSubField = 20                    # Number of rows in the modelled sub-field [pixels]

#--- End configuration parameters


# Plato has 4 groups, of each 6 telescopes. Each quarter of the year, the platform
# (not the telescopes!) is rotated along its roll axis to repoint the solar panels
# towards the Sun.

# Select which camera from the arguments with which the script is called 

group = int(sys.argv[1])
telescope = int(sys.argv[2])
quarter = int(sys.argv[3])


print("Configuring PlatoSim for quarter {0} of camera {1} of group {2}".format(quarter, telescope, group))

# Output will be stored in e.g. Run1_Q1_group2_camera7.hdf5

outputFilePrefix = outputPrefix + "_Q{0:1d}_group{1:1d}_camera{2:1d}".format(quarter, group, telescope)
sim = Simulation(outputFilePrefix, inputFile)
sim.outputDir = outputDir

# Set the simulation parameters that are the same for any quarter and for any telescope

sim["ObservingParameters/RApointing"] = np.rad2deg(raPlatform)
sim["ObservingParameters/DecPointing"] = np.rad2deg(decPlatform)
sim["SubField/NumColumns"] = numColumnsSubField
sim["SubField/NumRows"] = numRowsSubField

# Set the telescope group ID, this is needed for the subfield calculations later on.
  
sim["Telescope/GroupID"] = group

# Set the focal plane angle different per group. 45.0 is to ensure that the platform pointing axis
# falls on a CCD rather than just in between CCDs. The (group - 1) * 90.0 is convenient so that 
# on the sky the left/right/top/down labeling of the different CCDs is the same for each group.

sim["Camera/FocalPlaneOrientation/Source"] = "ConstantValue"
sim["Camera/FocalPlaneOrientation/ConstantValue"] = 45.0 + (group - 1) * 90.0


# Set the quarter specific parameters

sim["RandomSeeds/JitterSeed"] = 2033429158 + 100 * quarter
sim["Platform/SolarPanelOrientation"] = math.fmod(quarter * 90., 360.)         # 0, 90, 180, and 270 degrees for Q1, Q2, Q3, and Q4

cycleTime = sim["ObservingParameters/CycleTime"]
readoutTime, dummy = sim.getReadoutTime()
numExposures = 10
#numExposures = int(365.25 / 4 * 86400 / cycleTime)
sim["ObservingParameters/NumExposures"] = numExposures
sim["ObservingParameters/BeginExposureNr"] = (quarter-1) * numExposures  

# Attempt to set a subfield around the specified coordinates on one of the 4 CCDs of the telescope.
# This will fail (return value == False) if the subfield is not visible by any of the 4 CCDs or
# that the subfield is too large to entirely fit on a CCD. 
# If successful, the function sets the CCD and subfield parameters correctly in the 'sim' object.

isSuccessful =  sim.setSubfieldAroundCoordinates(raCenter, decCenter, numColumnsSubField, numRowsSubField, normal=True)

if isSuccessful:
    # Make sure that the following random seeds differ for each telescope and for each quarter
    # We assume a maximum of 8 quarter and 4 camera groups

    randomSeedOffset = 1000 * 8 * quarter +  10 * 4 * group + telescope

    sim["RandomSeeds/ReadOutNoiseSeed"] = 1424949740 + randomSeedOffset
    sim["RandomSeeds/PhotonNoiseSeed"]  = 1533320336 + randomSeedOffset
    sim["RandomSeeds/FlatFieldSeed"]    = 1633320381 + randomSeedOffset 
    sim["RandomSeeds/DriftSeed"]        = 1733429158 + randomSeedOffset  

    # Run the PlatoSim simulator 
    # logLevel can 1 (least verbose) to 3 (most verbose)

    print("Launching PlatoSim3 for {0} exposures".format(numExposures))
    simFile = sim.run(logLevel=1)

else:
    print("Sub-field does not lay entirely on any of the CCDs of telescope {0} of group {1} in quarter Q{2}".format(telescope, group, quarter))

