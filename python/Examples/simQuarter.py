
"""
Usage: simQuarter.py <inputfile> <cameraGroupNr> <cameraNr> <quarterNr> [<logLevel>]

inputfile: Platosim yaml inputfile 
cameraGroupNr: either 1,2,3, or 4
cameraNr: either 1,2,3,4,5 or 6
quarterNr: either 1,2,3,4,5,6,7 or 8
logLevel: either 1,2, or 3. Least verbose: 1 (default), most verbose: 3.

Example: $ python3 simQuarter.py inputfile.yaml 2 5 6
"""


import os
import sys
import math

import numpy as np
from platosim.simulation import Simulation
from platosim.referenceFrames import getCCDandPixelCoordinates
from platosim.referenceFrames import platformToTelescopePointingCoordinates
from platosim.referenceFrames import sunSkyCoordinatesAwayfromPlatformPointing
from platosim.referenceFrames import CCD



inputDir = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"


#--- Check the number of input arguments

if (len(sys.argv) < 5) or (len(sys.argv) > 6):
    print("Usage:   $ python3 simQuarter.py <inputfile> <cameraGroupNr> <cameraNr> <quarterNr> [<logLevel>]")
    print("Example: $ python3 simQuarter.py inputfile.yaml 2 5 6")
    print("Example: $ python3 simQuarter.py inputfile.yaml 2 4 5 3")
    exit(1)


if len(sys.argv) == 6:
    logLevel = int(sys.argv[5])
else:
    logLevel = 1


#--- Configuration parameters

inputFile = sys.argv[1]
outputDir = os.getcwd() + "/"
outputPrefix = "Run1"

#--- End configuration parameters


# Plato has 4 groups, of each 6 telescopes. Each quarter of the year, the platform
# (not the telescopes!) is rotated along its roll axis to repoint the solar panels
# towards the Sun.

# Select which camera from the arguments with which the script is called 

group = int(sys.argv[2])
telescope = int(sys.argv[3])
quarter = int(sys.argv[4])

print("Using " + inputFile + " as inputfile")
print("Configuring PlatoSim for quarter {0} of camera {1} of group {2}".format(quarter, telescope, group))
print("Writing output to " + outputDir + outputPrefix + "_group{0}_camera{1}_Q{2}.hdf5".format(group, telescope, quarter))


# Output will be stored in e.g. Run1_Q1_group2_camera7.hdf5

outputFilePrefix = outputPrefix + "_group{0:1d}_camera{1:1d}_Q{2:1d}".format(group, telescope, quarter)
sim = Simulation(outputFilePrefix, inputFile)
sim.outputDir = outputDir

# Set the simulation parameters that are the same for any quarter and for any telescope
# The subfield will be selected so that it's right on the platform pointing axis.
# This ensures that it's visible by all cameras.

raPlatform = sim["ObservingParameters/RApointing"]                      # [deg]
decPlatform = sim["ObservingParameters/DecPointing"]                    # [deg]
raCenter = np.deg2rad(raPlatform)                                       # [rad]
decCenter = np.deg2rad(decPlatform)                                     # [rad]

numColumnsSubField = sim["SubField/NumColumns"]
numRowsSubField = sim["SubField/NumRows"]

# Set the telescope group ID, this is needed for the subfield calculations later on.
  
sim["Telescope/GroupID"] = group

# Set the quarter specific parameters

sim["RandomSeeds/JitterSeed"] = 2033429158 + 100 * quarter
sim["Platform/SolarPanelOrientation"] = math.fmod(quarter * 90., 360.)         # 0, 90, 180, and 270 degrees for Q1, Q2, Q3, and Q4

cycleTime = sim["ObservingParameters/CycleTime"]
numExposuresCoveringOneQuarter = 90. * 86400. / cycleTime                      # One quarter is 90 days
numExposures = (90. - 2.) * 86400. / cycleTime                                 # Two days lost because of platform roll + thermal stabilisation
#numExposures = 100                                                             # For testing only
sim["ObservingParameters/NumExposures"] = int(numExposures)
sim["ObservingParameters/BeginExposureNr"] = (quarter-1) * int(numExposuresCoveringOneQuarter)

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
    simFile = sim.run(logLevel=logLevel)

else:
    print("Sub-field does not lay entirely on any of the CCDs of telescope {0} of group {1} in quarter Q{2}".format(telescope, group, quarter))

