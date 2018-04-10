
# Example script how to use the SimFile and Simulation classes
# Run this script using:
#  $ ipython demoCreateStarCatalog.py

import os
import numpy as np

from simfile import SimFile
from simulation import Simulation



# Specify the absolute paths of some of the input files and the output folder.

inputDir    = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"
inputFile   = inputDir + "/inputfile.yaml"

outputDir   = os.getcwd()
outputFile  = "DemoOutput"

starCatalogFileName = inputDir + "/myCatalog.txt"

# Set up a Simulation object

sim = Simulation(outputFile, inputFile)
sim.outputDir = outputDir

sim["Platform/UseJitter"] = "no"
sim["PSF/Model"] = "MappedGaussian"


# Specify the orientation of the platform, telescope, focal plane, etc.
# Do this before creating the star catalog.

sim["ObservingParameters/RApointing"]             = 180.0 
sim["ObservingParameters/DecPointing"]            = -70.0
sim["Telescope/AzimuthAngle"]                     =   0.0
sim["Telescope/TiltAngle"]                        =   0.0
sim["Camera/FocalPlaneOrientation/Source"]        = "ConstantValue"
sim["Camera/FocalPlaneOrientation/ConstantValue"] =   5.0
sim["CCD/OriginOffsetX"]                          =   0.0
sim["CCD/OriginOffsetY"]                          =   0.0
sim["CCD/Orientation"]                            =  90.0
sim["Camera/IncludeFieldDistortion"]              =  "no"
sim["Camera/FocalLength/Source"]                  = "ConstantValue"
sim["Camera/FocalLength/ConstantValue"]           = 0.24712595 

# Specify the pixel coordinates (of the CCD, not of the subfield) of your stars

row = np.array([7.0, 10.0, 15.0])
col = np.array([40., 45.0, 50.0])
magnitude = np.array([11.0, 10.0, 12.0])
starID = [100, 101, 102]


# Create the star catalog file: an ascii file will be written with the columns
# ra, dec, and magnitude.

sim.createStarCatalogFileFromPixelCoordinates(row, col, magnitude, starID, starCatalogFileName)


# Make sure the simulation object uses this star catalog

sim["ObservingParameters/StarCatalogFile"] = starCatalogFileName


# Set a subfield around the stars, large enough to contain all of them

sim["SubField/ZeroPointRow"]    = 5
sim["SubField/ZeroPointColumn"] = 30
sim["SubField/NumColumns"]      = 25
sim["SubField/NumRows"]         = 25


# Run the simulation

simFile = sim.run()
