
# Example script how to include variable stars in your simulation. 
# Note: the simulation concerns 100,000 exposures which takes a little while.
# Run this script using:
#  $ python demo_variableStars.py

import os
import numpy as np

from platosim.simfile    import SimFile
from platosim.simulation import Simulation



# Specify the absolute paths of some of the input files and the output folder.

inputDir    = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"
inputFile   = inputDir + "/inputfile.yaml"

outputDir   = os.getcwd()
outputFile  = "DemoOutput"

starCatalogFileName = inputDir + "/myCatalog.txt"

# Set up a Simulation object

sim = Simulation(outputFile, inputFile)
sim.outputDir = outputDir


# Set the simulation parameters

sim["ObservingParameters/NumExposures"]           = 100000
sim["ObservingParameters/RApointing"]             = 180.0 
sim["ObservingParameters/DecPointing"]            = -70.0
sim["ObservingParameters/StarCatalogFile"]        = starCatalogFileName

sim["Sky/IncludeVariableSources"]                 = "yes"
sim["Sky/VariableSourceList"]                     = "python/Examples/VariableStars/variableStarList.txt"
sim["Sky/IncludeCosmicsInSubField"]               = "no"
sim["Sky/IncludeCosmicsInSmearingMap"]            = "no"
sim["Sky/IncludeCosmicsInBiasMap"]                = "no"   

sim["Platform/UseJitter"]                         = "yes"
sim["Platform/JitterSource"]                      = "FromRedNoise"

sim["Telescope/GroupID"]                          = 1
sim["Telescope/UseDrift"]                         = "no"

sim["Camera/IncludeFieldDistortion"]              = "yes"

sim["PSF/Model"]                                  = "AnalyticNonGaussian"

sim["CCD/Position"]                               = 3
sim["CCD/IncludeConvolution"]                     = "no"

sim["SubField/ZeroPointRow"]                      = 2000
sim["SubField/ZeroPointColumn"]                   = 2000
sim["SubField/NumRows"]                           = 50
sim["SubField/NumColumns"]                        = 50


# Specify the pixel coordinates of the subfield (not CCD) of the stars
# Two stars, equal magnitudes, separations of 0.5, 1, 3, and 4 pixels.

col = [ 4.1,  4.6, 14.3, 15.3, 24.5, 27.5, 36.0, 42.0]
row = [ 4.1,  4.1,  4.3,  4.3,  4.7,  4.7,  4.5,  4.5]
mag = [12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0, 12.0]
starID = [100, 101, 102, 103, 104, 105, 106, 107]

# Two stars, mag diff of 4.0, separations of 0.5, 1, 3, and 4 pixels.

col = col + [ 6.1,  6.6, 16.3, 17.3, 26.5, 29.5, 38.0, 44.0]
row = row + [20.1, 20.1, 20.3, 20.3, 20.7, 20.7, 20.5, 20.5]
mag = mag + [ 9.0, 13.0,  9.0, 13.0,  9.0, 13.0,  9.0, 13.0]
starID = starID + [108, 109, 110, 111, 112, 113, 114, 115]

# Two stars, mag diff of 7.0, separations of 0.5, 1, 3, and 4 pixels.

col = col + [ 8.1,  8.6, 18.3, 19.3, 28.5, 31.5, 40.0, 46.0]
row = row + [41.1, 41.1, 41.3, 41.3, 41.7, 41.7, 41.5, 41.5]
mag = mag + [ 8.0, 15.0,  8.0, 15.0,  8.0, 15.0,  8.0, 15.0]
starID = starID + [116, 117, 118, 119, 120, 121, 122, 123]

# Convert from subfield to CCD pixel coordinates
 
row = np.array(row) + sim["SubField/ZeroPointRow"] 
col = np.array(col) + sim["SubField/ZeroPointColumn"] 
mag = np.array(mag)


# Create the star catalog file: an ascii file will be written with the columns
# ra, dec, and magnitude.

sim.createStarCatalogFileFromPixelCoordinates(row, col, mag, starID, starCatalogFileName)


# Run the simulation

simFile = sim.run()
