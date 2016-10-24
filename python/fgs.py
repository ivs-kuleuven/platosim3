# Example script how to use the SimFile and Simulation classes
# Run this script using:
#  $ ipython fgs.py

import os
import numpy as np

from simfile import SimFile
from simulation import Simulation
from referenceFrames import setSubfieldAroundCoordinates

# Specify the absolute paths of some of the input files and the output folder.
# The following default values will always work, but we advice you not to use
# them, but make your own input and output folders (and specify their paths here),
# so that you don't pollute the PlatoSim base inputfiles/ or python/ folders.

inputDir    = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"

inputFile   = inputDir + "/inputfgs.yaml"
starCatalog = inputDir + "/guide_stars_EQ.txt"
jitterFile  = inputDir + "/PlatoJitter_Airbus.txt"
psfFile     = inputDir + "/psf.hdf5"

outputDir   = os.getcwd()
outputFilePrefix = "/GuideStarThalesFine"

# Read the guide star catalog

ra, dec, V = np.loadtxt(starCatalog, unpack=True)
NguideStars = len(ra)

# For each guide star, center a subfield around it, and run the simulator

for n in range(NguideStars):

    print("Running the simulator for guide star {0}".format(n))
    print("Guide Star Coordinates [deg]: {}, {}".format(ra[n], dec[n]))

    # Set up a Simulation object

    sim = Simulation(outputFilePrefix + "{0:02d}".format(n), inputFile)
    sim.outputDir = outputDir

    # Make sure it uses the right starCatalog, jitter file, and PSF file

    sim["ObservingParameters/StarCatalogFile"] = starCatalog
    sim["Platform/JitterFileName"] = jitterFile
    sim["PSF/FromFile/Filename"] = psfFile 

    # Center the subfield around the current guide star
    # First extract the required information from the yaml input file.
    # Note that for this simulation, we want to use the fast cams, not the nominal ones.
    #
    # This function sets the following configuration parameters:
    # 
    # Camera/FieldDistortion/IncludeFieldDistortion = True
    # 
    # CCD/OriginOffsetX
    # CCD/OriginOffsetY
    # CCD/Orientation
    # CCD/NumColumns
    # CCD/NumRows
    # CCD/ReadoutTime
    # 
    # SubField/ZeroPointRow
    # SubField/ZeroPointColumn
    # SubField/NumRows
    # SubField/NumColumns
    # 
    # ObservingParameters/ExposureTime
    # 

    subfieldSizeX = 9     # column width [pixels]
    subfieldSizeY = 9     # row width [pixels]
    normalCamera = False

    hasCcdCode = setSubfieldAroundCoordinates(sim, np.deg2rad(ra[n]), np.deg2rad(dec[n]), subfieldSizeX, subfieldSizeY, normalCamera)

    # If the star does not fall on a CCD, or is too close to the edge, skip it.

    if not hasCcdCode:
        continue

 
    # Make sure that the random seeds are set different for each guide star,
    # so that they have all different noise realisations.

    sim["RandomSeeds/PhotonNoiseSeed"]  = 1433237514 + n
    sim["RandomSeeds/ReadOutNoiseSeed"] = 1424949740 + n
    sim["RandomSeeds/FlatFieldSeed"]    = 1425284070 + n
    sim["RandomSeeds/JitterSeed"]       = 1424967476 + n
    sim["RandomSeeds/CTESeed"]          = 1424949740 + n
    sim["RandomSeeds/DriftSeed"]        = 1433826961 + n


    # Run the simulation without the flux extraction, with an HDF5 file as output

    simFile = sim.run()
 
 
