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

inputFile   = inputDir + "/inputfile.yaml"
starCatalog = inputDir + "/guide_stars_EQ.txt"
jitterFile  = inputDir + "/PlatoJitter_Airbus.txt"
psfFile     = inputDir + "/psf.hdf5"

outputDir   = os.getcwd()
outputFilePrefix = "/GuideStarThalesFine"

# Read the guide star catalog

ra, dec, V = np.loadtxt(starCatalog, unpack=True)
NguideStars = len(ra)

# Set the platform pointing coordinates (not the optical axis)

RA_PLATFORM = 86.7987057905419
DEC_PLATFORM = -46.395950854582

# For each guide star, center a subfield around it, and run the simulator

for n in range(NguideStars):

    print("Running the simulator for guide star {0}".format(n))
    print("Guide Star Coordinates [deg]: {}, {}".format(ra[n], dec[n]))

    # Set up a Simulation object

    sim = Simulation(outputFilePrefix + "{0:02d}".format(n), inputFile)
    sim.outputDir = outputDir

    # Point the spacecraft. The coordinates refer to location of the spacecraft roll axis,
    # not the optical axis of the telescope

    sim["ObservingParameters/RApointing"]  = RA_PLATFORM
    sim["ObservingParameters/DecPointing"] = DEC_PLATFORM


    # Specify the orientation of the telescope on the platform

    azimuthTelescope = 1.0     # [deg]
    tiltTelescope    = 1.0     # [deg]


    # Center the subfield around the current guide star
    # First extract the required information from the yaml input file.
    # Note that for this simulation, we want to use the fast cams, not the nominal ones.

    raPlatform      = np.deg2rad(RA_PLATFORM)
    decPlatform     = np.deg2rad(DEC_PLATFORM)
    focalPlaneAngle = np.deg2rad(float(sim["Camera/FocalPlaneOrientation"]))
    focalLength     = float(sim["Camera/FocalLength"]) * 1000   # [mm]
    pixelSize       = int(sim["CCD/PixelSize"])    # [micron]
    plateScale      = float(sim["Camera/PlateScale"])   # [arcsec/micron]

    subfieldSizeX = 9     # column width [pixels]
    subfieldSizeY = 9     # row width [pixels]

    nominalCamera = False
    includeFieldDistortion = False

    # This function sets the following configuration parameters:
    # 
    # Camera/FieldDistortion/IncludeFieldDistortion = True
    # 
    # CCD/OriginOffsetX
    # CCD/OriginOffsetY
    # CCD/Orientation
    # CCD/NumColumns
    # CCD/NumRows
    # 
    # SubField/ZeroPointRow
    # SubField/ZeroPointColumn
    # SubField/NumRows
    # SubField/NumColumns
    # 
    # ObservingParameters/ExposureTime
    # 

    hasCcdCode = setSubfieldAroundCoordinates(sim, np.deg2rad(ra[n]), np.deg2rad(dec[n]), subfieldSizeX, subfieldSizeY, 
                                              focalLength, plateScale, pixelSize, raPlatform, decPlatform, focalPlaneAngle, 
                                              np.deg2rad(azimuthTelescope), np.deg2rad(tiltTelescope), includeFieldDistortion, 
                                              nominalCamera)

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

    # Set some other input parameters specific to this simulation

    sim["ObservingParameters/NumExposures"] = 10
    sim["ObservingParameters/SkyBackground"] = 150
    sim["ObservingParameters/ExposureTime"] = 2.25
    sim["ObservingParameters/StarCatalogFile"] = starCatalog

    sim["Telescope/AzimuthAngle"] = azimuthTelescope           # [deg]
    sim["Telescope/TiltAngle"] = tiltTelescope                 # [deg]
    sim["Telescope/TransmissionEfficiency"] = 0.5

    sim["CCD/Gain"] = 58
    sim["CCD/QuantumEfficiency"] = 0.5
    sim["CCD/FullWellSaturation"] = 1243000
    sim["CCD/DigitalSaturation"] = 16383
    sim["CCD/ReadoutNoise"] = 130
    sim["CCD/ElectronicOffset"] = 100
    sim["CCD/FlatfieldPtPNoise"] = 0.016
    sim["CCD/ReadoutTime"] = 0.25
    sim["CCD/IncludePhotonNoise"]       = "yes"

    sim["SubField/NumBiasPrescanRows"] = 0
    sim["SubField/NumSmearingOverscanRows"] = 0
    sim["SubField/SubPixels"] = 8

    sim["PSF/Model"] = "Gaussian"
    sim["PSF/Gaussian/Sigma"] = .25
    sim["PSF/Gaussian/NumberOfPixels"] = 8

    sim["Platform/UseJitter"] = "yes"
    sim["Platform/UseJitterFromFile"] = "yes"
    sim["Platform/JitterFileName"] = jitterFile


    # Run the simulation without the flux extraction, with an HDF5 file as output

    simFile = sim.run()
 
 
