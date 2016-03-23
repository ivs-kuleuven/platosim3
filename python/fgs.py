import numpy as np

from simfile import SimFile
from simulation import Simulation
from referenceFrames import setSubfieldAroundCoordinates

# Specify the absolute paths of some of the input files and the output folder

myInputs    = "/Users/joris/Development/Cpp/PlatoSim3/inputfiles"
myInputs    = "/Users/rik/Work/PLATO/myInputs"

inputFile   = myInputs + "/myInputfile.yaml"

starCatalog = myInputs + "/guide_stars_EQ.txt"
jitterFile  = myInputs + "/PlatoJitter_Airbus.txt"
psfFile     = myInputs + "/psf.hdf5"

outputDir   = "/Users/joris/Development/Cpp/PlatoSim3/python"
outputDir   = "/Users/rik/Work/PLATO/Simulations"
outputFilePrefix = "GuideStarThalesFine"

# Read the guide star catalog

ra, dec, V = np.loadtxt(starCatalog, unpack=True)
NguideStars = len(ra)


RA_OPTICAL_AXIS = 86.7987057905419
DEC_OPTICAL_AXIS = -46.395950854582

# For each guide star, center a subfield around it, and run the simulator

for n in range(NguideStars):

    print("Running the simulator for guide star {0}".format(n))
    print("Guide Star Coordinates [deg]: {}, {}".format(ra[n], dec[n]))

    # Set up a Simulation object

    sim = Simulation(outputFilePrefix + "{0:02d}".format(n), inputFile)
    sim.outputDir = outputDir

    # Point the spacecraft. The coordinates refer to location of the optical axis.

    sim["ObservingParameters/RApointing"] = RA_OPTICAL_AXIS
    sim["ObservingParameters/DecPointing"] = DEC_OPTICAL_AXIS


    # Center the subfield around the current guide star
    # First extract the required information from the XML.
    # Note that for this simulation, we want to use the fast cams, not the nominal ones.

    raOpticalAxis   = np.deg2rad(RA_OPTICAL_AXIS)
    decOpticalAxis  = np.deg2rad(DEC_OPTICAL_AXIS)
    focalPlaneAngle = np.deg2rad(float(sim["Camera/FocalPlaneOrientation"]))
    focalLength     = float(sim["Camera/FocalLength"]) * 1000.0
    pixelSize       = int(sim["CCD/PixelSize"])    # [micron]
    plateScale      = float(sim["Camera/PlateScale"])   # [arcsec/micron]

    subfieldSizeX = 9     # [pixels]
    subfieldSizeY = 9     # [pixels]

    nominalCamera = False

    hasCcdCode = setSubfieldAroundCoordinates(sim, np.deg2rad(ra[n]), np.deg2rad(dec[n]), 
                                              subfieldSizeX, subfieldSizeY, focalLength, plateScale, pixelSize, \
                                              raOpticalAxis, decOpticalAxis, focalPlaneAngle, nominalCamera)

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

    sim["CCD/IncludePhotonNoise"]       = "yes"

    # Set some other input parameters specific to this simulation

    sim["ObservingParameters/NumExposures"] = 10
    sim["ObservingParameters/SkyBackground"] = 150
    sim["ObservingParameters/ExposureTime"] = 2.25
    sim["ObservingParameters/StarCatalogFile"] = starCatalog

    sim["Telescope/TransmissionEfficiency"] = 0.5

    sim["CCD/Gain"] = 58
    sim["CCD/QuantumEfficiency"] = 0.5
    sim["CCD/FullWellSaturation"] = 1243000
    sim["CCD/DigitalSaturation"] = 16383
    sim["CCD/ReadoutNoise"] = 130
    sim["CCD/ElectronicOffset"] = 100
    sim["CCD/FlatfieldPtPNoise"] = 0.016
    sim["CCD/ReadoutTime"] = 0.25

    sim["SubField/NumBiasPrescanRows"] = 0
    sim["SubField/NumSmearingOverscanRows"] = 0
    sim["SubField/SubPixels"] = 128

    sim["PSF/UseGauss"] = "yes"
    sim["PSF/Sigma"] = .025
    sim["PSF/NumberOfPixels"] = 8
    sim["PSF/NumberOfSubPixels"] = 128
    sim["PSF/Filename"] = psfFile

    sim["Platform/UseJitter"] = "yes"
    sim["Platform/UseJitterFromFile"] = "yes"
    sim["Platform/JitterFileName"] = jitterFile


    # Run the simulation without the flux extraction, with an HDF5 file as output

    simFile = sim.run()
 
 
