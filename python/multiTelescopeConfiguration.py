"""
Script to run the Plato Simulator for a multi-telescope configuration.  The 32 telescopes are arranged in four groups of 8 telescopes.  All telescopes in the 
same group have the same FOV and the lines of sight of the four groups are offset by an angle of 9.2 degrees from the PLM z-axis.

For each of the telescopes a simulation will be performed (i.e. a sub-field will be modelled with the Plato Simulator).  The sub-field has the same dimensions 
(in pixels for each telescope) and is always centred on the same coordinates (raCenter, decCenter).

Differences between the telescopes are:

- azimuth angle (the same within a group)
- seed for the readout noise
- seed for the photon noise
- seed for the flatfield map
- seed for the CTE map

All other configuration parameters are the same for the individual telescopes, incl.

- all platform-related parameters
- the PSF
- the field distortion polynomial

For each telescope, the script will automatically determine on which CCD the sub-field will be positioned.  A simulation will only be performed if the sub-field
falls on the CCD entirely.
"""

#########
# Imports
#########

import os
import math

from simulation import Simulation
from referenceFrames import getCCDandPixelCoordinates
from referenceFrames import platformToTelescopePointingCoordinates
from referenceFrames import drawCCDsInSky
from referenceFrames import CCD




################
# Input & Output
################

inputDir = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"
inputFile   = inputDir + "/inputfile.yaml"

outputDir   = os.getcwd()
outputPrefix = "MultiTelescopeConfiguration"





##########################
# Configuration parameters
##########################


# Observing parameters
######################

numExposures = 10                                           # Number of exposures (for each telescope)
exposureTime = 23                                           # Exposure time [s]
raPointing = 180.0                                          # Platform right ascension pointing coordinate [degrees]
decPointing = -70.0                                          # Platform declination pointing coordinate [degrees]
flux_m0 = 1.00238e8                                         # Photon flux of a V = 0 G2V star [photons / s / m^2 / nm]
skyBackground = 220.0                                       # Stellar + zodiacal background level [photons / pixel / s]
starCatalog = inputDir + "/starField_RA180Dec-70.txt"       # Filename of the star catalogue

#raCenter = 180.0                                            # Right ascension of the point about which to centre the sub-field [degrees]
#decCenter = 70.0                                            # Declination of the point about which to centre the sub-field [degrees] 

#A
#Position of the centre of the CCD:
#195.446447763, -13.4117113619
#B
#Position of the centre of the CCD:
#282.679384245, -34.3213999026
#C
#Position of the centre of the CCD:
#112.473919726, -13.4117274315
#D
#Position of the centre of the CCD:
#25.2409685516, -34.3214188293


#raCenter = 153.960227088 
#decCenter = -75.0766240288

raCenter = 227.747586623
decCenter = -63.9874188222


# Platform parameters (the same for all telescopes)
#####################

useJitter = "yes"                           # Do you want to account for the platform jitter?
useJitterFromFile = "no"                    # Do you want to read the jitter from a file?

jitterSeed = 1433320381                     # Random seed for jitter
jitterYawRms = 2.3                           # Jitter yaw RMS [arcsec]
jitterPitchRms = 2.3                        # Jitter pitch RMS [arcsec]
jitterRollRms = 2.3                         # Jitter roll RMS [arcsec]
jitterTimescale = 3600.0                    # Jitter timescale [s]

jitterFilename = inputDir + "/jitter.txt"   # Name of the jitter file (full path)


# Telescope parameters
######################

azimuthAngles = [45, 135, -135, -45]    # Azimuth angles of the telescopes (same within a group) [degrees]
tiltAngle = 9.2                         # Tilt angle of the telescopes (same for the 4 groups) [degrees]

lightCollectingArea = 113.1             # Effective area of a single telescope [cm^2]
transmissionEfficiency =  0.757         # Transmission efficiency                                             <--- use the same value for all the telescopes for now
driftYawRms = 2.3                       # Telescope drift yaw RMS [arcsec]                                    <--- use the same value for all the telescopes for now
driftPitchRms = 2.3                     # Telescope drift pitch RMS [arcsec]                                  <--- use the same value for all the telescopes for now
driftRollRms = 2.3                      # Telescope drift roll RMS [arcsec]                                   <--- use the same value for all the telescopes for now
driftTimescale = 3600.0                 # Telescope drift timescale [s]                                       <--- use the same value for all the telescopes for now

driftSeed = 1433429158                  # Random seed for telescope drift                                     <--- different for each telescope!


# Camera parameters
###################

focalPlaneAngle = 0                                                                                                 # Focal plane angle [degrees]
plateScale = 0.8333                                                                                                 # Plate scale [arcsec / micron]
focalLength = 0.24712595                                                                                            # Focal length as recovered from ZEMAX model [m]
throughputBandwidth = 400                                                                                           # FWHM of the throughput passband[nm]
throughputCentralWavelength = 600                                                                                   # Central wavelength of the throughput passband [nm]
includeFieldDistortion = "no"                                                                                      # Do you want to include field distortion?

fieldDistortionType = "Polynomial1D"                                                                                # Field distortion implementation (1D/2D polynomial)
fieldDistortionDegree = 3                                                                                           # Degree of the field distortion polynomial
fieldDistortionCoefficients = "[-0.0036696919678, 1.0008542317, -4.12553764817e-05, 5.7201219949e-06]"              # Coefficients of the field distortion polynomial
fieldDistortionInverseCoefficients = "[-0.00458067036444, 1.00110311283, -5.61136295937e-05, -4.311925329e-06]"     # Coefficients of the inverse field distortion polynomial


# PSF parameters    <--- in principle different for each telescope, but we have only one PSF
################

psfModel = "Gaussian"                           # Do you want to use a Gaussian PSF or read it from a file?

gaussianPsfSigma = 0.50                         # Standard deviation of the Gaussian PSF [pixels]
gaussianPsfNumPixels = 8                        # Number of pixels for which the Gaussian PSF is generated, in both directions

psfFromFileFilename = inputDir + "/psf.hdf5"    # Name of the file holding the PSF
psfFromFileDistanceToOA = -1                    # Auto-complete distance to the optical axis
psfFromFileRotationAngle = -1                   # Auto-complete rotation angle of the PSF w.r.t. the x-axis of the focal plane
psfFromFileNumPixels = 8                        # Number of pixel for which the PSF was generated, in both directions


# CCD parameters
################

# Automatically determined:
#    - origin offset (x, y)
#    - orientation
#    - size

pixelSize = 18                          # Pixel size [micron]
gain = 16                               # Detector gain [e- / ADU ]                                                 <--- use the same value for all the telescopes for now
quantumEfficiency = 0.8745              # Quantum efficiency                                                        <--- use the same value for all the telescopes for now
fullWellSaturation = 1000000            # Full-well saturation limit [e- / pixel]
digitalSaturation = 65535               # Digital saturation limit [ADU / pixel]
readoutNoise = 28                       # Readout noise [e- / pixel]                                                <--- use the same value for all the telescopes for now
electronicOffset = 100                  # Electronic offset [ADU]                                                   <--- use the same value for all the telescopes for now
readoutTime = 2                         # Readout time [s]
flatfieldP2PNoise = 0.016               # Flatfield peak-to-peak pixel noise                                        <--- use the same value for all the telescopes for now
cte = 0.99999                           # Mean CTE                                                                  <--- use the same value for all the telescopes for now

includeFlatfield = "yes"                # Do you want to account for the flatfield?
includePhotonNoise = "yes"              # Do you want to account for the photon noise?
includeReadoutNoise = "yes"             # Do you want to account for the readout noise?
includeCtiEffects = "yes"               # Do you want to account for the CTI effects?
includeOpenShutterSmearing = "yes"      # Do you want to account for the open-shutter smearing?
includeVignetting = "yes"               # Do you want to account for vignetting?
includeConvolution = "yes"              # Do you want to convolve with the PSF?
includeFullWellSaturation = "yes"       # Do you want to account for the full-well saturation (i.e. blooming)?
includeDigitalSaturation = "yes"        # Do you want to account for the digital saturation?
writeSubPixelImagesToHDF5 = "no"        # Do you want to store the sub-pixels images?

readoutNoiseSeed =  1424949740          # Random seed for the readout noise                                         <--- different for each telescope!
photonNoiseSeed = 1433320336            # Random seed for the photon noise                                          <--- different for each telescope!
flatfieldSeed = 1433320381              # Random seed for the flatfield                                             <--- different for each telescope!
cteSeed = 1424949740                    # Random seed for the CTE                                                   <--- different for each telescope!


# Sub-field parameters
######################

# Zeropoint depends on the pointing of the telescope

numColumnsSubField = 10         # Number of columns in the modelled sub-field [pixels]
numRowsSubField = 10            # Number of rows in the modelled sub-field [pixels]
numBiasPreScanRows = 5          # Number of rows in the pre-scan strip to determine the bias
numSmearingOverScanRows = 5     # Number of rows in the over-scan strip to determine the smearing
numSubPixelsPerPixel = 8        # Number of sub-pixels per pixels, in both directions, in the modelled sub-field





############
# Simulation
############

numTelescopeGroups = 1
numTelescopesPerGroup = 1

# Loop over all groups of telescopes

for group in range(numTelescopeGroups):
    
    # Loop over all telescopes in the current group
    
    for telescope in range(numTelescopesPerGroup):
        
        telescopeIndex = numTelescopesPerGroup * group + telescope
        
        # Output will be stored in /MultiTelescopeConfiguration/group<group>/telescope<telescope>
        
        outputFilePrefix = outputPrefix + "_group" +  "{0:02d}".format(group + 1) + "_telescope" + "{0:02d}".format(telescope + 1)
        sim = Simulation(outputFilePrefix, inputFile)
        sim.outputDir = outputDir
        
        # Compute the telescope pointing, based on the platform pointing, and the tilt and azimuth angle of the telescope
        
        raTelescope, decTelescope = platformToTelescopePointingCoordinates(math.radians(raPointing), math.radians(decPointing), math.radians(azimuthAngles[group]), math.radians(tiltAngle))

        
        print "Platform pointing: " + str(raPointing) + ", " + str(decPointing)
        print "Telescope pointing: " + str(math.degrees(raTelescope)) + ", " + str(math.degrees(decTelescope))
        
        # Determine on which CCD (A, B, C, or D) the coordinates (raCenter, decCenter) are positioned and at which location
        # (in pixel coordinates)
        
        includeFieldDistortionAsBoolean = (includeFieldDistortion == "yes")
        ccdCode, columnCenter, rowCenter = getCCDandPixelCoordinates(math.radians(raCenter), math.radians(decCenter), raTelescope, decTelescope, math.radians(focalPlaneAngle), focalLength, plateScale, pixelSize, includeFieldDistortionAsBoolean, nominal=True)
        
        # Check whether the sub-field falls entirely on the detector
        print ccdCode, columnCenter, rowCenter
        #if (ccdCode != None) and (rowCenter - numRowsSubField / 2 >= 0) and (rowCenter + numRowsSubField / 2 < CCD[ccdCode]["NRows"]) and (columnCenter - numColumnsSubField / 2 >= 0) and (columnCenter + numColumnsSubField / 2 < CCD[ccdCode]["NCols"]):
        if ccdCode != None:
            print "Processing simulation for telescope " + str(telescope + 1) + " of group " + str(group + 1)
        
            # Observing parameters
        
            sim["ObservingParameters/NumExposures"] = numExposures
            sim["ObservingParameters/ExposureTime"]  = exposureTime
            sim["ObservingParameters/RApointing"] = raPointing
            sim["ObservingParameters/DecPointing"] = decPointing
            sim["ObservingParameters/Fluxm0"] = flux_m0
            sim["ObservingParameters/SkyBackground"] = skyBackground
            sim["ObservingParameters/StarCatalogFile"] = starCatalog 
        
            # Platform parameters
        
            sim["Platform/UseJitter"] = useJitter
            sim["Platform/UseJitterFromFile"] = useJitterFromFile 
            sim["Platform/JitterYawRms"] = jitterYawRms 
            sim["Platform/JitterPitchRms"] = jitterPitchRms 
            sim["Platform/JitterRollRms"] = jitterRollRms 
            sim["Platform/JitterTimeScale"] = jitterTimescale
            sim["Platform/JitterFileName"] = jitterFilename  
        
            # Telescope parameters
        
            sim["Telescope/AzimuthAngle"] = azimuthAngles[group]
            sim["Telescope/TiltAngle"] = tiltAngle
            sim["Telescope/LightCollectingArea"] = lightCollectingArea
            sim["Telescope/TransmissionEfficiency"] = transmissionEfficiency 
            sim["Telescope/DriftYawRms"] = driftYawRms 
            sim["Telescope/DriftPitchRms"] = driftPitchRms 
            sim["Telescope/DriftRollRms"] = driftRollRms 
            sim["Telescope/DriftTimeScale"] = driftTimescale
        
            # Camera parameters
        
            sim["Camera/FocalPlaneOrientation"] = focalPlaneAngle
            sim["Camera/PlateScale"] = plateScale
            sim["Camera/FocalLength"] = focalLength
            sim["Camera/ThroughputBandwidth"] = throughputBandwidth
            sim["Camera/ThroughputLambdaC"] = throughputCentralWavelength
            sim["Camera/IncludeFieldDistortion"] = includeFieldDistortion
        
            # PSF parameters
        
            sim["PSF/Model"] = psfModel
            sim["PSF/Gaussian/Sigma"] = gaussianPsfSigma
            sim["PSF/Gaussian/NumberOfPixels"] = gaussianPsfNumPixels
            sim["PSF/FromFile/Filename"] = psfFromFileFilename
            sim["PSF/FromFile/DistanceToOA"] = psfFromFileDistanceToOA
            sim["PSF/FromFile/RotationAngle"] = psfFromFileRotationAngle
            sim["PSF/FromFile/NumberOfPixels"] = psfFromFileNumPixels
        
            # CCD parameters
        
            sim["CCD/OriginOffsetX"] = CCD[ccdCode]["zeroPointXmm"]
            sim["CCD/OriginOffsetY"] = CCD[ccdCode]["zeroPointYmm"]
            sim["CCD/Orientation"] = CCD[ccdCode]["angle"]
            sim["CCD/NumColumns"] = CCD[ccdCode]["Ncols"]
            sim["CCD/NumRows"] = CCD[ccdCode]["Nrows"]
        
            sim["CCD/PixelSize"] = pixelSize
            sim["CCD/Gain"] = gain
            sim["CCD/QuantumEfficiency"] = quantumEfficiency
            sim["CCD/FullWellSaturation"] = fullWellSaturation
            sim["CCD/DigitalSaturation"] = digitalSaturation
            sim["CCD/ReadoutNoise"] = readoutNoise
            sim["CCD/ElectronicOffset"] = electronicOffset
            sim["CCD/ReadoutTime"] = readoutTime
            sim["CCD/FlatfieldPtPNoise"] = flatfieldP2PNoise
            sim["CCD/CTEMean"] = cte
        
            sim["CCD/IncludeFlatfield"] =  includeFlatfield
            sim["CCD/IncludePhotonNoise"] =  includePhotonNoise
            sim["CCD/IncludeReadoutNoise"] = includeReadoutNoise
            sim["CCD/IncludeCTIeffects"] = includeCtiEffects
            sim["CCD/IncludeOpenShutterSmearing"] = includeOpenShutterSmearing
            sim["CCD/IncludeVignetting"] = includeVignetting
            sim["CCD/IncludeConvolution"] = includeConvolution
            sim["CCD/IncludeFullWellSaturation"] = includeFullWellSaturation
            sim["CCD/IncludeDigitalSaturation"] = includeDigitalSaturation
            sim["CCD/WriteSubPixelImagesToHDF5"] = writeSubPixelImagesToHDF5
        
            # Sub-field parameters
        
            sim["SubField/ZeroPointRow"] = rowCenter - numRowsSubField / 2
            sim["SubField/ZeroPointColumn"] = columnCenter - numColumnsSubField / 2
        
            sim["SubField/NumColumns"] = numColumnsSubField
            sim["SubField/NumRows"] = numRowsSubField
            sim["SubField/NumBiasPrescanRows"] = numBiasPreScanRows
            sim["SubField/NumSmearingOverscanRows"] = numSmearingOverScanRows
            sim["SubField/SubPixels"] = numSubPixelsPerPixel
        
            # Seed parameters
        
            sim["RandomSeeds/ReadOutNoiseSeed"] = readoutNoiseSeed + telescopeIndex
            sim["RandomSeeds/PhotonNoiseSeed"] = photonNoiseSeed + telescopeIndex 
            sim["RandomSeeds/JitterSeed"] = jitterSeed 
            sim["RandomSeeds/FlatFieldSeed"] = flatfieldSeed + telescopeIndex 
            sim["RandomSeeds/CTESeed"] = cteSeed + telescopeIndex 
            sim["RandomSeeds/DriftSeed"] = driftSeed + telescopeIndex  
        
            simFile = sim.run()
            print "Done"
            
        else:
            print "Sub-field centred on (" + str(raCenter) + ", " + str(decCenter) + ") does not lay entirely on a CCD for telescope " + str(telescope + 1) + " of group " + str(group + 1)
