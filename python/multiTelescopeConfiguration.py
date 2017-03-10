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

raPointing = 180.0                      # Platform right ascension pointing coordinate [degrees]
decPointing = -70.0                     # Platform declination pointing coordinate [degrees]


raCenter = 227.747586623                # Right ascension on which to centre the sub-field [degrees]
decCenter = -63.9874188222              # Declination on which to centre the sub-field [degrees]

# Telescope parameters
######################

azimuthAngles = [45, 135, -135, -45]    # Azimuth angles of the telescopes (same within a group) [degrees]
tiltAngle = 9.2                         # Tilt angle of the telescopes (same for the 4 groups) [degrees]

driftSeed = 1433429158                  # Random seed for telescope drift                                           <--- different for each telescope!


# CCD parameters
################

# Automatically determined:
#    - origin offset (x, y)
#    - orientation
#    - size

readoutNoiseSeed =  1424949740          # Random seed for the readout noise                                         <--- different for each telescope!
photonNoiseSeed = 1433320336            # Random seed for the photon noise                                          <--- different for each telescope!
flatfieldSeed = 1433320381              # Random seed for the flatfield                                             <--- different for each telescope!
cteSeed = 1424949740                    # Random seed for the CTE                                                   <--- different for each telescope!


# Sub-field parameters
######################

# Zeropoint depends on the pointing of the telescope

numColumnsSubField = 10                 # Number of columns in the modelled sub-field [pixels]
numRowsSubField = 10                    # Number of rows in the modelled sub-field [pixels]





############
# Simulation
############

numTelescopeGroups = 4
numTelescopesPerGroup = 8

# Loop over all groups of telescopes

for group in range(numTelescopeGroups):
    
    # Loop over all telescopes in the current group
    
    for telescope in range(numTelescopesPerGroup):
        
        print "Processing telescope " + str(telescope + 1) + " of group " + str(group + 1)
        
        telescopeIndex = numTelescopesPerGroup * group + telescope
        
        # Output will be stored in /MultiTelescopeConfiguration_group<group>_telescope<telescope>
        
        outputFilePrefix = outputPrefix + "_group" +  "{0:02d}".format(group + 1) + "_telescope" + "{0:02d}".format(telescope + 1)
        sim = Simulation(outputFilePrefix, inputFile)
        sim.outputDir = outputDir
        
        # Compute the telescope pointing, based on the platform pointing, and the tilt and azimuth angle of the telescope
        
        raTelescope, decTelescope = platformToTelescopePointingCoordinates(math.radians(raPointing), math.radians(decPointing), math.radians(azimuthAngles[group]), math.radians(tiltAngle))

        
        print "Platform pointing: " + str(raPointing) + ", " + str(decPointing)
        print "Telescope pointing: " + str(math.degrees(raTelescope)) + ", " + str(math.degrees(decTelescope))
        
        includeFieldDistortion = (str(sim["Camera/IncludeFieldDistortion"] == "yes"))  or (str(sim["Camera/IncludeFieldDistortion"] == "1"))        # Whether or not to include field distortion
        focalPlaneAngle = sim["Camera/FocalPlaneOrientation"]                                                                                       # Focal-plane orientation [degrees]
        focalLength = sim["Camera/FocalLength"] * 1000                                                                                              # Focal length [mm]
        plateScale = sim["Camera/PlateScale"]                                                                                                       # Plate scale [arcsec / micron]
        pixelSize = sim["CCD/PixelSize"]                                                                                                            # Pixel size [micron / pixel]
        
        # Determine on which CCD (A, B, C, or D) the coordinates (raCenter, decCenter) are positioned and at which location
        # (in pixel coordinates)
        
        ccdCode, columnCenter, rowCenter = getCCDandPixelCoordinates(math.radians(raCenter), math.radians(decCenter), raTelescope, decTelescope, math.radians(focalPlaneAngle), focalLength, plateScale, pixelSize, includeFieldDistortion, nominal=True)
        
        # Check whether the sub-field falls entirely on the detector
        
        if (ccdCode != None) and (rowCenter - numRowsSubField / 2 >= 0) and (rowCenter + numRowsSubField / 2 < CCD[ccdCode]["Nrows"]) and (columnCenter - numColumnsSubField / 2 >= 0) and (columnCenter + numColumnsSubField / 2 < CCD[ccdCode]["Ncols"]):
            
            print "CCD " + ccdCode + " selected"
            
            # Observing parameters
        
            sim["ObservingParameters/RApointing"] = raPointing
            sim["ObservingParameters/DecPointing"] = decPointing
        
            # Telescope parameters
                        
            sim["Telescope/AzimuthAngle"] = azimuthAngles[group]
            sim["Telescope/TiltAngle"] = tiltAngle
    
            # CCD parameters
        
            sim["CCD/OriginOffsetX"] = CCD[ccdCode]["zeroPointXmm"]
            sim["CCD/OriginOffsetY"] = CCD[ccdCode]["zeroPointYmm"]
            sim["CCD/Orientation"] = CCD[ccdCode]["angle"]
            sim["CCD/NumColumns"] = CCD[ccdCode]["Ncols"]
            sim["CCD/NumRows"] = CCD[ccdCode]["Nrows"]
            
            # Sub-field parameters
        
            sim["SubField/ZeroPointRow"] = int(rowCenter - numRowsSubField / 2)
            sim["SubField/ZeroPointColumn"] = int(columnCenter - numColumnsSubField / 2)
        
            sim["SubField/NumColumns"] = numColumnsSubField
            sim["SubField/NumRows"] = numRowsSubField
        
            # Seed parameters
        
            sim["RandomSeeds/ReadOutNoiseSeed"] = readoutNoiseSeed + telescopeIndex
            sim["RandomSeeds/PhotonNoiseSeed"] = photonNoiseSeed + telescopeIndex
            sim["RandomSeeds/FlatFieldSeed"] = flatfieldSeed + telescopeIndex 
            sim["RandomSeeds/CTESeed"] = cteSeed + telescopeIndex 
            sim["RandomSeeds/DriftSeed"] = driftSeed + telescopeIndex  
        
            simFile = sim.run()
            
        else:
            print "Sub-field centred on (" + str(raCenter) + ", " + str(decCenter) + ") does not lay entirely on a CCD for telescope " + str(telescope + 1) + " of group " + str(group + 1)
