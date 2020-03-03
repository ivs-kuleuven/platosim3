"""

Usage: $ python3 demo_multiTelescopeConfiguration.py

Example script showing how to run the Plato Simulator for a multi-telescope configuration.  The 32 telescopes are arranged in four groups of 8 telescopes.
All telescopes in the same group have the same FOV and the lines of sight of the four groups are offset by an angle of 9.2 degrees from the PLM z-axis.

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
from referenceFrames import sunSkyCoordinatesAwayfromPlatformPointing
from referenceFrames import CCD




################
# Input & Output
################

inputDir = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"
inputFile   = inputDir + "/inputfile.yaml"

outputDir   = os.getenv("PLATO_WORKDIR")
outputPrefix = "MultiTelescopeConfiguration"





##########################
# Configuration parameters
##########################



# Telescope parameters
######################

driftSeed = 1433429158                  # Random seed for telescope drift                                           <--- different for each telescope!


# CCD parameters
################

# Automatically determined:
#    - origin offset (x, y)
#    - orientation
#    - size

# Different for each telescope

readoutNoiseSeed =  1424949740          # Random seed for the readout noise                                         <--- different for each telescope!
photonNoiseSeed = 1433320336            # Random seed for the photon noise                                          <--- different for each telescope!
flatfieldSeed = 1433320381              # Random seed for the flatfield                                             <--- different for each telescope!


# Sub-field parameters
######################

# Zeropoint depends on the pointing of the telescope

numColumnsSubField = 40                 # Number of columns in the modelled sub-field [pixels]
numRowsSubField = 40                    # Number of rows in the modelled sub-field [pixels]





############
# Simulation
############

numTelescopeGroups = 4
numTelescopesPerGroup = 8

# Loop over all groups of telescopes

for group in range(numTelescopeGroups):
    
    # Loop over all telescopes in the current group
    
    for telescope in range(numTelescopesPerGroup):
        
        print("Processing telescope {0} of group {1}".format(telescope + 1, group + 1))
        
        telescopeIndex = numTelescopesPerGroup * group + telescope
        
        # Output will be stored in /MultiTelescopeConfiguration_group<group>_telescope<telescope>
        
        outputFilePrefix = outputPrefix + "_group" +  "{0:02d}".format(group + 1) + "_telescope" + "{0:02d}".format(telescope + 1)
        sim = Simulation(outputFilePrefix, inputFile)
        sim.outputDir = outputDir
        
        azimuthAngles = sim["CameraGroups/AzimuthAngle"]
        tiltAngles = sim["CameraGroups/TiltAngle"]
       
        # Get the Platform pointing from the yaml file

        raPlatform = sim["ObservingParameters/RApointing"]        # [deg]
        decPlatform = sim["ObservingParameters/DecPointing"]      # [deg]

        # Compute the telescope pointing, based on the platform pointing, and the tilt and azimuth angle of the telescope
        
        raSun, decSun = sunSkyCoordinatesAwayfromPlatformPointing(math.radians(raPlatform), math.radians(decPlatform))
        raTelescope, decTelescope = platformToTelescopePointingCoordinates(math.radians(raPlatform), math.radians(decPlatform), raSun, decSun, math.radians(azimuthAngles[group]), math.radians(tiltAngles[group]))

        print("Platform pointing: {0}, {1}".format(raPlatform, decPlatform))
        print("Telescope pointing: {0}, {1}".format(math.degrees(raTelescope), math.degrees(decTelescope)))
        
        includeFieldDistortion = (sim["Camera/IncludeFieldDistortion"] == "no")  
        distortionCoefficients = sim["Camera/FieldDistortion/ConstantCoefficients"]
        focalLength = sim["Camera/FocalLength/ConstantValue"] * 1000                   # Focal length [mm]
        plateScale  = sim["Camera/PlateScale"]                                         # Plate scale [arcsec / micron]
        pixelSize   = sim["CCD/PixelSize"]                                             # Pixel size [micron / pixel]
       

        focalPlaneAngle = 45.0 + group * 90.0                                                       # [deg]
        solarPanelOrientation = math.radians(float(sim["Platform/SolarPanelOrientation"]))          # [rad]

        # Set the subfield to the same pointing as the platform
        # This alwas falls on a CCD, provided that you properly turn the focal plane (see below)

        raCenter = raPlatform
        decCenter = decPlatform

        # Determine on which CCD (A, B, C, or D) the coordinates (raCenter, decCenter) are positioned and at which location
        # (in pixel coordinates)
        
        ccdCode, columnCenter, rowCenter = getCCDandPixelCoordinates(math.radians(raCenter), math.radians(decCenter), math.radians(raPlatform), math.radians(decPlatform), solarPanelOrientation, math.radians(tiltAngles[group]), math.radians(azimuthAngles[group]), math.radians(focalPlaneAngle), focalLength, pixelSize, includeFieldDistortion, distortionCoefficients, normal = True)
        
        # Check whether the sub-field falls entirely on the detector
        
        if (ccdCode != None) and (rowCenter - numRowsSubField / 2 >= 0) and (rowCenter + numRowsSubField / 2 < CCD[ccdCode]["Nrows"]) and (columnCenter - numColumnsSubField / 2 >= 0) and (columnCenter + numColumnsSubField / 2 < CCD[ccdCode]["Ncols"]):
            
            print("CCD #{0} selected".format(ccdCode) +"\n")
            
            # Observing parameters
        
            sim["ObservingParameters/RApointing"] = raPlatform
            sim["ObservingParameters/DecPointing"] = decPlatform

            # Camera & Telescope parameters
            
            sim["Telescope/GroupID"] = group + 1
            sim["Camera/FocalPlaneOrientation/ConstantValue"] = focalPlaneAngle
    
            # CCD parameters
            
            sim["CCD/Position"] =  ccdCode
            
            # Sub-field parameters
        
            sim["SubField/ZeroPointRow"] = int(rowCenter - numRowsSubField / 2)
            sim["SubField/ZeroPointColumn"] = int(columnCenter - numColumnsSubField / 2)
        
            sim["SubField/NumColumns"] = numColumnsSubField
            sim["SubField/NumRows"] = numRowsSubField
        
            # Seed parameters
        
            sim["RandomSeeds/ReadOutNoiseSeed"] = readoutNoiseSeed + telescopeIndex
            sim["RandomSeeds/PhotonNoiseSeed"] = photonNoiseSeed + telescopeIndex
            sim["RandomSeeds/FlatFieldSeed"] = flatfieldSeed + telescopeIndex 
            sim["RandomSeeds/DriftSeed"] = driftSeed + telescopeIndex  
        
            simFile = sim.run()
            
        else:
            print("Sub-field centred on ({0}, {1}) does not lay entirely on a CCD for telescope {2} of group {3}".format(raCenter, decCenter, telescope + 1, group + 1))
