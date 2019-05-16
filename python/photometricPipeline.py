"""
Photometric pipeline.

Example usage:

    >>> from photometricPipeline import *
    >>> photometricPipeline = PhotometricPipeline(runName, configurationFile, simulationFile, outputDir)
    >>> photometricPipeline.run()
with
    - runName: Arbitrary name for this run of the photometric pipeline (will be used as name for the output HDF5 file)
    - configurationFile: Configuration file for the photometric pipeline (e.g. /inputfiles/photometricPipeline.yaml)
    - simulationFile: Filename of the HDF5 output file of PlatoSim
    - outputDir: Output directory
"""

#########
# Imports
#########

from simfile import SimFile
import numpy as np
import math
import yaml
import pyaml 
import inspect
import subprocess
import os
import sys
import ast
import numpy as np
import h5py
import astropy
from analyticpsf import AnalyticPSF






#####################
# Auxiliary functions
#####################

def applyShiftedDataAlgorithm(data):

    """
    PURPOSE: Shifted data algorithm to compute the mean and variance for the given data.  This method avoids 
             loss of significance with big numbers when computing the variance.

    INPUT:
        - data: Data array for which to calculate the mean and the variance

    OUTPUT:
        - mean: Mean value of the given data array, as computed with the shifted data algorithm
        - variance: Variance of the given data array, as computed with the hifted data algorithm
        - numUsefulDatapoints: Number of datapoints
    """

    # Discard flagged data (if any)

    approximationOfMean = data[0]                            # Use the 1st value as approximation of the mean (arbitrary choice)
    data = data - approximationOfMean                        # Shift the data (by subtracting the approximation of the mean)

    sumShiftedData = np.sum(data)                            # Sum of the shifted data
    squaredSumShiftedData = np.sum(np.power(data, 2))        # Sum of the squares of the shifted data

    numUsefulDatapoints = len(data)                          # Number of datapoints

    mean = (sumShiftedData + numUsefulDatapoints * approximationOfMean) / numUsefulDatapoints                            # Mean
    variance = (squaredSumShiftedData - (pow(sumShiftedData, 2) / numUsefulDatapoints)) / (numUsefulDatapoints - 1)      # Variance

    return mean, variance, numUsefulDatapoints





def applyMadClippingAroundMedian(data, threshold, index = None):

    """
    PURPOSE: Outlier detection, based on Median Absolute Deviation (MAD) clipping around the median. In
             case an index is specified, a boolean is returned, indicating whether or not the datapoint
             with the specified index is an outlier.  Otherwise an outlier flag for the given data is 
             returned.
    INPUT:
        - data: Data for which to flag outliers
        - threshold: Threshold for MAD clipping around the median
        - index: If specified, it is calculated whether the datapoint at this position is an outlier

    OUTPUT:
        - in case index was specified: boolean indicating whether or not the datapoint at the given
          position is an outlier
        - otherwise: boolean truth values, where 0 means "no outlier" and 1 means "outlier".  The truth 
          value 1 is equivalent to a flag.
    """

    mad = astropy.stats.median_absolute_deviation(data)       # MAD

    if index == None:

        return np.where(data > threshold * mad)               # Outlier flag for all datapoints

    else:

        return data[index] / mad > threshold                  # Is the datapoint at the given index an outlier?










class PhotometricPipeline(object):
    
    def __init__(self, runName, configurationFile, simulationFile, outputDir = None):

        """
        PURPOSE: Initialisation of the class variables and readong of the default input files.

        INPUT:
            - runName: Name of the run of the photometric pipeline (will be used as name for the output file)
            - configurationFile: Name of the configuration file to run the photometric pipeline
            - simulationFile: Name of the output file of PlatoSim (this file contains the serial pre-scan,
                              parallel over-scan, and fluxes in a simulated sub-field on the CCD)
            - outputDir: Directory in which the results of the photometric pipeline will be stored
        """

        # Run name

        self.runName = runName

        # Configuration file

        if configurationFile != None:

            self.readConfigurationFile(configurationFile)
        
        else:

            self.readConfigurationFile(os.environ["PLATO_PROJECT_HOME"] + "/inputfiles/photometricPipeline.yaml")

        # Simulation file

        self.simFile = SimFile(simulationFile)

        # Output directory
        # (flag to check if output directory has been specified before running the photometric pipeline)

        self.hasTargetLocation = False
        self.targetOutputFilesLocation = None

        if outputDir != None:
            self.outputDir = outputDir
        




    @property
    def outputDir(self):

        """
        PURPOSE: Return the location of the output file.

        OUTPUT: 
            - Location of the output file
        """
        return self.targetOutputFilesLocation





    @outputDir.setter
    def outputDir(self, path):

        """
        PURPOSE: Specify the absolute path for the output directory. This directory will contain
                 the HDF5 output file for this run of the photometric pipeline.  If the output path 
                 doesn't exist, it is created.
        
        INPUT: 
            - path: Path to the directory to use as output directory
        """

        if not os.path.exists(path):

            self.createDirectory(path)
        
        self.targetOutputFilesLocation = path
        self.hasTargetLocation = True





    def readConfigurationFile(self, filename):

        """
        PURPOSE: Read the YAML input configuration file. 

        INPUT:
            - filename: Path to the configuration file
        """

        self.configurationFilename = filename
        
        with open(filename, 'r') as stream:

            try:

                self.yamlDocument = yaml.load(stream)

            except yaml.YAMLError as exc:
                
                print(exc)





    def getYamlConfiguration(self):

        """
        PURPOSE: Return the YAML configuration as a dictionary.

        OUTPUT:
            - YAML configuration as a dictionary
        """

        return self.yamlDocument





    def __contains__(self, key):

        """
        PURPOSE: Return True if the input parameter (key) is known/exists; False otherwise.

        INPUT:
            - key: string containing the parameter name or "Group/ParameterName" combination

        OUTPUT:
            - Boolean indicating whether or not the given key exists in the configuration file
        """

        if key.find('/') == -1:
        
            nodeNames = [key]
        
        else:

            nodeNames = key.split("/")

        node = self.yamlDocument

        for nodeName in nodeNames:

            print("> {}, {}".format(nodeName, type(node)))

            try:

                node = node[nodeName]
            
            except:
            
                return False

        return True





    def __getitem__(self, key):

        """
        PURPOSE: Return the value of the input parameter (key).

        INPUT: 
            - key: a string containing the parameter name or "Group/ParameterName" combination

        OUTPUT:
            - Value of the parameter
        """
        
        # Split the path into node names
        # E.g. "PSF/MappedGaussian/Sigma" into ["PSF", "MappedGaussian", "Sigma"]

        if key.find('/') == -1:

            # parentNodeName, nodeName = key, None
            print ("usage: the given parameter name (key) should include the group name of the group that contains the parameter.")
            print ("       E.g in 'Camera/PlateScale', Camera is the group, PlateScale is the parameter.")

            return None

        else:

            nodeNames = key.split("/")

        # Navigate to the deepest node, starting from the document root

        node = self.yamlDocument 

        for nodeName in nodeNames:

            if nodeName in node:

                node = node[nodeName]
            
            else:

                print("ERROR: The group '{}' was not found in the yaml inputfile '{}'.".format(key, self.configurationFilename))

                return None

        # node is a string, so cast it to its proper value

        try:

            value = ast.literal_eval(node)

        except ValueError:
            
            value = node
        
        # Return the value of the deepest node

        return value





    def __setitem__(self, key, item):
        
        """
        PURPOSE: Update a specific node.

        INPUT:
            - key: string with parent node name and node name seperated by a slash
            - item: string with the new node value, if not a string the value is converted using str()

        OUTPUT: 
            - True if node could be updated; False otherwise
        """
        
        # Ensure that the given item is a string

        item = str(item)

        # Split the path into node names
        # E.g. "PSF/MappedGaussian/Sigma" into ["PSF", "MappedGaussian", "Sigma"]

        if key.find('/') == -1:

            print ("usage: the given parameter name (key) should include the group name of the group that contains the parameter.")
            print ("       E.g in 'Camera/PlateScale', Camera is the group, PlatScale is the parameter.")

            return None
        
        else:

            nodeNames = key.split("/")

        # Check whether the parent node is in the document. If not, complain

        if nodeNames[0] not in self.yamlDocument:

             print("Error: no node with the name {0} found in input yaml file".format(nodeNames[0]))

             return False

        # If there is only 1 node in the path, we're finished after setting its value

        if len(nodeNames) == 1:

            self.yamlDocument[nodeNames[0]] = item

            return True

        # If we arrive here, there are at least 2 node in the path, check if 2nd parent node exists

        if nodeNames[1] not in self.yamlDocument[nodeNames[0]]:

             print("Error: no node with the name {0} found in input yaml file".format(nodeNames[0]+"/"+nodeNames[1]))

             return False

        # If there are only 2 nodes in the path, we're finished after setting its value

        if len(nodeNames) == 2:

            self.yamlDocument[nodeNames[0]][nodeNames[1]] = item

            return True

        # If we arrive here, there are at least 3 nodes in the path, check if 3rd parent node exists

        if nodeNames[2] not in self.yamlDocument[nodeNames[0]][nodeNames[1]]:

             print("Error: no node with the name {0} found in input yaml file".format(nodeNames[0]+"/"+nodeNames[1]+"/"+nodeNames[2]))

             return False

        # If there are only 3 nodes in the path, we're finished after setting its value

        if len(nodeNames) == 3:

            self.yamlDocument[nodeNames[0]][nodeNames[1]][nodeNames[2]] = item

            return True

        print("Error: detected more than 3 nodes in the path {0}".format(key))

        return False





    def createDirectory(self, path):

        """
        PURPOSE: Create the given directory.

        INPUT:
            - path: Directory to create
        """

        try: 
            os.makedirs(path)

        except OSError as ose:

            print (ose)

            if not os.path.isdir(path):
                raise Exception("Couldn't create directory {}".format(path))

        return





    def writeYamlConfigurationFile(self, filename):

        """
        PURPOSE: Write the modified configuration to output file location. This configuration will 
                 be loaded by the photometric pipeline when the run() method is executed.
        
        INPUT:
            - filename: Filename
        """

        with open(filename, 'w') as outfile:

            outfile.write( pyaml.dump(self.yamlDocument, indent=4, width=120) )





    def __str__(self):

        """
        PURPOSE: Return a listing of all settings from the photometric pipeline.  If a parameter value has been updated 
                 for this run, [updated] will be printed after the value.
        """

        root = self.yamlDocument
        msg = "YAML Configuration:\n"
        msg += pyaml.dump(root, indent=4)

        return msg





    def run(self, removeOutputFile = False):

        """
        PURPOSE: Run the photometric pipeline.

        INPUT:
            - removeOutputFile: Boolean indicating whether or not the output file should be delete in case it
                                already exists before the run started
        """

        if not self.hasTargetLocation:

            raise Exception("Output location not set for this Simulation. Set the outputDir before executing the run() method.")

        inputFilename = "{}/{}.yaml".format(self.targetOutputFilesLocation, self.runName)
        outputFilename = "{}/{}.hdf5".format(self.targetOutputFilesLocation, self.runName)

        if removeOutputFile:
        
            try:

                os.remove(outputFilename)
            
            except OSError:
            
                pass
        
        self.writeYamlConfigurationFile(inputFilename)



        # Configuration:
        #   - read the required input parameters
        #   - make placeholders for all arrays that will be stored

        self.beginExposureNr = self.simFile.getInputParameter("ObservingParameters/", "BeginExposureNr")
        self.numExposures = self.simFile.getInputParameter("ObservingParameters", "NumExposures")

        self.configure()


        self.outputFile = h5py.File(outputFilename,"w")

        # Process all exposures

        photometryGroup = self.outputFile.create_group("/Photometry");
        time = np.array(self.simFile.getTime())
        photometryGroup.create_dataset("Time", data=time)


        masks = []

        for exposure in range(self.beginExposureNr, self.beginExposureNr + self.numExposures):

            exposureGroupName = "/Photometry/Exposure{0:06d}".format(exposure)
            exposureGroup = self.outputFile.create_group(exposureGroupName) 

            masks = self.processExposure(exposure, exposureGroup, masks)

        

        # Write results to HDF5 file

        self.writeOutput()
        self.outputFile.close()
    




    def configure(self):

        """
        PURPOSE: Configuration:
                    - read all required input parameters
                    - make placeholders for all arrays that will be stored
        """


        # self.time = self.simFile["/StarPositions/Time"]
        self.numRowsSubField = self.simFile.getInputParameter("SubField", "NumRows")
        self.numColumnsSubField = self.simFile.getInputParameter("SubField", "NumColumns")
        self.flatfield = self.simFile.getPRNU()
        self.varianceReadoutNoise = math.pow(self.simFile.getInputParameter("FEE", "ReadoutNoise"), 2) + math.pow(self.simFile.getInputParameter("CCD", "ReadoutNoise"), 2)

        # Parameters that are specific for the detector half on which the star windows are located
        # (we assume the whole sub-field is on the same detector half)

        self.configureDetectorHalf()

        # Parameters that are specific for the offset calculation

        self.configureOffsetCalculation()

        # Parameters that are specific for the smearing pattern calculation
        
        self.configureSmearingPatternCalculation()

        # Parameters that are specific for the sky background

        self.configureSkyBackground()

        # Parameters that are specific for the mask

        self.configureMaskUpdate()

        # Parameters that are specific for the star windows

        self.configureStarWindows()

        # Parameters that are specific for the flux and COB calculation

        self.configureFluxAndCobCalculation()

        # # Parameters that are specific for light curve outlier detection

        # self.configureLightCurveOutlierDetection()

        # # Parameters that are specific for time averaging
        
        # self.configureTimeAveraging()





    def configureDetectorHalf(self):

        """
        PURPOSE: Configuration of the parameters that are specific for the detector half 
                 on which the star windows are located.  Assumed is that the whole sub-field
                 is on the same detector half.
        """

        self.subfieldZeropointRow = self.simFile.getInputParameter("SubField", "ZeroPointRow")           # Row coordinate of the sub-field zeropoint on the detector [pixels]
        self.subfieldZeropointColumn = self.simFile.getInputParameter("SubField", "ZeroPointColumn")     # Column coordinate of the sub-field zeropoint on the detector [pixels]

        

        if self.subfieldZeropointColumn >= self.simFile.getInputParameter("CCD", "NumColumns") / 2:
        
            self.detectorHalf = "Right"
            self.gain = 1.0 / (self.simFile.getInputParameter("FEE/Gain", "RefValueRight") * self.simFile.getInputParameter("CCD/Gain", "RefValueRight"))  # Total gain for the left detector half [e- / ADU]
        
        else:
            
            self.detectorHalf = "Left"
            self.gain = 1.0 / (self.simFile.getInputParameter("FEE/Gain", "RefValueLeft") * self.simFile.getInputParameter("CCD/Gain", "RefValueLeft"))    # Total gain for the left detector half [e- / ADU]
            




    def configureOffsetCalculation(self):

        """
        PURPOSE: Configuration of the parameters that are specific for the offset
                 calculation (and the corresponding outlier detection) and making 
                 placeholders for the information that will be stored.
        """
        
        self.includeOffsetOutlierDetection = self["Offset/IncludeOutlierDetection"]                     # Enable/disable outlier detection
        self.offsetOutlierDetectionNumSkippedElementsBothEnds = self["Offset/OutlierDetection/k"]       # Number of largest and smallest values to flag as outliers

        if self.includeLightCurveOutlierDetection:

            numBiasElements = self.simFile.getInputParameter("SubField", "NumBiasPrescanRows") * self.simFile.getInputParameter("SubField", "NumBiasPrescanColumns")

            if 2 * self.offsetOutlierDetectionNumSkippedElementsBothEnds <= numBiasElements:

                raise Exception("For the outlier detection in the offset calculation, the number of skipped pixels must be larger than the number of pixels in the serial pre-scan (i.e. bias register map)")
        
        self.offsetValueArrayFastCadence = np.array([])
        # self.offsetVarianceArrayFastCadence_array = np.array([])
    




    def configureSmearingPatternCalculation(self):

        """
        PURPOSE: Configuration of the parameters that are specific for the smearing
                 pattern calculation (and the corresponding outlier detection) and
                 making placeholders for the information that will be stored.
        """

        self.smearingCoefficientA0Array = np.empty(self.simFile.getInputParameter("SubField", "NumColumns")) 
        self.smearingCoefficientA0Array.fill(np.array(self["Smearing/Coefficients/a"])[0])                      # Coefficient a0 (one entry per column)
        self.smearingCoefficientsA = np.array(self["Smearing/Coefficients/a"])[1:]                              # Coefficients [a1, a2, a3]
        self.smearingCoefficientsB = np.array(self["Smearing/Coefficients/b"])                                  # Coefficients [b0, b1, b2, b3]
        self.smearingNumRowsSkipped = self["Smearing/NumRowsSkipped"]                                           # Number of rows to skip for CTI correction (1st rows may be affected by bright stars at the top of the detector) 
        self.includeSmearingOutlierDetection = self["Smearing/IncludeOutlierDetection"]                         # Enable/disable outlier detection
        self.smearingOutlierDetectionThreshold = self["Smearing/OutlierDetection/Threshold"]                    # Threshold for outlier detection
        self.smearingRegularization = self["Smearing/Regularization"]                                           # Epsilon regularisation parameter

        self.smearingPatternArrayFastCadence = np.array([])                                                     # Smearing pattern for the different columns at fast cadence (25s)
        self.smearingPatternArrayLongCadence = np.array([])                                                     # Smearing pattern for the difference columns at long cadence (600s)

        if self.smearingNumRowsSkipped <= self.simFile.getInputParameter("SubField", "NumSmearingOverscanRows"):

            raise Exception("You cannot skip all rows in the parallel over-scan (i.e. smearing map)")





    def configureSkyBackground(self):

        """
        PURPOSE: Making placeholders for the information that will be stored
                 when calculating the sky background.
        """

        self.skyBackgroundArrayFastCadence = np.array([])

    



    def configureStarWindows(self):

        """
        PURPOSE: Configuration of the parameters that are specific for the star windows.
        """

        self.contaminationRadius = self["StarWindows/ContaminationRadius"]
        self.contaminantIds = self["StarWindows/ContaminantIds"]
        self.hasKnownContaminants = (len(self.contaminantIds) > 0)
        self.windowDimensions = self["StarWindows/WindowDimension"]
        self.windowOffset = self.windowDimensions // 2 + 1

        self.targetIds = self["StarWindows/TargetIds"]
        self.numTargets = len(self.targetIds)

        self.targetRowsArray = np.empty(self.numTargets)
        self.targetColumnsArray = np.empty(self.numTargets)





    def configureMaskUpdate(self):

        """
        PURPOSE: Configuration of the parameters that are specific for updating
                 the (nominal) mask.
        """

        self.numExposuresMaskUpdate = self["StarWindows/MaskUpdateInterval"] * 24 * 60 * 60 / 25         # Number of exposures after which the mask needs to be updated 
        self.numMaskUpdates = math.ceil(self.numExposures / self.numExposuresMaskUpdate)

        path                = "../inputfiles/psfallv3.txt"                                                                  # File comprising the parameters characterising the analytic non-Gaussian PSF
        sigmaPSF            = self.simFile.getInputParameter("PSF/AnalyticNonGaussian/Sigma", "ConstantValue")              # Width of the analytic non-Gaussian PSF [pixels]
        sigmaDiffusion      = self.simFile.getInputParameter("PSF/AnalyticNonGaussian", "ChargeDiffusionStrength")          # Width of Gaussian diffusion kernel, modelling the charge diffusion [pixels]
        focalLength         = self.simFile.getInputParameter("Camera/FocalLength", "ConstantValue") * 1000                  # Focal length [mm]
        ccdOrientation      = self.simFile.getInputParameter("CCD", "Orientation") / 180 * np.pi                            # Orientation angle of the CCD w.r.t. the orientation of the focal plane [radians]
        ccdZeropointRow     = self.simFile.getInputParameter("CCD", "OriginOffsetY")                                        # Row coordinate of the sub-field origin on the detector[mm]
        ccdZeropointColumn  = self.simFile.getInputParameter("CCD", "OriginOffsetX")                                        # Column coordinate of the sub-field origin on the detector[mm]
        pixelSize           = self.simFile.getInputParameter("CCD", "PixelSize")                                            # Pixel size [micron]
        
        self.analyticPSF = AnalyticPSF(path, sigmaPSF, sigmaDiffusion, focalLength, ccdOrientation, ccdZeropointColumn, ccdZeropointRow, pixelSize)

        



    def configureFluxAndCobCalculation(self):

        """
        PURPOSE: Making placeholders for the information that will be stored when 
                 calculating the flux and the COB.
        """

        self.cobOffsetX = np.arange(0.5, self.windowDimensions).repeat(self.windowDimensions).reshape(self.windowDimensions, self.windowDimensions).transpose()   # Matrix to account for the pixel offset in the x-direction when calculating the COB
        self.cobOffsetY = np.arange(0.5, self.windowDimensions).repeat(self.windowDimensions).reshape(self.windowDimensions, self.windowDimensions)               # Matrix to account for the pixel offset in the y-direction when calculating the COB




        
    # def configureLightCurveOutlierDetection(self):

    #     """
    #     PURPOSE: Configuration of the parameters that are specific for the light curve
    #              outlier detection and making placeholders for the information that will
    #              be stored.
    #     """

    #     self.includeLightCurveOutlierDetection = self["LightCurve/IncludeOutlierDetection"]          # Enable/disable outlier detection
    #     self.fluxOutlierDetectionHalfBinWidth = self["LightCurve/OutlierDetection/b"]                # Number of past and future datapoints needed to decide whether or not a datapoint is an outlier
    #     self.lightCurveOutlierDetectionThreshold = self["LightCurve/OutlierDetection/Threshold"]     # Outlier detection threshold


    #     if self.includeLightCurveOutlierDetection:

    #         self.fluxFlagFastCadence = np.zeros(self.fluxOutlierDetectionHalfBinWidth)
        
    #     else:

    #         self.fluxFlagFastCadence = None





    # def configureTimeAveraging(self):

    #     """
    #     PURPOSE: Configuration of the parameters that are specific for time averaging
    #              and providing placeholders for the information that will be stored.
    #     """

    #     self.timeAveragingCadence = self["LightCurve/TimeAveraging/Cadence"]                  # Choose between short and long cadence

    #     self.numExposuresShortCadence = self["LightCurve/TimeAveraging/NumSamples/Short"]     # Short cadence: 50s (i.e. 2 samples)
    #     self.numExposuresLongCadence = self["LightCurve/TimeAveraging/NumSamples/Long"]       # Long cadence: 600s (i.e. 24 samples)

    #     if self.timeAveragingCadence == "Short":
            
    #         self.numExposuresTimeAveraging = self.numExposuresShortCadence

    #         self.fluxArrayShortCadence = np.array([])
    #         self.cobArrayShortCadence = np.array([])
    #         self.fluxOutlierFlagShortCadence = np.array([])

    #     elif self.timeAveragingCadence == "Long":
            
    #         self.numExposuresTimeAveraging = self.numExposuresShortCadence

    #         self.fluxArrayLongCadence = np.array([])
    #         self.fluxVarianceLongCadence = np.array([])
    #         self.fluxOutlierFlagLongCadence = np.array([])





    def processExposure(self, exposure, exposureGroup, masks):

        """
        PURPOSE: Process the exposure with the given index.

        INPUT:
            - exposure: Index of the exposure to process
            - exposureGroup: Group in the output HDF5 file in which the flux and COB will be stored
            - masks: Masks for all target stars, before processing the given exposure

        OUTPUT:
            - Masks for all target stars, after processing the given exposure
        """

        # Calculate the offset for the current exposure (fast cadence)

        if self.detectorHalf == "Left":
            serialPreScan = self.simFile.getBiasMapLeft(exposure)

        elif self.detectorHalf == "Right":
            serialPreScan = self.simFile.getBiasMapRight(exposure)

        offsetValueFastCadence = self.calculateOffset(serialPreScan)[0]
        self.offsetValueArrayFastCadence = np.append(self.offsetValueArrayFastCadence, offsetValueFastCadence)



        # Calculate the smearing pattern for the current exposure (fast cadence)

        parallelOverScan = self.simFile.getSmearingMap(exposure)
        stdDevPrevious = np.empty(parallelOverScan.shape[1])
        stdDevPrevious.fill(9999)

        smearingPatternFastCadence, stdDevPrevious = self.calculateSmearing(parallelOverScan, offsetValueFastCadence, stdDevPrevious)
        self.smearingPatternArrayFastCadence = np.append(self.smearingPatternArrayFastCadence, smearingPatternFastCadence)



        # Calculate the background for the current exposure (fast cadence)
        # Real photometric pipeline: ~100 background windows (nominally 4x4) per CCD half,
        # for which we know the fluxes and the position on the CCD.  The latter is needed to 
        # subtract the correct smearing pattern from each background window.  For each of the
        # background windows the mean and variance are Calculated, which are used in an
        # interpolation schema (based on radial basis functions) to determine the background
        # at the position of the target star.
        # It seems unfeasible to make (long-term) simulations for all these background windows, so
        # we will have to look for an alternative way to estimate the background at the target
        # location.


        skyBackgroundFastCadence = self.simFile.getSkyBackground(exposure - self.beginExposureNr)                       # [Photons]
        self.skyBackgroundArrayFastCadence = np.append(self.skyBackgroundArrayFastCadence, skyBackgroundFastCadence)
        

        image, throughput = self.getStarPhotons(exposure)



        

        # Update the mask when required

        if exposure % self.numExposuresMaskUpdate == 0:

            masks = self.calculateMask(exposure, image, throughput)


        # Calculate the flux & COB (nominal & extended mask) for the current exposure
        # -> add new datapoint to the light curve (fast cadence)

        fluxArrayFastCadence = np.array([])
        cobRowArrayFastCadence = np.array([])
        cobColumnArrayFastCadence = np.array([])

        for targetIndex in range(self.numTargets):

            fluxFastCadence, cobFastCadence = self.calculateFluxAndCob(exposure, offsetValueFastCadence, smearingPatternFastCadence, targetIndex, masks[targetIndex])
            fluxArrayFastCadence = np.append(fluxArrayFastCadence, fluxFastCadence)
            cobRowArrayFastCadence = np.append(cobRowArrayFastCadence, cobFastCadence[0])
            cobColumnArrayFastCadence = np.array(cobColumnArrayFastCadence, cobFastCadence[1])

        starIds, starRows, starColumns, inputFlux = (self.simFile.getStarCoordinates(exposure)[index] for index in [0, 1, 2, 5])
        inputFlux = inputFlux[starIds == self.targetIds]
        targetRows = starRows[starIds == self.targetIds]
        targetColumns = starColumns[starIds == self.targetIds]

        for targetIndex in self.numTargets:

            targetRowAsInt = int(targetRows[targetIndex])
            targetColumnAsInt = int(targetColumns[targetIndex])
            inputFlux[targetIndex] *= throughput[targetRowAsInt, targetColumnAsInt] * self.flatfield[targetRowAsInt, targetColumnAsInt]


        exposureGroup.create_dataset("FluxFastCadence", data = fluxArrayFastCadence)
        exposureGroup.create_dataset("InputFlux", data = inputFlux)
        exposureGroup.create_dataset("CobRowFastCadence", data = cobRowArrayFastCadence)
        exposureGroup.create_dataset("CobColumnFastCadence", data = cobColumnArrayFastCadence)
        



        # # Smearing time averaging (fast -> long cadence)

        # if (exposure + 1) % self.numExposuresLongCadence == 0:

        #     smearingPatternLongCadence = self.timeAverageSmearing()
        #     self.smearingPatternArrayLongCadence = np.append(self.smearingPatternArrayLongCadence, smearingPatternLongCadence)



        # # Light curve outlier detection
        # # Note that you need b datapoints in the past and b datapoints in the future!

        # if self.includeLightCurveOutlierDetection:

        #     if (exposure + 1) > 2 * self.fluxFlagFastCadence:
            
        #         self.fluxFlagFastCadence = self.detectFluxOutliers(self.fluxArrayFastCadence, self.fluxFlagFastCadence)
        #         # TODO What at the end of the time series (when there are no future datapoints)?



        # # Flux & COB time averaging

        # if (exposure + 1) % self.numExposuresTimeAveraging == 0:

        #     self.timeAverageFluxAndCob()

        return masks





    def writeOutput(self):

        """
        PURPOSE: Write results of the photometric pipeline to the HDF5 output file.

        INPUT:
            - filename: Name of the output file
        """

        # Offset

        offsetGroup = self.outputFile.create_group("Offset")
        offsetGroup.create_dataset("OffsetValueFastCadence", dtype = "float32", data = self.offsetValueArrayFastCadence)
        # offsetGroup.create_dataset("offsetVarianceFastCadence", dtype = "float32", data = self.offsetVarianceArrayFastCadence)

        # Smearing pattern, fast cadence

        smearingGroup = self.outputFile.create_group("SmearingPattern")
        smearingGroup.create_dataset("SmearingPatternFastCadence", dtype = "float32", data = self.smearingPatternArrayFastCadence)

        # Sky background, fast cadence

        skyBackgroundGroup = self.outputFile.create_group("SkyBackground")
        skyBackgroundGroup.create_dataset("SkyBackgroundFastCadence", dtype = "float32", data = self.skyBackgroundArrayFastCadence)

        # # Smearing pattern, long cadence

        # smearingGroup.create_dataset("smearingPatternLongCadence", dtype = "float32", data = self.smearingPatternArrayLongCadence)

        # # Star window, fast cadence

        # starWindowGroupFastCadence = outputFile.create_group("starWindowFastCadence")
        # starWindowGroupFastCadence.create_dataset("fluxFastCadence", dtype = "float32", data = self.fluxArrayFastCadence)
        # starWindowGroupFastCadence.create_dataset("cobFastCadence", dtype = "float32", data = self.cobArrayFastCadence)

        # # Star window, short cadence

        # if self.timeAveragingCadence == "Short":

        #     starWindowGroupShortCadence = outputFile.create_group("starWindowShortCadence")
        #     starWindowGroupShortCadence.create_dataset("flux", dtype = "float32", data = self.fluxArrayShortCadence)
        #     starWindowGroupShortCadence.create_dataset("cob", dtype = "float32", data = self.cobArrayShortCadence)
        #     # starWindowGroupShortCadence.create_dataset("FX_EXPOSURE_ERROR_SC_ARRAY", dtype = "float32", data = self.fluxOutlierFlagShortCadence)

        
        # elif self.timeAveragingCadence == "Long":

        #     starWindowGroupLongCadence = outputFile.create_group("starWindowLongCadence")
        #     starWindowGroupLongCadence.create_dataset("flux", dtype = "float32", data = self.fluxArrayLongCadence)
        #     starWindowGroupLongCadence.create_dataset("fluxVariance", dtype = "float32", data = self.fluxVarianceArrayLongCadence)
        #     # starWindowGroupLongCadence.create_dataset("FX_EXPOSURE_ERRORLongCadence_ARRAY", dtype = "float32", data = self.fluxOutlierFlagLongCadence)





    ####################
    # Offset calculation
    ####################

    def calculateOffset(self, biasMap):

        """
        PURPOSE: Offset calculation as explained in PLATO-LESIA-PDC-DD-005 
                 (PLATO: N_DPU Onboards Offset Calculation ATBD).
                 Algorithm name: ONB-OFFCAL-010.

        INPUT:
            - biasMap: Serial pre-scan (i.e. bias register map) that is used to calculate the offset [ADU]
    
        OUTPUT:
            - offsetValueFastCadence: Mean of the values in the serial pre-scan after discarding the 
                                      outliers [ADU]
            - offsetVarianceFastCadence: Variance of the values in the serial pre-scan after discarding 
                                         the outliers [ADU]
        """

        # Outlier detection

        if self.includeOffsetOutlierDetection:

            # Flag: 
            #   - 0 means "no outlier"
            #   - 1 means "outlier"

            offsetFlag = self.detectOffsetOutliers(biasMap)
    
            # Shifted data algorithm

            offsetValueFastCadence, offsetVarianceFastCadence = applyShiftedDataAlgorithm(biasMap[~offsetFlag])[:2]
        


        # No outlier detection

        else:

            offsetValueFastCadence, offsetVarianceFastCadence = applyShiftedDataAlgorithm(biasMap)[:2]

        return offsetValueFastCadence, offsetVarianceFastCadence


    


    def detectOffsetOutliers(self, biasMap):
  
        """
        PURPOSE: Outlier detection in the serial pre-scan as explained in PLATO-MPSSR-PDC-DD-0002
                 (PLATO On-Board Offset & Prescan Outlier Detection Algorithm Theoretical Baseline
                 Document).
                 Algorithm name: ONB-OFFOUTDET-010.
    
        INPUT:
            - biasMap: Serial pre-scan that is used to calculate the offset [ADU]
    
        OUTPUT:
            - offsetFlag: Boolean truth values, where 0 means "no outlier" and 1 means
                          "outlier".  The truth value 1 is equivalent to a flag.
        """

        biasMap1d = np.ravel(biasMap)               # 2D -> 1D
        biasMap1dSorted = np.sort(biasMap1d)        # Sort

        if len(biasMap) < 2 * self.offsetOutlierDetectionNumSkippedElementsBothEnds:
            raise Exception("Number of entries, {0}, in the serial pre-scan (bias register map) should exceed 2k = {1}".format(len(biasMap1d), 2 * self.offsetOutlierDetectionNumSkippedElementsBothEnds))

        # Flag the offsetOutlierDetectionNumSkippedElementsBothEnds smallest values and the offsetOutlierDetectionNumSkippedElementsBothEnds largest values

        minOffset = biasMap1dSorted[self.offsetOutlierDetectionNumSkippedElementsBothEnds]
        maxOffset = biasMap1dSorted[-(self.offsetOutlierDetectionNumSkippedElementsBothEnds + 1)]

        offsetFlag = np.logical_or((biasMap < minOffset), (biasMap > maxOffset))

        return offsetFlag





    ######################
    # Smearing calculation
    ######################

    def calculateSmearing(self, smearingMap, offsetValueFastCadence, stdDevPrevious):

        """
        PURPOSE: Smearing calculation as explained in PLATO-LESIA-PDC-TN- (Parallel overscan rows:
                 correction of the CTI) and PLATO-LESIA-PDC-DD-006 (PLATO: N-DPU Onboard Smearing
                 Calculation ATBD). Coefficient a0 will be updated.
                 Algorithm name: ONB-SMRCAL-010.

        INPUT:
            - smearingMap: Parallel over-scan that is used to calculate the smearing [ADU]
            - offsetValueFastCadence: Electronic offset as calculated from the serial pre-scan [ADU]
            - stdDevPrevious: Standard deviation of the previous measurement of the parallel over-scan
                              (one entry per column of the parallel over-scan) [e-]

        OUTPUT:
            - smearingPatternFastCadence: One smearing row for this CCD half [e-]
            - stdDevPrevious: Standard deviation of the current measurement of this column of the parallel 
                              over-scan
        """

        numRowsSmearingMap = smearingMap.shape[0]
        smearingMap = (smearingMap - offsetValueFastCadence) * self.gain   # [ADU] -> [e-]

        # Placeholders for:
        #   - the smearing map after correction of the CTI (note that the first couple of rows will not be corrected (see further))
        #   - the smearing pattern (one entry per column in the sub-field)

        ctiCorrectedSmearingMap = np.zeros(smearingMap.shape)
        smearingPatternFastCadence = np.zeros(smearingMap.shape[1])
    
        # Smearing map at pixel (i,j) after correction for the CTI:
        #   Ic(i,j) = I(i,j) - a0(i,j) * [exp(-b0 * i) + a1 * exp(-b1 * i) + a2 * exp(-b2 * i) + a3 * exp(-b3 * i)]
        #           = I(i,j) - a0(i,j) * tau(i)
        #           = I(i,j) - a0(i,j) * [exp(-b0)^i + a1 * exp(-b1)^i + a2 * exp(-b2)^i + a3 * exp(-b3)^i]
        #           = I(i,j) - a0(i,j) * [u0^i + (a1 * u1^i) + (a2 * u2^i) + (a3 * u3^i)]
        #           = I(i,j) - a0(i,j) * [(c0 * f0(i)) + (c1 * f1(i)) +  (c2 * f2(i)) + (c3 * f3(i))]
        #   with
        #       uk = exp(-bk)
        #       fk(i) = uk^i        

        a1, a2, a3 = self.smearingCoefficientsA[0], self.smearingCoefficientsA[1], self.smearingCoefficientsA[2]                            # a0 is not in the array (will be updated)
        b0, b1, b2, b3 = self.smearingCoefficientsB[0], self.smearingCoefficientsB[1], self.smearingCoefficientsB[2], self.smearingCoefficientsB[3]
        u0, u1, u2, u3 = math.exp(-b0), math.exp(-b1), math.exp(-b2), math.exp(-b3)                     # Eq. (12) in PLATO-LESIA-PDC-TN-



        for column in range(smearingMap.shape[1]):
        
            # CTI correction

            f0, f1, f2, f3 = 1, 1, 1, 1     # Will be updated iteratively -> initialisation for i = 0
            tau = (1 + a1 + a2 + a3)        # Will be updated iteratively -> initialisation for i = 0 -> Eq. (3) in PLATO-LESIA-PDC-TN-

            for row in range(numRowsSmearingMap):

                if row >= self.smearingNumRowsSkipped:

                    # Correct the value

                    ctiCorrectedSmearingMap[row][column] = smearingMap[row][column] - self.smearingCoefficientA0Array[column] * tau
                
                # Update f0, f1, f2, f3, and tau
            
                f0 *= u0
                f1 *= u1
                f2 *= u2
                f3 *= u3

                tau = f0 + a1 * f1 + a2 * f2 + a3 * f3



            # Calculation of the mean smearing (first n0 measurements are excluded) 
        
            if self.includeSmearingOutlierDetection:
                
                smearingFlag = self.detectSmearingOutliers(ctiCorrectedSmearingMap[self.smearingNumRowsSkipped:, column], stdDevPrevious[column])

                if(np.sum(smearingFlag) == np.size(smearingFlag)):

                    smearingPatternFastCadence[column] = 0
                    continue
                
                else: 

                    smearingPatternFastCadence[column], stdDevPrevious[column] = applyShiftedDataAlgorithm(ctiCorrectedSmearingMap[self.smearingNumRowsSkipped:, column][~smearingFlag])[:2]

            else:

                smearingPatternFastCadence[column], stdDevPrevious[column] = applyShiftedDataAlgorithm(ctiCorrectedSmearingMap[self.smearingNumRowsSkipped:, column])[:2]



            # Update a0 for the current column and exposure n:
            #       a0(j; n) = [chi + epsilon * rho * ao(j; n - 1)] / [rho * (1 + epsilon)]
            # with
            #       chi = sum_i [I(i,j) - S(j)] * tau(i)
            #       rho = sum_i [tau(i)]^2

            chi, rho = 0, 0
            f0, f1, f2, f3 = 1, 1, 1, 1 # Will be updated iteratively -> initialisation for i = 0
            tau = (1 + a1 + a2 + a3)    # Will be updated iteratively -> initialisation for i = 0 -> Eq. (3) in PLATO-LESIA-PDC-TN-

            for row in range(numRowsSmearingMap):
            
                if (row >= self.smearingNumRowsSkipped) and (not self.includeSmearingOutlierDetection or (smearingFlag[row - self.smearingNumRowsSkipped] == 0)):

                    chi += (smearingMap[row][column] - smearingPatternFastCadence[column]) * tau
                    rho += pow(tau, 2)

                # Update f0, f1, f2, f3, and tau
            
                f0 *= u0
                f1 *= u1
                f2 *= u2
                f3 *= u3

                tau = f0 + a1 * f1 + a2 * f2 + a3 * f3

            self.smearingCoefficientA0Array[column] = (chi + self.smearingRegularization * rho * self.smearingCoefficientA0Array[column]) / (rho * (1 + self.smearingRegularization)) 
    
        return smearingPatternFastCadence, stdDevPrevious





    def detectSmearingOutliers(self, ctiCorrectedSmearingColumn, stdDev):

        """
        PURPOSE: Outlier detection in the parallel over-scan as explained in PLATO-MPSSR-PDC-PT-0003
                 (PLATO Onboards Overscan Outlier Detection Algorithm Theoretical Baseline Document).
                 Algorithm name: ONB-OVEROUTDET-010.
    
        INPUT:
            - ctiCorrectedSmearingColumn: Column from the parallel over-scan, after CTI correction (and discarding
                                          rows to avoid contamination by bright sources at the top of the detector)
            - stdDev: Standard deviation of the previous measurement of this column of the parallel over-scan
    
        OUTPUT:
            - smearingFlag: Boolean truth values, where 0 means "no outlier" and 1 means "outlier".  The 
                            truth value 1 is equivalent to a flag.
        """

        median = np.median(ctiCorrectedSmearingColumn)

        # Sigma-clipping around the median (we use the standard deviation of the previous measurement
        # of this column as sigma)

        smearingFlag = (ctiCorrectedSmearingColumn - median >= self.smearingOutlierDetectionThreshold * stdDev)

        return smearingFlag




    ##################
    # Image extraction
    ##################

    def getStarPhotons(self, exposure):

        """
        PURPOSE: Extract the contribution of the photons coming from the stars
                 in the simulated sub-field.
    
        INPUT:
            - exposure: Index of the exposure for which to return the background subtracted
                        flux, expressed in photons

        OUTPUT:
            - image: Fluxes in the simulated sub-field for the given exposure, after correcting 
                     for the offset and the smearing pattern, converting from ADU to photons, 
                     and subtracting the sky background.
            - throughput: Throughput map for the given exposure
        """

        # Subtract the offset and the smearing pattern, and account for the gain

        image = self.simFile.getImage(exposure)                                                 # [ADU]
        image = image - self. offsetValueArrayFastCadence[-1]                                   # Subtract the offset [ADU]
        image = image * self.gain                                                               # Multiply with gain [ADU] -> [e-]
        image = (image.transpose() - self.smearingPatternArrayFastCadence[-1]).transpose()      # Subtract smearing pattern (per column) -> was already multiplied with gain [e-]



        # Divide the exposure by the throuhgput map and the flatfield map

        throughput = np.array(self.simFile.getThroughputMap(exposure))
        
        image = np.divide(image, throughput, out = np.zeros_like(image), where = (throughput != 0))      # [Photons]
        image /= self.flatfield



        # Background subtraction

        image -= self.skyBackgroundArrayFastCadence[-1]



        return image, throughput





    ###################################
    # Calculation of the (nominal) mask
    ####################################                                     

    def calculateMask(self, exposure, image, throughput):

        """
        PURPOSE: Calculation of the (nominal) mask.

        INPUT:
            - exposure: Index of the exposure for which to update the (nominal) mask
            - image: Fluxes in the simulated sub-field for the given exposure, after correcting 
                     for the offset and the smearing pattern, converting from ADU to photons, 
                     and subtracting the sky background.
            - throughput: Throughput map for the given exposure

        OUTPUT:
            - masks: Nominal masks for all target stars
        """

        # Foresee a group in the output HDF5 file for the masks (for all targets) for the current exposure
        # (this will only done for the exposures for which the masks are updated!)

        maskGroupName = "/Masks/Mask{0:06d}".format(exposure)
        maskGroup = self.outputFile.create_group(maskGroupName)
        maskGroup.create_dataset("TargetIDs", data = self.targetIds)



        # Make placeholders for the information that will be stored (1st index = target index)

        masks = np.empty((self.numTargets, self.windowDimensions, self.windowDimensions), dtype = np.int)
        maskRowsInSubField = np.zeros((self.numTargets, int(math.pow(self.windowDimensions, 2))), dtype = np.int)        # Store the row coordinates of the pixels in the masks for each target in the output file (in the group for the current masks)
        maskColumnsInSubField = np.zeros((self.numTargets, int(math.pow(self.windowDimensions, 2))), dtype = np.int)     # Store the column coordinates of the pixels in the masks for each target in the output file (in the group for the current masks)
        maskSizes = np.zeros(self.numTargets)                                                                            # Store the size of the masks



        # Information for *all* sources (not only the target stars):
        #   - IDs
        #   - row coordinates [pixels]
        #   - column coordinates [pixels]
        #   - focal-plane x-coordinates [mm] (not needed)
        #   - focal-plane y-coordinates [mm] (not needed)
        #   - flux as derived from the catalogue magnitude [photons / exposure]

        starIds, starRows, starColumns, inputFlux = (self.simFile.getStarCoordinates(exposure)[index] for index in [0, 1, 2, 5])
    



        # Loop over all target stars and update their mask

        for targetIndex in range(self.numTargets):                                  # Index of the current target star in the list of target stars

            starIndex = np.where(starIds == self.targetIds[targetIndex])[0]         # Index of the current target star in the list with all stars

            if(len(starIndex) == 0):

                raise Exception("No star with ID", self.targetIds[targetIndex], "found")

            targetRow = starRows[starIndex]                                         # Row coordinate of the current target star in the simulated sub-field [pixels]
            targetColumn = starColumns[starIndex]                                   # Column coordinate of the current target star in the simulated sub-field [pixels]
            targetRowAsInt, targetColumnAsInt = int(targetRow), int(targetColumn)   # Round the target coordinates to the centre of the pixels

            self.targetRowsArray[targetIndex] = targetRow
            self.targetColumnsArray[targetIndex] = targetColumn

            if (targetRowAsInt <= self.windowOffset) or (targetRowAsInt + self.windowOffset >= self.numRowsSubField) or (targetColumnAsInt <= self.windowOffset) or (targetColumnAsInt + self.windowOffset >= self.numColumnsSubField):

                raise Exception("Mask around target", self.targetIds[targetIndex], "falls (at least partially) off the sub-field")


            # We want to use a square mask, centred around the current target position.  We only want
            # to include the pixels in the mask, such that they contribute (together) more to the 
            # signal of the current target star than to the noise.  The noise consists of the contaminant 
            # stars (within the specified distance), readout noise, and sky background.

            # Create a theoretical sub-field image as if only the current target star was on the CCD

            targetInputFlux = inputFlux[starIndex] * throughput[targetRowAsInt, targetColumnAsInt] * self.flatfield[targetRowAsInt, targetColumnAsInt]      # Flux of the current target star [e- / exposure]
            targetMap = self.analyticPSF.getPSF(targetRow, targetColumn, targetInputFlux, self.subfieldZeropointRow, self.subfieldZeropointColumn, self.numRowsSubField, self.numColumnsSubField)

            # Create a theoretical sub-field as if only the contaminants of the current target star were on the CCD
            # (a star cannot not be its own contaminant!)
            
            contaminantIndices = np.where((np.abs(starRows - targetRow) <= self.contaminationRadius) & (np.abs(starColumns - targetColumn) <= self.contaminationRadius) & (starIds != self.targetIds[targetIndex]))[0]  # Other stars within the contaminant radius

            if self.hasKnownContaminants:

                knowContaminantIndices = np.where(starIds == self.contaminantIds)[0]
                contaminantIndices = np.intersect1d(contaminantIndices, knowContaminantIndices)     # Known contaminants within the contaminant radius

            # numContaminants = len(contaminantIndices)

            contaminantsMap = np.zeros_like(image)

            for contaminantIndex in contaminantIndices:

                inputFluxContaminant = inputFlux[contaminantIndex] * throughput[int(starRows[contaminantIndex]), int(starColumns[contaminantIndex])]    # Flux of the current contaminant [e- / exposure]
                contaminantsMap += self.analyticPSF.getPSF(starRows[contaminantIndex], starColumns[contaminantIndex], inputFluxContaminant, self.subfieldZeropointRow, self.subfieldZeropointColumn, self.numRowsSubField, self.numColumnsSubField)  # Add the current contaminant to the contaminants map
            


            # Given the contaminants map for the current target, decide which pixels 
            # should be included in the mask for the current target

            minRow = max(0, targetRowAsInt - self.windowOffset)                                 # Minimum row index of the mask
            maxRow = min(self.numRowsSubField, targetRowAsInt + self.windowOffset)              # Maximum row index in the mask (incl.)
            minColumn = max(0, targetColumnAsInt - self.windowOffset)                           # Minimum column index in the mask
            maxColumn = min(self.numColumnsSubField, targetColumnAsInt + self.windowOffset)     # Maximum column index in the mask (incl.)

            varianceMap = targetMap + contaminantsMap + self.varianceReadoutNoise + self.skyBackgroundArrayFastCadence[-1] * throughput * self.flatfield
            noiseToSignalRatio = (np.sqrt(varianceMap) / targetMap)[minRow : maxRow + 1, minColumn : maxColumn + 1]

            # Sort the pixels in the NSR map in ascending order

            maskRows, maskColumns = np.unravel_index(np.argsort(noiseToSignalRatio.ravel()), noiseToSignalRatio.shape)  # Coordinates in the mask

            # Initialize with the first pixel (with the lowest NSR)

            pixelRow, pixelColumn                 = maskRows[0], maskColumns[0]
            maskRowsInSubField[targetIndex, 0]    = pixelRow
            maskColumnsInSubField[targetIndex, 0] = pixelColumn
            aggregatedVariance                    = varianceMap[pixelRow, pixelColumn]
            aggregatedTheoreticalTargetFlux       = targetMap[pixelRow, pixelColumn]
            aggregatedObservedTargetFlux          = image[pixelRow, pixelColumn]
            aggregatedNSR                         = noiseToSignalRatio[pixelRow, pixelColumn]
            maskSize = 1

 
            # Then add other pixels

            for maskPixelIndex in range(1, len(maskRows)):
                
                pixelRow, pixelColumn = maskRows[maskPixelIndex], maskColumns[maskPixelIndex]
                
                temp = np.sqrt(aggregatedVariance + varianceMap[pixelRow, pixelColumn]) / (aggregatedTheoreticalTargetFlux + targetMap[pixelRow, pixelColumn])
                
                if temp < aggregatedNSR:

                    maskRowsInSubField[targetIndex, maskPixelIndex] = pixelRow
                    maskColumnsInSubField[targetIndex, maskPixelIndex] = pixelColumn
                    maskRowsInSubField[targetIndex, maskPixelIndex] = pixelRow
                    maskColumnsInSubField[targetIndex, maskPixelIndex] = pixelColumn
                    aggregatedVariance              += varianceMap[pixelRow, pixelColumn]
                    aggregatedTheoreticalTargetFlux += targetMap[pixelRow, pixelColumn]
                    aggregatedObservedTargetFlux    += image[pixelRow, pixelColumn]
                    aggregatedNSR = temp
                    maskSize += 1

                else:

                    # Copy the exact pixels of the mask, so that we can persist them in the HDF5 file

                    targetMask = np.zeros((self.windowDimensions, self.windowDimensions))

                    for index in range(maskSize):
                        targetMask[maskRows[index]][maskColumns[index]] = 1                  

                    # Disregard all other pixels of the window around the target star:
                    # they all contribute more to the noise than to the signal.

                    break

            maskRowsInSubField[targetIndex] += minRow
            maskColumnsInSubField[targetIndex] += minColumn

            masks[targetIndex] = targetMask

        maskGroup.create_dataset("MaskPixelRows", data = maskRowsInSubField)
        maskGroup.create_dataset("MaskPixelColumns", data = maskColumnsInSubField)
        maskGroup.create_dataset("MaskSize", data = maskSizes)

            
        
        return masks






    #############################################
    # Flux & COB calculation using (nominal) mask
    #############################################

    def calculateFluxAndCob(self, exposure, offsetValueFastCadence, smearingPatternFastCadence, targetIndex, mask):

        """
        PURPOSE: Flux and COB calculation as explained in PLATO-LESIA-PDC-DD-008 
                 (PLATO Onboards Flux & COB Calculation ATBD).
                 Algorithm name: ONB-FXCOBCAL-010.

        INPUT:
            - exposure: Index of the exposure to process
            - offsetValueFastCadence: Electronic offset as calculated from the serial pre-scan
            - windowSmearingPatternFastCadence: Smearing pattern for the window as calculated from the parallel over-scan
            - targetIndex: Index of the target source to process
            - mask: Mask (nominal)

        OUTPUT:
            - fluxFastCadence: Flux computed using the (nominal) mask
            - cobFastCadence: COB computed using the (nominal) mask
        """

        # Extract the star window and the corresponding smearing pattern

        windowZeropointRow = int(self.targetRowsArray[targetIndex]) - self.windowOffset
        windowZeropointColumn = int(self.targetColumnsArray[targetIndex]) - self.windowOffset

        starWindow = self.simFile.getImage(exposure)[windowZeropointRow - self.subfieldZeropointRow : windowZeropointRow + self.windowDimensions - self.subfieldZeropointRow, windowZeropointColumn - self.subfieldZeropointColumn : windowZeropointColumn + self.windowDimensions - self.subfieldZeropointColumn]
        windowSmearingPatternFastCadence = smearingPatternFastCadence[windowZeropointColumn - self.subfieldZeropointColumn : windowZeropointColumn + self.windowDimensions - self.subfieldZeropointColumn]

        maskedWindow = starWindow - offsetValueFastCadence                                            # Subtract the offset
        maskedWindow = maskedWindow * self.gain                                                       # Multiply with gain [ADU] -> [e]

        maskedWindow = (maskedWindow.transpose() - windowSmearingPatternFastCadence).transpose()      # Subtract smearing pattern (per column)
        maskedWindow = np.multiply(maskedWindow, mask)                                                # Multiply with mask (element-wise)

        fluxFastCadence = np.sum(maskedWindow)

        cobFastCadence = [np.sum(np.multiply(mask, self.cobOffsetX)) / fluxFastCadence, np.sum(np.multiply(mask, self.cobOffsetY)) / fluxFastCadence]
        
        return fluxFastCadence, cobFastCadence





    # #########################
    # # Smearing time averaging
    # #########################

    # def timeAverageSmearing(self):

    #     """
    #     PURPOSE: Smearing time averaging for long cadence as explained in PLATO-LESIA-PDC-DD-007
    #              (PLATO Ob-board Smearing Pattern Time Averaging ATBD).
    #              Algorithm name: ONB-SMRAVG-010
    
    #     OUTPUT:
    #         - smearing _patternLongCadence: Smearing pattern averaged over the last 24 samples (i.e. 600s) [e-]
    #     """

    #     # Average out the last 24 samples per column

    #     smearingPatternArrayLastLongCadence = self.smearingPatternArrayFastCadence[-int(self.numExposuresLongCadence) :]     # Last 24 samples

    #     shape = (1, self.numExposuresLongCadence, smearingPatternArrayLastLongCadence.shape[1], 1)
    #     smearingPatternLongCadence = smearingPatternArrayLastLongCadence.reshape(shape).mean(-1).mean(1)

    #     return smearingPatternLongCadence




    
    # ####################################
    # # Outlier detection over light curve
    # ####################################

    # def detectFluxOutliers(self, fluxArrayFastCadence, fluxFlagShortCadence):

    #     """
    #     PURPOSE: Outlier detection over light curve as explained in PLATO-MPSSR-PDC-DD-0001
    #              (PLATO OnboardLongCadence Outlier Detection ATBD).
    #              Algorithm name: ONBLongCadenceOUTDET-010.

    #     INPUT:
    #         - fluxArrayFastCadence: Flux time series processed so far. This can either be the flux 
    #                                 obtained in the nominal mask (fx_fc) or the difference in flux 
    #                                 between the extended and the nominal mask (dfx_fc).
    #         - fluxFlagShortCadence: Flag for the time series processed so far, apart from
    #                                    the last 2 * b + 1 samples
        
    #     OUTPUT:
    #         - fluxFlagShortCadence: Flag for the time series processed so far, including
    #                                 the last 2 * b + 1
    #     """

    #     fluxArrayLastBin = fluxArrayFastCadence[-2 * self.fluxOutlierDetectionHalfBinWidth - 1:]

    #     # Step 1 - 4

    #     isOutlier = applyMadClippingAroundMedian(fluxArrayLastBin, self.lightCurveOutlierDetectionThreshold, self.fluxOutlierDetectionHalfBinWidth)

    #     # Step 5

    #     if fluxFlagShortCadence[-2] and fluxFlagShortCadence[-1]:

    #         fluxFlagShortCadence[-2] = False

    #         if not isOutlier:

    #             fluxFlagShortCadence[-1] = False

    #     fluxFlagShortCadence = np.append(fluxFlagShortCadence, isOutlier)

    #     return fluxFlagShortCadence





    # ##########################
    # # Flux & COB time averaging
    # ###########################

    # def timeAverageFluxAndCob(self):

    #     """
    #     PURPOSE: Flux time averaging for short cadence (50s) as explained in PLATO-LESIA-PDC-DD-009
    #              (PLATO Onboard Flux and COB Short Cadence Time Averaging ATBD) or long cadence as 
    #              explained in PLATO-LESIA-PDC-DD-010 (PLATO Onboard Flux and COB Long Cadence Time Averaging ATBD).
    #              Algorithm names: ONB-FXAGV-011 and ONB-FXAGV-012.
    #     """

    #     # Short cadence

    #     if(self.timeAveragingCadence == "Short"):
            
    #         fluxValueShortCadence, cobShortCadence, numUsefulDatapointsShortCadence, numFlaggedDatapointsShortCadence = self.timeAverageFluxAndCobShortCadence(self.fluxArrayFastCadence, self.cobArrayFastCadence, self.fluxFlagFastCadence)           # Select duration of short cadence (50s)

    #         if numUsefulDatapointsShortCadence != 0:

    #             self.fluxArrayShortCadence = np.append(self.fluxArrayShortCadence, fluxValueShortCadence)
    #             self.cobArrayShortCadence = np.append(self.cobArrayShortCadence, cobShortCadence)
    #             self.fluxOutlierFlagShortCadence = np.append(self.fluxOutlierFlagShortCadence, numFlaggedDatapointsShortCadence)

        

    #     # Long cadence

    #     elif self.timeAveragingCadence == "Long":

    #         fluxValueLongCadence, fluxVarianceLongCadence, numUsefulDatapointsLongCadence, numFlaggedDatapointsLongCadence = self.timeAverageFluxAndCobLongCadence(self.fluxArrayFastCadence, self.fluxFlagFastCadence)         # Select duration of long cadence (600s)

    #         if numUsefulDatapointsLongCadence != 0:
            
    #             self.fluxArrayLongCadence = np.append(self.fluxArrayLongCadence, fluxValueLongCadence)
    #             self.fluxVarianceArrayLongCadence = np.append(self.fluxVarianceArrayLongCadence, fluxVarianceLongCadence)
    #             self.fluxOutlierFlagLongCadence = np.append(self.fluxOutlierFlagLongCadence, numFlaggedDatapointsLongCadence)





    # def timeAverageFluxAndCobShortCadence(self, fluxArrayFastCadence, cobArrayFastCadence, fluxFlagShortCadence = None):

    #     """
    #     PURPOSE: Flux time averaging for short cadence (50s) as explained in PLATO-LESIA-PDC-DD-009
    #              (PLATO Onboard Flux and COB Short Cadence Time Averaging ATBD).
    #              Algorithm name: ONB-FXAGV-011.

    #     INPUT:
    #         - fluxArrayFastCadence: Flux time series processed so far. This can either be the flux 
    #                                 obtained in the nominal mask (fx_fc) or the difference in flux 
    #                                 between the extended and the nominal mask (dfx_fc).
    #         - cobArrayFastCadence: COB time series processed so far.  This can either be the COB
    #                                obtained from the nominal mask or from the extended mask.
    #         - fluxFlagShortCadence: Outlier detection flags for the flux time series processed so far.  If None, 
    #                                outlier detection on the light curve was disabled.

    #     OUTPUT:
    #         - fluxValueShortCadence: Mean of the given flux time series over the last 50s, Calculated with the non-flagged datapoints only
    #         - cobShortCadence: Mean of the given COB time series over the last 50s, Calculated with the non-flagged datapoints only
    #         - numUsefulDatapointsShortCadence: Number of non-flagged datapoints over the last 50s
    #         - numFlaggedDatapointsShortCadence: Number of flagged datapoints over the last 50s
    #     """

    #     # Select the last 50s (i.e. duration of short cadence) in the flux and COB time series

    #     fluxArrayLastShortCadence = fluxArrayFastCadence[-int(self.numExposuresShortCadence) :]
    #     cobArrayLastShortCadence = cobArrayFastCadence[-int(self.numExposuresShortCadence) :]

    #     # Outlier detection enabled

    #     if self.includeLightCurveOutlierDetection:

    #         # The flag that has been derived during the outlier detection on the flux values will also
    #         # be used as flag for the COB here

    #         flagLastShortCadence = fluxFlagShortCadence[-int(self.numExposuresShortCadence) :]                       # Flag for the last 50s
    #         numFlaggedDatapointsShortCadence = np.sum(flagLastShortCadence)                                          # Number of flagged datapoints over the last 50s
    #         numUsefulDatapointsShortCadence = self.numExposuresShortCadence - numFlaggedDatapointsShortCadence       # Number of non-flagged datapoints over the last 50s

    #         fluxValueShortCadence = np.sum(fluxArrayLastShortCadence[~flagLastShortCadence]) / numUsefulDatapointsShortCadence      # Mean of the non-flagged flux values over the last 50s 

    #         cobShortCadenceX = np.sum(cobArrayLastShortCadence[~flagLastShortCadence][0]) / numUsefulDatapointsShortCadence         # Mean of the non-flagged x-coordinates of the COB over the last 50s
    #         cobShortCadenceY = np.sum(cobArrayLastShortCadence[~flagLastShortCadence][1]) / numUsefulDatapointsShortCadence         # Mean of the non-flagged y-coordnates of the COB over the last 50s
    #         cobShortCadence = [cobShortCadenceX, cobShortCadenceY]                                                                  # Couple the x- and y-coordinates

    #     # Outlier detection disabled

    #     else:

    #         # No flagging

    #         numFlaggedDatapointsShortCadence = 0
    #         numUsefulDatapointsShortCadence = self.numExposuresShortCadence

    #         fluxValueShortCadence = np.sum(fluxArrayLastShortCadence) / numUsefulDatapointsShortCadence         # Mean of the flux values over the last 50s

    #         cobShortCadenceX = np.sum(cobArrayLastShortCadence[0]) / numUsefulDatapointsShortCadence            # Mean over the x-coordinates of the COB over the last 50s
    #         cobShortCadenceY = np.sum(cobArrayLastShortCadence[1]) / numUsefulDatapointsShortCadence            # Mean over the y-coordinates of the COB over the last 50s
    #         cobShortCadence = [cobShortCadenceX, cobShortCadenceY]                                              # Couple the x- and y-coordinates

    #     return fluxValueShortCadence, cobShortCadence, numUsefulDatapointsShortCadence, numFlaggedDatapointsShortCadence





    # def timeAverageFluxAndCobLongCadence(self, fluxArrayFastCadence, fluxFlagShortCadence = None):

    #     """
    #     PURPOSE: Flux time averaging for long cadence (600s) as explained in PLATO-LESIA-PDC-DD-010
    #              (PLATO Onboard Flux and COB Long Cadence Time Averaging ATBD).
    #              Algorithm name: ONB-FXAGV-012.

    #     INPUT:
    #         - fluxArrayFastCadence: Flux time series processed so far. This can either be the flux 
    #                                 obtained in the nominal mask (fx_fc) or the difference in flux 
    #                                 between the extended and the nominal mask (dfx_fc).
    #         - fluxFlagShortCadence: Outlier detection flags for the flux time series processed so far.  If None, 
    #                                 outlier detection on the light curve was disabled.
    #     OUTPUT:
    #         - fluxValueLongCadence: Mean of the given flux time series over the last 600s, Calculated with the non-flagged datapoints only
    #         - fluxVarianceLongCadence: Variance of the given flux time series over the last 600s, Calculated with the non-flagged datapoints only
    #         - numUsefulDatapointsLongCadence: Number of non-flagged datapoints over the last 600s
    #         - numFlaggedDatapointsLongCadence: Number of flagged datapoints over the last 600s
    #     """

    #     # Select the last 600s (i.e. duration of long cadence) in the flux time series

    #     fluxArrayLastLongCadence = fluxArrayFastCadence[-int(self.numExposuresLongCadence) :]

    #     # Outlier detection enabled

    #     if self.includeLightCurveOutlierDetection:

    #         flagLastLongCadence = fluxFlagShortCadence[-int(self.numExposuresLongCadence) :]        # Flag for the last 600s
            
    #         fluxValueLongCadence, fluxVarianceLongCadence, numUsefulDatapointsLongCadence = applyShiftedDataAlgorithm(fluxArrayLastLongCadence, flagLastLongCadence)
    #         numFlaggedDatapointsLongCadence = np.sum(flagLastLongCadence)

    #     # Outlier detection disabled

    #     else:
            
    #         fluxValueLongCadence, fluxVarianceLongCadence, numUsefulDatapointsLongCadence = applyShiftedDataAlgorithm(fluxArrayLastLongCadence)
    #         numFlaggedDatapointsLongCadence = 0

    #     return fluxValueLongCadence, fluxVarianceLongCadence, numUsefulDatapointsLongCadence, numFlaggedDatapointsLongCadence










def getPhotometryTimeSeries(filename, targetId):

    """
     PURPOSE: Extract the flux time series of star with a given star ID.

     INPUT: 
        - filename: Name of the HDF5 output file written by the photometric pipeline
        - targetId:  Target identifier (integer, e.g. 9789)

     OUTPUT: 
        - time: Time points [s]
        - inputFlux: Input flux of the target star as derived from the input catalogue
        - outputFlux: Flux of the target star as calculated by the photometric pipeline

     REMARK: To find out which star identifiers are in the photometry file, look in the HDF5 simulation
             output file of PlatoSim: 
             allStarIDs = np.array(platosimOutputFile["StarCatalog/starIDs"])
    """

    photFile = h5py.File(filename)
    time =  np.array(photFile["/Photometry/time"])
    numExposures = len(time)

    targetIds = photFile[""]

    inputFlux = []
    outputFlux = []

    for exposure in range(numExposures):

        allStarIDsInImage = np.array(photFile["/Photometry/Exposure{0:06d}/starID".format(exposure)])
        
        if targetId in allStarIDsInImage:

            inputFluxAllTargets = np.array(photFile["Photometry/Exposure{0:06d}/InputFlux".format(exposure)])
            outputFluxAllTargets = np.array(photFile["Photometry/Exposure{0:06d}/FluxFastCadence".format(exposure)])

            inputFlux.append(inputFluxAllTargets[np.where(targetIds == targetId)][0])
            outputFlux.append(outputFluxAllTargets[np.where(targetIds == targetId)][0])

    inputFlux = np.array(inputFlux)
    outputFlux = np.array(outputFlux)

    photFile.close()

    return time, inputFlux, outputFlux
 