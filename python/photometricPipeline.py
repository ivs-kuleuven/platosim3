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



###############################
# Background window calculation
###############################

# def processBackgroundWindow(background_window, offset_value_fc, smearing_pattern_fc, half_ccd_gain, outlier_detection_threshold):

#     """
#     PURPOSE: Background calculation as explained in PLATO-LESIA-PDC-DD-011 (PLATO Onboard Background
#              Window Calculation ATBD).
#              Algorithm name: ONB-BKGCAL-010.

#     INPUT:
#         - background_window: values in the background window [ADU]
#         - offset_value:_fc Electronic offset as calculated from the serial pre-scan
#         - smearing_pattern_fc: Smearing as calculated from the parallel over-scan
#         - half_ccd_gain: CCD gain for the detector half to which the given parallel over-scan
#                          corresponds
#         - outlier_detection_threshold: Threshold for outlier detection
    
#     OUTPUT:
#         - background_value_fc: Mean of the values in the background window [e-]
#         - background_variance_fc: Variance of the values in the background window [e-]
#     """

#     window_smearing_pattern_fc = window_smearing_pattern(smearing_pattern_fc)

#     # Subtract offset, multiply with gain, and subtract smearing pattern

#     background_window = background_window - offset_value_fc

#     background_window *= half_ccd_gain

#     for column in range(background_window.shape[1]):

#         background_window[:][column] = background_window[:][column] - window_smearing_pattern_fc[column]
    


#     # Outlier detection

#     flag = background_outlier_detection(background_window, outlier_detection_threshold)[0]

#     # Shifted data algorithm

#     background_value_fc, background_variance_fc, n_useful = shifted_data_algorithm(background_window[~flag])

#     return background_value_fc, background_variance_fc





# # def background_outlier_detection(background_window, threshold):

#     """
#     PURPOSE: Detection of outliers in the background window as explained in PLATO-MPSSR-PDC-DD-004
#              (PLATO On-board Background Window Outlier Detection Algorithm Theoretical Baseline
#              Document).
#              Algorithm name: ONB-BACKOUTDET-010.

#     INPUT:
#         - background_window: Array of background window pixels [-e]
#         - threshold: Threshold for outlier detection (number of median absolute deviations used to flag
#                      outliers)

#     OUTPUT:
#         - raw_back_flags: Boolean truth values, where 0 means "no outlier" and 1 means "outlier".  The
#                           truth value 1 is equivalent to a flag
#     """

#     raw_back_flags = mad_median_clipping(background_window, threshold)

#     return raw_back_flags





#####################
# Auxiliary functions
#####################

def shifted_data_algorithm(data):

    """
    PURPOSE: Shifted data algorithm to compute the mean and variance for the given data.  This method
             avoids loss of significance with big numbers when computing the variance.

    INPUT:
        - data: Data array for which to calculated the mean and the variance

    OUTPUT:
        - mean: Mean value of the given data array, as computed with the shifted data algorithm
        - variance: Variance of the given data array, as computed with the shifted data algorithm
        - n_useful: Number of non-flagged datapoints
    """

    k = data[0]         # Use the 1st non-flagged value as approximation of the mean (arbitrary choice)
    data = data - k     # Shift the data

    data_shifted_sum = np.sum(data)                             # Sum of the shifted data
    data_shifted_squared_sum = np.sum(np.power(data, 2))        # Sum of the squares of the shifted data

    n_useful = len(data)

    mean = (data_shifted_sum + n_useful * k) / n_useful                                                 # Mean
    variance = (data_shifted_squared_sum - (pow(data_shifted_sum, 2) / n_useful)) / (n_useful - 1)      # Variance

    return (mean, variance, n_useful)





def mad_median_clipping(data, threshold, index = None):

    """
    PURPOSE: Outlier detection, based on Median Absolute Deviation (MAD) clipping around the median.

    INPUT:
        - data: Data for which to flag outliers
        - threshold: Threshold for MAD clipping around the median
        - index: If specified, it is calculated whether the datapoint at this position is an outlier.

    OUTPUT:
        - in case index was specified: boolean indicating whether or not the datapoint at the given
          position is an outlier; otherwise: boolean truth values, where 0 means "no outlier" and 1 means
                                 "outlier".  The truth value 1 is equivalent to a flag.
    """

    median = np.median(data)         # Median
    y = np.fabs(data - median)       # Subtract median from background window
    mad = np.median(y)               # MAD

    if index == None:

        flag = np.where(y > threshold * mad)
        return flag

    else:
        return y[index] / mad > threshold





def window_smearing_pattern(smearing_pattern_fc):

    # TODO We will need extra parameters (spanned columns) to extract the required smearing
    # pattern for the window.

    return smearing_pattern_fc










class PhotometricPipeline(object):
    
    def __init__(self, runName, configurationFile, simulationFile, outputDir = None, debug = False):

        """
        PURPOSE: Initialisation of the class variables and readong of the default input files.

        INPUT:
            - runName: Name of the run of the photometric pipeline (will be used as name for the output file)
            - configurationFile: Name of the configuration file to run the photometric pipeline
            - simulationFile: Name of the output file of PlatoSim (this file contains the serial pre-scan,
                              parallel over-scan, and fluxes in a simulated sub-field on the CCD)
            - outputDir: Directory in which the results of the photometric pipeline will be stored
            - debug: Indicates whether or not debugging messages should be printed
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

        # Debugging

        self.debug = debug
        




    @property
    def outputDir(self):

        """
        PURPOSE: Return location of the output file.

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
            - Path to the directory to use as output directory
        """

        if not os.path.exists(path):
            if self.debug:

                print("DEBUG: creating output directory {}.".format(path))

            self.createDirectory(path)
        
        self.targetOutputFilesLocation = path
        self.hasTargetLocation = True

        if self.debug:

            print("DEBUG: output dir set to {}.".format(path))





    def readConfigurationFile(self, filename):

        """
        PURPOSE: Read the YAML input configuration file. 

        INPUT:
            - Path to the configuration file
        """

        self.configurationFilename = filename

        if self.debug:

            print ("Parsing YAML configuration file {}.".format(filename))
        
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

            parentNodeName, nodeName = key, None
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
    
        if self.debug:
            print ("Writing the Yaml configuration file {}.".format(filename))
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
        #logFilename = "{}/{}.log".format(self.targetOutputFilesLocation, self.runName)

        if removeOutputFile:
        
            try:
                os.remove(outputFilename)
            except OSError:
                pass
        
        self.writeYamlConfigurationFile(inputFilename)

        

        # Parameters that are specific for the star window

        self.window_zp_row = self["StarWindow/ZeropointRow"]            # Row coordinate of the star window zeropoint on the detector [pixels]
        self.window_zp_column = self["StarWindow/ZeropointColumn"]      # Column coordinate of the star window zeropoint on the detector [pixels]
        self.window_num_rows = self["StarWindow/NumRows"]               # Height of the star window [pixels]
        self.window_num_columns = self["StarWindow/NumColumns"]         # Width of the star window [pixels]



        # Parameters that are specific for the detector half on which the star window is located

        self.subfield_zp_row = self.simFile.getInputParameter("SubField", "ZeroPointRow")           # Row coordinate of the sub-field zeropoint on the detector [pixels]
        self.subfield_zp_column = self.simFile.getInputParameter("SubField", "ZeroPointColumn")     # Column coordinate of the sub-field zeropoint on the detector [pixels]

        self.detector_half = "Left"
        if self.window_zp_column >= self.simFile.getInputParameter("CCD", "NumColumns") / 2:
            self.detector_half = "Right"

        if self.detector_half == "Left":
            self.half_ccd_gain = 1.0 / (self.simFile.getInputParameter("FEE/Gain", "RefValueLeft") * self.simFile.getInputParameter("CCD/Gain", "RefValueLeft"))    # Total gain for the left detector half [e- / ADU]
        elif self.detector_half == "Right":
            self.half_ccd_gain = 1.0 / (self.simFile.getInputParameter("FEE/Gain", "RefValueRight") * self.simFile.getInputParameter("CCD/Gain", "RefValueRight"))  # Total gain for the left detector half [e- / ADU]



        # Parameters that are specific for the offset calculation

        self.offset_outlier_detection_enabled = self["Offset/OutlierDetection/Enabled"]                 # Enable/disable outlier detection
        self.offset_outlier_detection_k = self["Offset/OutlierDetection/k"]                             # Number of largest and smallest values to flag as outliers



        # Parameters that are specific for the smearing pattern calculations
        
        self.smearing_a0_array = np.empty(self.simFile.getInputParameter("SubField", "NumColumns")) 
        self.smearing_a0_array.fill(np.array(self["Smearing/Coefficients/a"])[0])
        self.smearing_a_coefficients = np.array(self["Smearing/Coefficients/a"])[1:]                    # Coefficients [a1, a2, a3]
        self.smearing_b_coefficients = np.array(self["Smearing/Coefficients/b"])                        # Coefficients [b0, b1, b2, b3]
        self.smearing_n0 = self["Smearing/NumRowsSkipped"]                                              # Number of rows to skip for CTI correction (1st rows may be affected by bright stars at the top of the detector) 
        self.smearing_outlier_detection_enabled = self["Smearing/OutlierDetection/Enabled"]             # Enable/disable outlier detection
        self.smearing_outlier_detection_threshold = self["Smearing/OutlierDetection/Threshold"]         # Threshold for outlier detection
        self.smearing_epsilon = self["Smearing/Regularization"]                                         # Epsilon

        self.smearing_pattern_fc_array = np.array([])                                                   # Smearing pattern for the different columns at fast cadence (25s)
        self.smearing_pattern_lc_array = np.array([])                                                   # Smearing pattern for the difference columns at long cadence (600s)



        # Parameters that are specific for the flux and COB calculation

        self.fx_fc_array = np.array([])         # Flux calculated with the nominal mask at fast cadence (25s)
        self.dfx_fc_array = np.array([])        # Difference in flux between the extended and the nominal mask at fast cadence (25s)
        self.ncob_fc_array = np.array([])       # COB calculated with the nominal mask at fast cadence (25s)
        self.ecob_fc_array = np.array([])       # COB calculated with the extende mask at fast cadence (25s)



        # Parameters that are specific for light curve outlier detection

        self.lc_outlier_detection_enabled = self["LightCurve/OutlierDetection/Enabled"]         # Enable/disable outlier detection
        self.lc_outlier_detection_b = self["LightCurve/OutlierDetection/b"]                     # Number of past and future datapoints needed to decide whether or not a datapoint is an outlier
        self.lc_outlier_detection_threshold = self["LightCurve/OutlierDetection/Threshold"]     # Outlier detection threshold


        if self.lc_outlier_detection_enabled:

            self.nflag_fc = np.zeros(self.lc_outlier_detection_b)
            self.eflag_fc = np.zeros(self.lc_outlier_detection_b)
        
        else:
            self.nflag_fc = None
            self.eflag_fc = None



        # Parameters that are specific for time averaging

        self.time_averaging_cadence = self["LightCurve/TimeAveraging/Cadence"]                  # Choose between short and long cadence

        self.num_exposures_sc = self["LightCurve/TimeAveraging/NumSamples/Short"]               # Short cadence: 50s (i.e. 2 samples)
        self.num_exposures_lc = self["LightCurve/TimeAveraging/NumSamples/Long"]                # Long cadence: 600s (i.e. 24 samples)

        if self.time_averaging_cadence == "Short":
            
            self.num_exposures_time_averaging = self.num_exposures_sc

        elif self.time_averaging_cadence == "Long":
            
            self.num_exposures_time_averaging = self.num_exposures_sc



        # Process all exposures

        numExposures = self.simFile.getInputParameter("ObservingParameters", "NumExposures")

        for exposure in range(numExposures):

            self.process_exposure(exposure)

        

        # Write results to HDF5 file

        self.write_output(outputFilename)
    




    def process_exposure(self, exposure):

        # Calculate the offset for the current exposure (fast cadence)

        if self.detector_half == "Left":
            serialPreScan = self.simFile.getBiasMapLeft(exposure)

        elif self.detector_half == "Right":
            serialPreScan = self.simFile.getBiasMapRight(exposure)

        offset_value_fc = self.offset_calculation(serialPreScan)[0]



        # Calculate the smearing pattern for the current exposure (fast cadence)

        parallelOverScan = self.simFile.getSmearingMap(exposure)
        std_dev_previous = 9999

        smearing_pattern_fc, std_dev_previous = self.smearing_calculation(parallelOverScan, offset_value_fc, std_dev_previous)
        self.smearing_pattern_fc_array = np.append(self.smearing_pattern_fc_array, smearing_pattern_fc)



        # Calculate the background for the current exposure (fast cadence)
        # TODO
        # Real photometric pipeline: ~100 background windows (nominally 4x4) per CCD half,
        # for which we know the fluxes and the position on the CCD.  The latter is needed to 
        # subtract the correct smearing pattern from each background window.  For each of the
        # background windows the mean and variance are calculated, which are used in an
        # interpolation schema (based on radial basis functions) to determine the background
        # at the position of the target star.
        # It seems unfeasible to make (long-term) simulations for all these background windows, so
        # we will have to look for an alternative way to estimate the background at the target
        # location.



        # Calculate the flux & COB (nominal & extended mask) for the current exposure
        # -> add new datapoint to the light curve (fast cadence)
        # TODO
        # This will fill fx_fc (flux in the nominal window at fast cadence), 
        # dfx_fc (flux difference between the extended and the nominal mask), 
        # ncob_fc (COB as obtained in the nominal mask), and ecob_fc (COB as
        # obtained in the extended mask).

        self.flux_cob_calculation(exposure, offset_value_fc,smearing_pattern_fc, nmask, emask)  # TODO Where do we get the masks from?



        # Smearing time averaging (fast -> long cadence)

        if (exposure + 1) % self.num_exposures_lc == 0:

            smearing_pattern_lc = self.smearing_time_averaging()
            self.smearing_pattern_lc_array = np.append(self.smearing_pattern_lc_array, smearing_pattern_lc)



        # Light curve outlier detection
        # Note that you need b datapoints in the past and b datapoints in the future!

        if self.lc_outlier_detection_enabled:

            if (exposure + 1) > 2 * self.lc_outlier_detection_b:
            
                self.nflag_fc = self.flux_cob_outlier_detection(self.fx_fc_array, self.nflag_fc)
                self.eflag_fc = self.flux_cob_outlier_detection(self.dfx_fc_array, self.eflag_fc)       # TODO Is this what we want to do?
                # TODO What at the end of the time series (when there are no future datapoints)?



        # Flux & COB time averaging

        if (exposure + 1) % self.num_exposures_time_averaging == 0:

            self.flux_cob_time_averaging()





    def write_output(self, filename):

        """
        PURPOSE: Write results of the photometric pipeline to an HDF5 file with the given name.

        INPUT:
            - filename: Name of the output file
        """

        # OFFSET_VALUE_FC
        # OFFSET_VARIANCE_FC
        # SMEARING_PATTERN_FC
        # BKG_WINDOWS
        # TARGET_STARS  -> long cadence
        # REFERENCE_STARS

        # {background window ID}_VALUE_FC
        # {background window ID}_VARIANCE_FC
        # {background window ID}_ERROR_NUMBER_FC
        # {star window ID}-FX_FC
        # {star window ID}-DFX_FC
        # {star window ID}-NCOB_FC
        # {star window ID}-ECOB_FC
        # -FX_SC
        # -DFX_SC
        # -NCOB_SC
        # -ECOB_SC
        # -FX_EXPOSURE_ERROR_SC
        # -FX_LC
        # -FXVAR_LC
        # -FX_EXPOSURE_ERROR_LC

        return None

    ####################
    # Offset calculation
    ####################

    def offset_calculation(self, offset_rows):

        """
        PURPOSE: Offset calculation as explained in PLATO-LESIA-PDC-DD-005 
                 (PLATO: N_DPU Onboards Offset Calculation ATBD).
                 Algorithm name: ONB-OFFCAL-010.

        INPUT:
            - offset_rows: Serial pre-scan that is used to calculate the offset [ADU]
    
        OUTPUT:
            - offset_value_fc: Mean of the values in the serial pre-scan after discarding the 
                               outliers [ADU]
            - offset_variance_fc: Variance of the values in the serial pre-scan after discarding 
                                  the outliers [ADU]
        """

        if self.offset_outlier_detection_enabled:

            # Outlier detection

            flag = self.offset_outliers_detection(self.offset_outliers_detection)
    
            # Shifted data algorithm

            offset_value_fc, offset_variance_fc = shifted_data_algorithm(offset_rows[~flag])[:2]
        
        else:

            offset_value_fc, offset_variance_fc = shifted_data_algorithm(offset_rows)[:2]

        return offset_value_fc, offset_variance_fc


    


    def offset_outliers_detection(self, offset_rows):
  
        """
        PURPOSE: Outlier detection in the serial pre-scan as explained in PLATO-MPSSR-PDC-DD-0002
                 (PLATO On-Board Offset & Prescan Outlier Detection Algorithm Theoretical Baseline
                 Document).
                Algorithm name: ONB-OFFOUTDET-010.
    
        INPUT:
            - offset_rows: Serial pre-scan that is used to calculate the offset [ADU]
    
        OUTPUT:
            - outliers_offset_array: Boolean truth values, where 0 means "no outlier" and 1 means
                                     "outlier".  The truth value 1 is equivalent to a flag.
        """

        offset_rows_1d = np.ravel(offset_rows)      # 2D -> 1D
        np.sort(offset_rows_1d)                     # Sort

        if len(offset_rows) < 2 * self.offset_outlier_detection_k:
            raise Exception("Number of entries, {0}, in the serial pre-scan (bias register map) should exceed 2k = {1}".format(len(offset_rows_1d), 2 * self.offset_outlier_detection_k))

        # Flag the k smallest values and the k largest values

        minOffset = offset_rows_1d[self.offset_outlier_detection_k]
        maxOffset = offset_rows_1d[-(self.offset_outlier_detection_k + 1)]

        outliers_offset_array = np.logical_or(offset_rows < minOffset, offset_rows > maxOffset)

        return outliers_offset_array





    ######################
    # Smearing calculation
    ######################

    def smearing_calculation(self, smearing_rows, offset_value_fc, std_dev_previous):

        """
        PURPOSE: Smearing calculation as explained in PLATO-LESIA-PDC-TN- (Parallel overscan rows:
                 correction of the CTI) and PLATO-LESIA-PDC-DD-006 (PLATO: N-DPU Onboard Smearing
                 Calculation ATBD). Coefficient a0 will be updated for the next exposure
                 Algorithm name: ONB-SMRCAL-010.

        INPUT:
            - smearing_rows: Parallel over-scan that is used to calculate the smearing [ADU]
            - offset_value_fc: Electronic offset as calculated from the serial pre-scan [ADU]
            - std_dev_previous: Standard deviation of the previous measurement of the parallel over-scan
                                (one entry per column of the parallel over-scan) [e-]

        OUTPUT:
            - smearing_pattern_fc: One smearing row for this CCD half [e-]
            
            - std_dev_previous: Standard deviation of the current measurement of this column of the parallel 
                               over-scan
        """

        n1 = smearing_rows.shape[0]

        smearing_rows = (smearing_rows - offset_value_fc) * self.half_ccd_gain   # [ADU] -> [e-]

        Ic = np.zeros(smearing_rows.shape)
        smearing_pattern_fc = np.zeros(smearing_rows.shape[1])
    
        a1, a2, a3 = self.smearing_a_coefficients[0], self.smearing_a_coefficients[1], self.smearing_a_coefficients[2]                            # a0 is not in the array (will be updated)
        b0, b1, b2, b3 = self.smearing_b_coefficients[0], self.smearing_b_coefficients[1], self.smearing_b_coefficients[2], self.smearing_b_coefficients[3]
        u0, u1, u2, u3 = math.exp(-b0), math.exp(-b1), math.exp(-b2), math.exp(-b3)                     # Eq. (12) in PLATO-LESIA-PDC-TN-

        for column in range(smearing_rows.shape[1]):
        
            # CTI correction

            f0, f1, f2, f3 = 1, 1, 1, 1     # Will be updated iteratively -> initialisation for i = 0
            tau = (1 + a1 + a2 + a3)        # Will be updated iteratively -> initialisation for i = 0 -> Eq. (3) in PLATO-LESIA-PDC-TN-

            for i in range(n1):

                if i >= self.smearing_n0:

                    # Correct the value

                    Ic[i][column] = smearing_rows[i][column] - self.smearing_a0_array[column] * tau
                
                # Update f0, f1, f2, f3, and tau
            
                f0 *= u0
                f1 *= u1
                f2 *= u2
                f3 *= u3

                tau = f0 + a1 * f1 + a2 * f2 + a3 * f3

            # Outlier detection
        
            if self.smearing_outlier_detection_enabled:
                
                flag = self.smearing_outlier_detection(Ic[self.smearing_n0:n1][column], std_dev_previous[column])
                smearing_pattern_fc[column], std_dev_previous[column] = shifted_data_algorithm(Ic[self.smearing_n0:][column][~flag])[:2]

                if(np.sum(flag) == np.size(flag)):

                    smearing_pattern_fc[column] = 0
                    continue
            else:

                # Calculation of the mean smearing (first n0 measurements are excluded)
                
                smearing_pattern_fc[column], std_dev_previous[column] = shifted_data_algorithm(Ic[self.smearing_n0:][column])[:2]

            # Update a0 for the current column

            chi, rho = 0, 0
            f0, f1, f2, f3 = 1, 1, 1, 1 # Will be updated iteratively -> initialisation for i = 0
            tau = (1 + a1 + a2 + a3)    # Will be updated iteratively -> initialisation for i = 0 -> Eq. (3) in PLATO-LESIA-PDC-TN-

            for i in range(n1):
            
                if (i >= self.smearing_n0) and (flag[i - self.smearing_n0] == 0):
                    chi += (smearing_rows[i][column] - smearing_pattern_fc[column]) * tau
                    rho += pow(tau, 2)

                # Update f0, f1, f2, f3, and tau
            
                f0 *= u0
                f1 *= u1
                f2 *= u2
                f3 *= u3

                tau = f0 + a1 * f1 + a2 * f2 + a3 * f3

            self.smearing_a0_array[column] = (chi + self.smearing_epsilon * rho * self.smearing_a0_array[column]) / (rho * (1 + self.smearing_epsilon))
    
        return smearing_pattern_fc, std_dev_previous





    def smearing_outlier_detection(self, cti_corrected_column, std_dev_previous):

        """
        PURPOSE: Outlier detection in the parallel over-scan as explained in PLATO-MPSSR-PDC-PT-0003
                 (PLATO Onboards Overscan Outlier Detection Algorithm Theoretical Baseline Document).
                 Algorithm name: ONB-OVEROUTDET-010.
    
        INPUT:
            - cti_corrected_column: Column from the parallel over-scan, after CTI correction (and discarding
                                    rows to avoid contamination by bright sources at the top of the detector)
            - std_dev_previous: Standard deviation of the previous measurement of this column of the parallel 
                                over-scan
    
        OUTPUT:
            - raw_overscan_flags: Boolean truth values, where 0 means "no outlier" and 1 means "outlier".  The 
                                  truth value 1 is equivalent to a flag.
        """

        median = np.median(cti_corrected_column)

        raw_overscan_flags = np.ones(len(cti_corrected_column))
        raw_overscan_flags[cti_corrected_column - median >= self.smearing_outlier_detection_threshold * std_dev_previous]

        return raw_overscan_flags                                           





    ########################################################
    # Flux & COB calculations using nominal & extended masks
    ########################################################

    def flux_cob_calculation(self, exposure, offset_value_fc,smearing_pattern_fc, nmask, emask):

        """
        PURPOSE: Flux and COB calculation as explained in PLATO-LESIA-PDC-DD-008 
                (PLATO Onboards Flux & COB Calculation ATBD).
                Algorithm name: ONB-FXCOBCAL-010.

        INPUT:
            - offset_value_fc: Electronic offset as calculated from the serial pre-scan
            - smearing_pattern_fc: Smearing as calculated from the parallel over-scan
            - nmask: Nominal mask (computed on-ground and uploaded)
            - emask: Extended mask (computed on-ground and uploaded)

        OUTPUT:
            - fx_fc: Flux computed using the nominal mask
            - dfx_fc: Flux difference between the extended and the nominal mask
            - ncob_fc: COB computed using the nominal mask
            - ecob_fc: COB computed using the extended mask
        """

        star_window = self.simFile.getImage(exposure)[self.window_zp_row - self.subfield_zp_row : self.window_zp_row + self.window_num_rows - self.subfield_zp_row][self.window_zp_column - self.subfield_zp_column : self.window_zp_column + self.window_num_columns - self.subfield_zp_column]
        window_smearing_pattern_fc = window_smearing_pattern(smearing_pattern_fc)

        rowRange = np.arange(star_window.shape[0])
        columnRange = np.arange(star_window.shape[1])

        pixel_center_offset = 0.5

        # Variable name conventions from Sect. 3.6.1 in PLATO-LESIA-PDC-DD-008

        nmask_column_integral = np.sum(nmask, axis = 0)
        emask_column_integral = np.sum(emask, axis = 0)

        nmask_integral = np.sum(nmask_column_integral)
        emask_integral = np.sum(emask_column_integral)

        nmask_smearing_integral = np.dot(window_smearing_pattern_fc, nmask_column_integral)
        emask_smearing_integral = np.dot(window_smearing_pattern_fc, emask_column_integral)

        nmasked_flux_pixel = np.multiply(star_window, nmask)
        emasked_flux_pixel = np.multiply(star_window, emask)

        nmask_flux_column_integral = np.sum(nmasked_flux_pixel, axis = 0)
        emask_flux_column_integral = np.sum(emasked_flux_pixel, axis = 0)

        nmasked_flux_integral = np.sum(nmask_flux_column_integral)
        emasked_flux_integral = np.sum(emask_flux_column_integral)

        nmasked_Xcob_sum = np.dot(columnRange, nmask_flux_column_integral)
        emasked_Xcob_sum = np.dot(columnRange, emask_flux_column_integral)

        nmasked_Ycob_sum = np.dot(np.sum(nmasked_flux_pixel, axis = 1), rowRange)
        emasked_Ycob_sum = np.dot(np.sum(emasked_flux_pixel, axis = 1), rowRange)

        nmask_column_Xweighted = np.multiply(columnRange, nmask_column_integral)
        emask_column_Xweighted = np.multiply(columnRange, emask_column_integral)

        nmask_column_Yweighted = np.sum(np.multiply(nmask, columnRange), axis = 0)
        emask_column_Yweighted = np.sum(np.multiply(emask, columnRange), axis = 0)

        nmask_Xweighted_integral = np.sum(nmask_column_Xweighted)
        emask_Xweighted_integral = np.sum(emask_column_Xweighted)

        nmask_Yweighted_integral = np.sum(nmask_column_Yweighted)
        emask_Yweighted_integral = np.sum(emask_column_Yweighted)

        nmask_smearing_Xweighted_integral = np.dot(window_smearing_pattern, nmask_column_Xweighted)
        emask_smearing_Xweighted_integral = np.dot(window_smearing_pattern, emask_column_Xweighted)

        nmask_smearing_Yweighted_integral = np.dot(window_smearing_pattern, nmask_column_Yweighted)
        emask_smearing_Yweighted_integral = np.dot(window_smearing_pattern, emask_column_Yweighted)

        # Calculate the flux in the nominal mask & the difference in flux between the extended and the nominal mask

        fx_fc = (nmasked_flux_integral - offset_value_fc * nmask_integral) * self.half_ccd_gain - nmask_smearing_integral    # Flux in nominal mask
        dfx_fc = (emasked_flux_integral - offset_value_fc * emask_integral) * self.half_ccd_gain - emask_smearing_integral   # Flux difference between extended & nominal mask

        # Calculate the COB for the nominal and for the extended mask

        nmasked_Xcob_sum_corrected = (nmasked_Xcob_sum - offset_value_fc * nmask_Xweighted_integral) * self.half_ccd_gain - nmask_smearing_Xweighted_integral
        nmasked_Ycob_sum_corrected = (nmasked_Ycob_sum - offset_value_fc * nmask_Yweighted_integral) * self.half_ccd_gain - nmask_smearing_Yweighted_integral

        ncob_fc = [nmasked_Xcob_sum_corrected / fx_fc + pixel_center_offset, nmasked_Ycob_sum_corrected / fx_fc + pixel_center_offset]

        emasked_Xcob_sum_corrected = (emasked_Xcob_sum - offset_value_fc * emask_Xweighted_integral) * self.half_ccd_gain - emask_smearing_Xweighted_integral
        emasked_Ycob_sum_corrected = (emasked_Ycob_sum - offset_value_fc * emask_Yweighted_integral) * self.half_ccd_gain - emask_smearing_Yweighted_integral

        masked_Xcob_sum_corrected = nmasked_Xcob_sum_corrected + emasked_Xcob_sum_corrected
        masked_Ycob_sum_crrected = nmasked_Ycob_sum_corrected + emasked_Ycob_sum_corrected
        efx_fc = fx_fc + dfx_fc

        ecob_fc = [masked_Xcob_sum_corrected / efx_fc + pixel_center_offset, masked_Ycob_sum_crrected / efx_fc + pixel_center_offset]
    
        return fx_fc, dfx_fc, ncob_fc, ecob_fc





    #########################
    # Smearing time averaging
    #########################

    def smearing_time_averaging(self):

        """
        PURPOSE: Smearing time averaging for long cadence as explained in PLATO-LESIA-PDC-DD-007
                 (PLATO Ob-board Smearing Pattern Time Averaging ATBD).
                 Algorithm name: ONB-SMRAVG-010
    
        OUTPUT:
            - smearing _pattern_lc: Smearing pattern averaged over the last 24 samples (i.e. 600s) [e-]
        """

        smearing_pattern_fc_array_last = self.smearing_pattern_fc_array[-int(self.num_exposures_lc) :]     # Last 24 samples

        shape = (1, self.num_exposures_lc, smearing_pattern_fc_array_last.shape[1], 1)
        smearing_pattern_lc = smearing_pattern_fc_array_last.reshape(shape).mean(-1).mean(1)

        return smearing_pattern_lc




    
    ####################################
    # Outlier detection over light curve
    ####################################

    def flux_cob_outlier_detection(self, fx_fc_array, fx_exposure_error_array):

        """
        PURPOSE: Outlier detection over light curve as explained in PLATO-MPSSR-PDC-DD-0001
                 (PLATO Onboard LC Outlier Detection ATBD).
                 Algorithm name: ONB-LCOUTDET-010.

        INPUT:
            - flux_fc_array: Flux time series processed so far. This can either be the flux 
                             obtained in the nominal mask (fx_fc) or the difference in flux 
                             between the extended and the nominal mask (dfx_fc).
            - fx_exposure_error_array: Flag for the time series processed so far, apart from
                                       the last 2 * b + 1 samples
        
        OUTPUT:
            - fx_exposure_error_array: Flag for the time series processed so far, including
                                       the last 2 * b + 1
        """

        fx_fc_array_last = fx_fc_array[-2 * self.lc_outlier_detection_b - 1:]

        # Step 1 - 4

        is_outlier = mad_median_clipping(fx_fc_array_last, self.lc_outlier_detection_threshold, self.lc_outlier_detection_b)

        # Step 5

        if fx_exposure_error_array[-2] and fx_exposure_error_array[-1]:

            fx_exposure_error_array[-2] = False

            if not is_outlier:

                fx_exposure_error_array[-1] = False

        fx_exposure_error_array = np.append(fx_exposure_error_array, is_outlier)

        return fx_exposure_error_array





    ##########################
    # Flux & COB time averaging
    ###########################

    def flux_cob_time_averaging(self):

        """
        PURPOSE: Flux time averaging for short cadence (50s) as explained in PLATO-LESIA-PDC-DD-009
                 (PLATO Onboard Flux and COB Short Cadence Time Averaging ATBD) or long cadence as 
                 explained in PLATO-LESIA-PDC-DD-010 (PLATO Onboard Flux and COB Long Cadence Time Averaging ATBD).
                 Algorithm names: ONB-FXAGV-011 and ONB-FXAGV-012.
        """

        # Short cadence

        if(self.time_averaging_cadence == "Short"):
            
            # Nominal mask

            fx_sc, ncob_sc, n_useful = self.flux_cob_time_averaging_sc(self.fx_fc_array, self.ncob_fc_array, self.nflag_fc)           # Select duration of short cadence (50s)

            if n_useful != 0:

                self.fx_sc_array = np.append(self.fx_sc_array, fx_sc)
                self.ncob_sc_array = np.append(self.ncob_sc_array, ncob_sc)

            # Extended mask

            dfx_sc, ecob_sc, n_useful = self.flux_cob_time_averaging_sc(self.dfx_fc_array, self.ecob_fc_array, self.eflag_fc)         # Select duration of short cadence (50s)

            if n_useful != 0:

                self.dfx_sc_array = np.append(self.dfx_sc_array, dfx_sc)
                self.ecob_sc_array = np.append(self.ecob_sc_array, ecob_sc)

        

        # Long cadence

        elif self.time_averaging_cadence == "Long":

            # Nominal mask

            fx_lc, fxvar_lc, n_useful = self.flux_cob_time_averaging_lc(self.fx_fc_array, self.nflag_fc)         # Select duration of long cadence (600s)

            if n_useful != 0:
            
                self.fx_lc_array = np.append(self.fx_lc_array, fx_lc)
                self.fxvar_lc_array = np.append(self.fxvar_lc_array, fxvar_lc)
            
            # Extended mask

            dfx_lc, dfxvar_lc, n_useful = self.flux_cob_time_averaging_lc(self.dfx_fc_array, self.eflag_fc)      # Select duration of long cadence (600s)
            
            if n_useful != 0:

                self.dfx_lc_array = np.append(self.dfx_lc_array, dfx_lc)
                self.dfxvar_lc_array = np.array(self.dfxvar_lc_array, dfxvar_lc)





    def flux_cob_time_averaging_sc(self, fx_fc_array, cob_fc_array, fx_exposure_error_array = None):

        """
        PURPOSE: Flux time averaging for short cadence (50s) as explained in PLATO-LESIA-PDC-DD-009
                 (PLATO Onboard Flux and COB Short Cadence Time Averaging ATBD).
                 Algorithm name: ONB-FXAGV-011.

        INPUT:
            - fx_fc_array: Flux time series processed so far. This can either be the flux 
                           obtained in the nominal mask (fx_fc) or the difference in flux 
                           between the extended and the nominal mask (dfx_fc).
            - cob_fc_array: COB time series processed so far.  This can either be the COB
                            obtained from the nominal mask or from the extended mask.
            - fx_exposure_error_array: Outlier detection flags for the flux time series 
                                       processed so far.  If None, outlier detection on the 
                                       light curve was disabled.

        OUTPUT:
            - fx_sc: Mean of the given flux time series over the last 50s, calculated with the non-flagged datapoints only
            - cob_sc: Mean of the given COB time series over the last 50s, calculated with the non-flagged datapoints only
            - n_useful_sc: Number of non-flagged datapoints in the time series
        """

        # Select the last 50s (i.e. duration of short cadence)

        fx_fc_array_last = fx_fc_array[-int(self.num_exposures_sc) :]
        cob_fc_array_last = cob_fc_array[-int(self.num_exposures_sc) :]

        # Outlier detection enabled

        if self.lc_outlier_detection_enabled:

            fx_exposure_error_array_last = fx_exposure_error_array[-int(self.num_exposures_sc) :]
            n_useful_sc = len(fx_exposure_error_array_last) - np.sum(fx_exposure_error_array_last)

            fx_sc = np.sum(fx_fc_array_last[~fx_exposure_error_array_last]) / n_useful_sc 

            cob_sc_x = np.sum(cob_fc_array_last[~fx_exposure_error_array_last][0]) / n_useful_sc
            cob_sc_y = np.sum(cob_fc_array_last[~fx_exposure_error_array_last][1]) / n_useful_sc
            cob_sc = [cob_sc_x, cob_sc_y]

        # Outlier detection disabled

        else:

            n_useful_sc = self.num_exposures_sc

            fx_sc = np.sum(fx_fc_array_last) / n_useful_sc

            cob_sc_x = np.sum(cob_fc_array_last[0]) / n_useful_sc
            cob_sc_y = np.sum(cob_fc_array_last[1]) / n_useful_sc
            cob_sc = [cob_sc_x, cob_sc_y]

        return fx_sc, cob_sc, n_useful_sc





    def flux_cob_time_averaging_lc(self, fx_fc_array, fx_exposure_error_array = None):

        """
        PURPOSE: Flux time averaging for long cadence (600s) as explained in PLATO-LESIA-PDC-DD-010
                 (PLATO Onboard Flux and COB Long Cadence Time Averaging ATBD).
                 Algorithm name: ONB-FXAGV-012.

        INPUT:
            - fx_fc_array: Flux time series processed so far. This can either be the flux 
                           obtained in the nominal mask (fx_fc) or the difference in flux 
                           between the extended and the nominal mask (dfx_fc).
            - fx_exposure_error_array: Outlier detection flags for the flux time series 
                                       processed so far.  If None, outlier detection on the 
                                       light curve was disabled.
        OUTPUT:
            - fx_lc: Mean of the given flux time series over the last 600s, calculated with the non-flagged datapoints only
            - fxvar_lc: Variance of the given flux time series over the last 600s, calculated with the non-flagged datapoints only
            - n_useful_lc: Number of non-flagged datapoints in the time series
        """

        # Select the last 600s (i.e. duration of long cadence)

        fx_fc_array_last = fx_fc_array[-int(self.num_exposures_lc) :]

        # Outlier detection enabled

        if self.lc_outlier_detection_enabled:

            fx_exposure_error_array_last = fx_exposure_error_array[-int(self.num_exposures_lc) :]
            
            fx_lc, fxvar_lc, n_useful_lc = shifted_data_algorithm(fx_fc_array_last[~fx_exposure_error_array_last])

        # Outlier detection disabled

        else:
            
            fx_lc, fxvar_lc, n_useful_lc = shifted_data_algorithm(fx_fc_array_last)

        return fx_lc, fxvar_lc, n_useful_lc
