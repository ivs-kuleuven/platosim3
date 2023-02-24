"""
Run the PLATO Simulator from Python.

The Simulation class provides the opportunity to interactively tune the input parameters
before the simulator is started. The parameters that are available can be inspected by
just printing the Simulation object, i.e. print (sim), which will dump all the parameters
and their current values on the command line.

For usage, see the tutorial Jupyter-notebooks available at "PlatoSim/docs/tutorials".
"""

import os
import sys
import ast
import math
import yaml
import pyaml
import inspect
import datetime
import subprocess
import numpy as np

import platosim.referenceFrames as rf
import platosim.instrument      as it
from platosim.simfile import SimFile


#==============================================================#
#                         BEGIN CLASS                          #
#==============================================================#


class Simulation(object):

    """Class for running PlatoSim simulations.

    Simulation class allows running the PLATO simulator interactively from Python
    and tuning the input parameters before each run. For more help, type:
    
    Example
    -------
    >>> import platosim.simulation as Simulation
    >>> print(Simulation)
    """


    def __init__(self, runName, configurationFile=None, outputDir=None, debug=False):

        """Initialise class variables.
        """

        self.debug   = debug
        self.runName = runName

        # Glag to check if output directory has been specified before running the Simulation

        self.hasTargetLocation = False
        self.targetOutputFilesLocation = None

        # Set output directory
        
        if outputDir is not None:
            self.outputDir = outputDir

        # Set simulation location
        
        self.setSimulatorLocation();

        # Read the YAML input file
        
        if configurationFile:
            self.readConfigurationFile( configurationFile )
        else:
            self.readConfigurationFile( self.originalInputFilesLocation + "/inputfile.yaml" )





    def setSimulatorLocation(self):

        """Set the location of the simulation.

        Given the location of the simfile.py module, try to find the build/ folder where
        the simulator executable should be. Then set the default locations for the platosim
        executable (i.e. build directory), and the original input files location.
        """

        # Find the absolute path of the simfile.py. This will allow us to located the other
        # default project directories. Build in a test to check that the build sub-directory
        # exists, otherwise this is probably not a correct default installation.

        path = os.environ["PLATO_PROJECT_HOME"]

        self.platosimLocation = path

        if os.path.exists(path + "/build"):
            self.platosimBuildLocation = self.platosimLocation + "/build"

        elif os.path.exists(path + "/bin"):
            self.platosimBuildLocation = self.platosimLocation + "/bin"

        else:
            raise Exception("Unexpected directory structure for this PLATO Simulator" +
                            "distribution: no build nor bin sub-directory in PLATO_PROJECT_HOME")

        # This is the location of the original input files as distributed by the PLATO Simulator

        self.originalInputFilesLocation  = self.platosimLocation + "/inputfiles"
        self.originalOutputFilesLocation = self.platosimLocation + "/outputfiles"





    @property
    def outputDir(self):

        """Return the output files location.
        """
        return self.targetOutputFilesLocation




    
    @outputDir.setter
    def outputDir(self, path):

        """Specify the absolute path for the output directory. 

        This directory will contain a copy of the modified input file and
        the HDF5 output file for this simulation. If the output path do not
        exist, it is created.
        """

        if not os.path.exists(path):
            if self.debug:
                print(f"DEBUG: creating output directory {path}")
            self.createDirectory(path)

        self.targetOutputFilesLocation = path
        self.hasTargetLocation = True

        if self.debug:
            print(f"DEBUG: output dir set to {path}")





    def readConfigurationFile(self, filename):

        """Read the YAML input configuration file.
        """
        
        self.configurationFilename = filename

        if self.debug:
            print(f"DEBUG: Parsing YAML configuration file {filename}")

        with open(filename, 'r') as stream:
            try:
                self.yamlDocument = yaml.safe_load(stream)
            except yaml.YAMLError as exc:
                print(exc)





    def getYamlConfiguration(self):

        """Return the YAML configuration as a dictionary.
        """
        
        return self.yamlDocument





    def __contains__(self, key):

        """Returns true if the input parameter (key) is known/exists.

        Usage: if "Group/ParameterName" in sim: # do something

        Parameters
        ----------
        key : str
            A string containing the parameter name or "Group/ParameterName" combination.

        Return
        ------
        bool : True if key exists, False otherwise
        """
        
        if key.find('/') == -1:
            nodeNames = [key]
        else:
            nodeNames = key.split("/")

        node = self.yamlDocument

        for nodeName in nodeNames:
            print(f"> {nodeName}, {type(node)}")
            try:
                node = node[nodeName]
            except:
                return False

        return True





    def __getitem__(self, key):

        """Returns the value of the input parameter (key).

        Parameters
        ----------
        key : str
            A string containing the parameter name or "Group/ParameterName" combination.

        Return
        ------
        value : int, float, ndarray, str
            The value of the requested parameter.
        """

        # Split the path into node names
        # E.g. "PSF/MappedGaussian/Sigma" into [PSF, MappedGaussian, Sigma]

        if key.find('/') == -1:
            parentNodeName, nodeName = key, None
            print("usage: the given parameter name (key) should " +
                  "include the group name of the group that contains the parameter.")
            print("E.g in 'Camera/PlateScale', Camera is the group, PlateScale is the parameter.")
            return None
        else:
            nodeNames = key.split("/")

        # Navigate to the deepest node, starting from the document root

        node = self.yamlDocument

        for nodeName in nodeNames:
            if nodeName in node:
                node = node[nodeName]
            else:
                print(f"ERROR: The group '{key}' was not found in the " +
                      f"yaml inputfile '{self.configurationFilename}'"
                return None

        # Node is a string, so cast it to its proper value

        try:
            value = ast.literal_eval(node)
        except ValueError:
            value = node

        # Return the value of the deepest node

        return value





    def __setitem__(self, key, item):

        """Update a specific node.

        Parameters
        ----------
        key : str
            A string with parent node name and node name seperated by a slash.
        item : str
            A string with the new node value, if not a string the value is converted using str().

        Return
        ------
        bool : True if node could be updated, False otherwise
        """

        # Ensure that the given item is a string

        if type(item) is np.ndarray:
            item = list(item)

        if type(item) != list:
            item = str(item)

        # Split the path into node names
        # E.g. "PSF/MappedGaussian/Sigma" into ["PSF", "MappedGaussian", "Sigma"]

        if key.find('/') == -1:
            print ("USAGE: the given parameter name (key) should include the " +
                   "group name of the group that contains the parameter.")
            print ("       E.g in 'Camera/PlateScale', Camera is the group," +
                   "PlatScale is the parameter.")
            return None
        else:
            nodeNames = key.split("/")

        # Check whether the parent node is in the document. If not, complain

        if nodeNames[0] not in self.yamlDocument:
             print(f"ERROR: no node with the name {nodeNames[0]} found in input yaml file")
             return False

        # If there is only 1 node in the path, we're finished after setting its value

        if len(nodeNames) == 1:
            self.yamlDocument[nodeNames[0]] = item
            return True

        # If we arrive here, there are at least 2 node in the path, check if 2nd parent node exists

        if nodeNames[1] not in self.yamlDocument[nodeNames[0]]:
             print("ERROR: no node with the name " +
                   f"{nodeNames[0]}/{nodeNames[1]} found in input yaml file")
             return False

        # If there are only 2 nodes in the path, we're finished after setting its value

        if len(nodeNames) == 2:
            self.yamlDocument[nodeNames[0]][nodeNames[1]] = item
            return True

        # If we arrive here, there are at least 3 nodes in the path, check if 3rd parent node exists

        if nodeNames[2] not in self.yamlDocument[nodeNames[0]][nodeNames[1]]:
             print("ERROR: no node with the name " +
                   f"{nodeNames[0]}/{nodeNames[1]}/{nodeNames[2]} found in input yaml file")
             return False

        # If there are only 3 nodes in the path, we're finished after setting its value

        if len(nodeNames) == 3:
            self.yamlDocument[nodeNames[0]][nodeNames[1]][nodeNames[2]] = item
            return True

        # If we arrive here, there are at least 4 nodes in the path, check if 4th parent node exists

        if nodeNames[3] not in self.yamlDocument[nodeNames[0]][nodeNames[1]][nodeNames[2]]:
             print("ERROR: no node with the name " +
                   f"{nodeNames[0]}/{nodeNames[1]}/{nodeNames[2]}/{nodeNames[3]} " +
                   "found in input yaml file")
             return False

        # If there are only 34nodes in the path, we're finished after setting its value

        if len(nodeNames) == 4:
            self.yamlDocument[nodeNames[0]][nodeNames[1]][nodeNames[2]][nodeNames[3]] = item
            return True

        # If we arrive here, there are at least 5 nodes in the path.
        # Issue a not-implemented error message.

        print(f"ERROR: detected 5 or more nodes in the path {key}")
        return False





    def createDirectory(self, path):

        """Create a directory.
        """

        try:
            os.makedirs(path)
        except OSError as ose:
            print (ose)
            if not os.path.isdir(path):
                raise Exception(f"Could not create directory {path}")

        return





    def writeYamlConfigurationFile(self, filename):

        """Write the modified configuration to output file location. 

        This configuration will be loaded by the PLATO Simulator when
        the run() method is executed.

        Parameter
        ---------
        filename : str
            Filename of the output YAML file.
        """
        
        if self.debug:
            print(f"DEBUG: writing the YAML configuration file {filename}")
        with open(filename, 'w') as outfile:
            outfile.write( pyaml.dump(self.yamlDocument, indent=4, width=120) )





    def run(self, removeOutputFile=False, executionTime=False, logLevel=3):

        """Run the PLATO Simulator.

        Parameters
        ----------
        removeOutputFile : bool
            If the outputfile already exists before the run started, simply delete it.
        logLevel : int
            Level of verbosoty: 1 (least verbose) to 3 (most verbose)
        
        Return
        ------
        When PlatoSim fails for some reason and returns an error code (!= 0),
        an Exception is raised.
        """

        if executionTime:
            tic = datetime.datetime.now()

        if not self.hasTargetLocation:
            raise Exception("Output location not set for this Simulation. " +
                            "Set the outputDir before executing the run() method.")

        inputFilename  = f"{self.targetOutputFilesLocation}/{self.runName}.yaml"
        outputFilename = f"{self.targetOutputFilesLocation}/{self.runName}.hdf5"
        logFilename    = f"{self.targetOutputFilesLocation}/{self.runName}.log"

        if removeOutputFile:
            try:
                os.remove(outputFilename)
            except OSError:
                pass

        self.writeYamlConfigurationFile(inputFilename)

        # The run() method was only introduced with Python 3.5.
        # Use the older call() method when running e.g. Python 2.7

        if sys.version_info < (3, 5):
            rc = subprocess.call([self.platosimBuildLocation + "/platosim", inputFilename,
                                  outputFilename, logFilename, str(logLevel)])
            if rc:
                raise Exception(f"Simulation.run(): PlatoSim returned with exit code {rc}.")
        else:
            completedProcess = subprocess.run([self.platosimBuildLocation + "/platosim",
                                               inputFilename, outputFilename, logFilename,
                                               str(logLevel)],
                                               stdout=subprocess.PIPE, stderr=subprocess.PIPE)

            # print(str(completedProcess.stdout.decode("utf-8")))
            # print(str(completedProcess.stderr.decode("utf-8")))

            if completedProcess.returncode:
                raise Exception("Simulation.run(): PlatoSim returned with " +
                                "exit code {completedProcess.returncode}.")

        # Print computation time
        
        if executionTime:
            toc = datetime.datetime.now()
            print(f"Execution time : {toc - tic} [hh:mm:ss]")

        return SimFile(outputFilename)





    def __str__(self):

        """Return a listing of all settings.
        
        This function fetch a listing of all settings from both the CCD
        and the Photometry input files. If a parameter value has been
        updated for this Simulation, [updated] will be printed after the
        value.
        """
        
        root = self.yamlDocument
        msg = "YAML Configuration:\n"
        msg += pyaml.dump(root, indent=4)

        return msg





    def writeAllOutputToHDF5(self, write=True):

        """Function to write all or nothing to the HDF5 file.
        """

        # Fetch names of ControlHDF5Content attributes

        group = "ControlHDF5Content"
        entries = []
        for name, dict_ in self.yamlDocument[group].items():
            entries.append(name)

        # Control the content

        for entry in entries:

            # Set all HDF5 content parameters to "yes"

            if write:
                self.__setitem__(f"{group}/{entry}", "yes")

            # Set all HDF5 content parameters to "no"

            elif write is False:
                self.__setitem__(f"{group}/{entry}", "no")

            else:
                print("ERROR: only 'True' or 'False' can be parsed as argument!")





    def useNominalCamera(self):

        """Change the input parameters to use the nominal camera's.

        The following parameters are updated:
        
            CCD/NumColumns = 4510
            CCD/NumRows    = 4510
            ObservingParameters/CycleTime = 25
        """

        self.__setitem__("CCD/NumColumns", "4510")
        self.__setitem__("CCD/NumRows",    "4510")
        self.__setitem__("ObservingParameters/CycleTime", "25")

        return





    def useFastCamera(self):

        """Change the input parameters to use the fast camera's.

        The following parameters are updated:

            CCD/NumColumns = 4510
            CCD/NumRows    = 2255
            ObservingParameters/CycleTime    = 2.5
            ObservingParameters/ExposureTime = 2.3
        """

        self.__setitem__("CCD/NumColumns", "4510")
        self.__setitem__("CCD/NumRows",    "2255")
        self.__setitem__("ObservingParameters/CycleTime", "2.5")

        return





    def setSubfieldAroundPixelCoordinates(self, ccdCode, xCCDpixel, yCCDpixel,
                                          subfieldSizeX, subfieldSizeY):

        """Set the subfield around pixel coordinates.
        
        This function calculate the location of the subField such that it
        is centered on the star with the given pixel coordinates.
        
        Parameters
        ----------
        ccdCode : str
            For nominal camera: either '1', '2', '3', '4'
            For fast camera: either '1F', '2F', '3F', '4F'
        xCCDpixel : int
            X-coordinate (column-number) of the star on the CCD [pixel/float]
        yCCDpixel : int
            Y-coordinate (row-number) of the star on the CCD [pixel/float]
        subfieldSizeX : int
            Width (i.e. number of columns) of the sub-field [pixels]
        SubfieldSizeY : int
            Feight (i.e. number of rows) of the sub-field [pixels]

        Return
        ------
        None
        """

        raStar, decStar = rf.pixelToSkyCoordinates(self, ccdCode, xCCDpixel, yCCDpixel)

        # TODO: determine nominal from the given ccdCode

        nominal = True

        success = self.setSubfieldAroundCoordinates(raStar, decStar,
                                                    subfieldSizeX, subfieldSizeY, True)

        if not success:
            print("WARNING: setSubfieldAroundPixelCoordinates() " +
                  "failed to set subField around the star.")

        return





    def setSubfieldAroundPixelRow(self, ccdCode, yCCDpixel, subfieldSizeY):

        """Set subfield around a pixel row.

        This function sets the location of the subfield so that its rows are
        centered around the given pixel coordinate.
        
        Paramters
        ---------
        ccdCode : str
            For nominal camera: either '1', '2', '3', '4'.
            For fast camera: either '1F', '2F', '3F', '4F'.
        yCCDpixel : int
            Y-coordinate (row-number) around which to center the subField  [pixel/float].
        subfieldSizeY : int
            Height (i.e. number of rows) of the sub-field [pixels].

        Return
        ------
        None

        Rmarks
        ------
        If it is not possible to center the subfield around the row, because the row
        is too close the the edge, the subfield will be set at that respective edge.
        """

        subfieldSizeX = int(self["SubField/NumColumns"])
        subfieldRowZero = int(self["SubField/ZeroPointColumn"])
        xCCDpixel = subfieldRowZero + subfieldSizeX / 2

        if not (0 <= yCCDpixel <  4510):
            print("Error: we expect input row coordinate in [0, 4510], " +
                  f"but value {yCCDpixel} was given.")
            return

        if not (1 <= subfieldSizeY <= 4510):
            print("Error: we expect size of the row subfield in [1, 4510], " +
                  f"but value {subfieldSizeY} was given.")
            return

        self.setSubfieldAroundPixelCoordinates(ccdCode, xCCDpixel, yCCDpixel, 1, 1)

        deltaY = int(subfieldSizeY / 2)
        self["SubField/NumColumns"] = subfieldSizeX
        self["SubField/NumRows"]    = subfieldSizeY

        self["SubField/ZeroPointColumn"] = subfieldRowZero
        self["SubField/ZeroPointRow"]    = min(max(0, self["SubField/ZeroPointRow"] - deltaY),
                                               4510 - subfieldSizeY)





    def setSubfieldAroundCoordinates(self, raStar, decStar, subfieldSizeX, subfieldSizeY,
                                     normal=True):

        """Set subfield around stellar coordinates
        
        Set the location of the sub-field such that it is centred on the star
        with the given sky coordinates.  Depending on the CCD (in nomincal mode:
        "1", "2", "3", or "4"; in fast mode: "1F", "2F", "3F", or "4F"), the
        configuration file for the given simulation is adapted.  These include the
        pre-defined CCD position, the dimensions of the CCD (and also of the subfield,
        although this is not affected by the calculations), the sub-field zeropoint
        and the exposure time.

        Notes
        -----
        - This function calls the calculateSubfieldAroundCoordinates() function in
          reference frames.
        - It is assumed that the configuration parameters in the sim object contains
          a correct (ra, dec)  of the platform, a correct (azimuth, tilt) of the telescope,
          a valid values for the focal length, the plate scale, the pixel size, and that
          the switch to include distortion or not is set correctly
        - The function does not set the exposure time, nor the focal length source, etc.

        Paramters
        ---------
        raStar : float
            Right ascension of the star [radians]
        decStar : float 
            Declination of the star [radians]
        subfieldSizeX : int
            Width (i.e. number of columns) of the subiield [pixels]
        subfieldSizeY : int
            Height (i.e. number of rows) of the sub-field [pixels]
        normal : bool
            True for the normal camera configuration, False for the fast cameras

        Return
        ------
        bool : True if the entire subfield fit on one of the 4 (pre-defined) CCDs, 
               False otherwise.
        """

        # Find out some instrumental characteristics from the sim object

        raPlatform  = np.deg2rad(float(self["ObservingParameters/RApointing"]))
        decPlatform = np.deg2rad(float(self["ObservingParameters/DecPointing"]))

        telescopeGroupID = self["Telescope/GroupID"]
        if telescopeGroupID == "Custom":
            azimuthTelescope = np.deg2rad(float(self["Telescope/AzimuthAngle"]))
            tiltTelescope    = np.deg2rad(float(self["Telescope/TiltAngle"]))
        elif telescopeGroupID == "Fast":
            azimuthTelescope = np.deg2rad(self["CameraGroups/AzimuthAngle"][4])
            tiltTelescope    = np.deg2rad(self["CameraGroups/TiltAngle"][4])
        else:
            azimuthTelescope = np.deg2rad(self["CameraGroups/AzimuthAngle"][telescopeGroupID-1])
            tiltTelescope    = np.deg2rad(self["CameraGroups/TiltAngle"][telescopeGroupID-1])

        solarPanelOrientation = np.deg2rad(float(self["Platform/SolarPanelOrientation"])) # [rad]
        focalLength     = float(self["Camera/FocalLength/ConstantValue"]) * 1000.0        # [m]->[mm]
        focalPlaneAngle = np.deg2rad(float(self["Camera/FocalPlaneOrientation/ConstantValue"]))
        pixelSize       = float(self["CCD/PixelSize"])

        # If the psf is MappedFromFile we need to include mapped field distortion

        if self["PSF/Model"] == "MappedFromFile":
            includeFieldDistortion = True
            mappedDistortion       = True
            pathToPsfFile          = self["PSF/MappedFromFile/Filename"]
            distortionCoefficients = None
        elif (self["Camera/IncludeFieldDistortion"] == "yes" or
              self["Camera/IncludeFieldDistortion"] == "1"   or
              self["Camera/IncludeFieldDistortion"] == True):
            includeFieldDistortion = True
            mappedDistortion       = False
            pathToPsfFile          = None
            distortionCoefficients = self["Camera/FieldDistortion/ConstantCoefficients"]
        else:
            includeFieldDistortion = False
            mappedDistortion       = False
            pathToPsfFile          = None
            distortionCoefficients = None

        # Compute the position of the subfield. xPix and yPix are the CCD coordinates
        # of the star, given a 4510x4510 CCD [colNumber, rowNumber]. The function below
        # also checks if the subfield fits entirely on the CCD. If not: ccdCode is None.

        ccdCode, xPix, yPix = rf.calculateSubfieldAroundCoordinates(subfieldSizeX, subfieldSizeY,
                                                                    raStar, decStar,
                                                                    raPlatform, decPlatform,
                                                                    solarPanelOrientation,
                                                                    tiltTelescope, azimuthTelescope,
                                                                    focalPlaneAngle,
                                                                    focalLength, pixelSize,
                                                                    includeFieldDistortion, normal,
                                                                    mappedDistortion,
                                                                    distortionCoefficients,
                                                                    pathToPsfFile)
        if ccdCode == None:
            return False

        CCDSizeX         = rf.CCD[ccdCode]["Ncols"]
        CCDSizeY         = rf.CCD[ccdCode]["Nrows"]
        CCDOriginOffsetX = rf.CCD[ccdCode]["zeroPointXmm"]
        CCDOriginOffsetY = rf.CCD[ccdCode]["zeroPointYmm"]
        CCDOrientation   = rf.CCD[ccdCode]["angle"]

        # If we arrive here, there is no problem accommodating the entire sufield on the CCD

        self["CCD/Position"]      = str(ccdCode)
        self["CCD/OriginOffsetX"] = str(CCDOriginOffsetX)
        self["CCD/OriginOffsetY"] = str(CCDOriginOffsetY)
        self["CCD/Orientation"]   = str(np.rad2deg(CCDOrientation))

        self["CCD/NumColumns"] = CCDSizeX
        self["CCD/NumRows"]    = CCDSizeY

        if telescopeGroupID == "Fast":
            self["CCD/FirstRowExposed"] = str(2255)
        else:
            self["CCD/FirstRowExposed"] = str(0)

        self["SubField/ZeroPointRow"]    = str(yPix - int(subfieldSizeY/2))
        self["SubField/ZeroPointColumn"] = str(xPix - int(subfieldSizeX/2))
        self["SubField/NumRows"]    = str(subfieldSizeY)
        self["SubField/NumColumns"] = str(subfieldSizeX)

        self["Telescope/AzimuthAngle"] = np.rad2deg(azimuthTelescope)
        self["Telescope/TiltAngle"]    = np.rad2deg(tiltTelescope)

        # That's it

        return True





    def createStarCatalogFileFromPixelCoordinates(self, rows, cols, magnitudes,
                                                  starIDs, starCatalogPath):

        """Create a star catalogue file from the pixel coordinates.
        
        Create a star catalog ascii file given the pixel coordinates 
        (row and column) of the stars. This requires the orientation
        of the spacecraft, telescopes, focal plane, hence it's a member
        function of the Simulation class.

        Paramters
        ---------
        rows : ndarray
            Array with fractional row coordinates of the stars (CCD, not subfield) [pix]
        cols : ndarray
            Array with fractional column coordinates of the stars (CCD, not subfield) [pix]
        magnitudes : ndarray
            Array with Johnson V magnitudes of the stars
        starIDs : ndarray
            Array with IDs of the star (integers)
        starCatalogPath : str
            Path of the star catalog file that will be written.

        Return
        ------
        A file will be saved, containing, ra, dec, and magnitude of the stars.
        The "ObservingParameters/StarCatalogFile" tag in the yaml tree will be
        changed to the given starCatalogPath
        """

        # Extract the needed information from the yaml input file
        # Note: groupIDs and ccdIDs start counting from 1...

        if self["Telescope/GroupID"] == "Custom":
            azimuthAngle    = np.deg2rad(self["Telescope/AzimuthAngle"])
            tiltAngle       = np.deg2rad(self["Telescope/TiltAngle"])
        else:
            groupID = int(self["Telescope/GroupID"])
            azimuthAngle    = np.deg2rad(self["CameraGroups/AzimuthAngle"][groupID-1])
            tiltAngle       = np.deg2rad(self["CameraGroups/TiltAngle"][groupID-1])

        if self["CCD/Position"] == "Custom":
            ccdZeroPointX   = self["CCD/OriginOffsetX"]
            ccdZeroPointY   = self["CCD/OriginOffsetY"]
            CCDangle        = np.deg2rad(self["CCD/Orientation"])
        else:
            ccdID = int(self["CCD/Position"])
            ccdZeroPointX   = self["CCDPositions/OriginOffsetX"][ccdID-1]
            ccdZeroPointY   = self["CCDPositions/OriginOffsetY"][ccdID-1]
            CCDangle        = np.deg2rad(self["CCDPositions/Orientation"][ccdID-1])

        pixelSize       = self["CCD/PixelSize"]  # [micron]
        raPlatform      = np.deg2rad(self["ObservingParameters/RApointing"])
        decPlatform     = np.deg2rad(self["ObservingParameters/DecPointing"])
        focalPlaneAngle = np.deg2rad(self["Camera/FocalPlaneOrientation/ConstantValue"])
        focalLength     = self["Camera/FocalLength/ConstantValue"] * 1000.0  # [m] -> [mm]

        if (self["PSF/Model"] == "MappedFromFile"):
            incldueFieldDistortion        = True
            mappedDistortion              = True
            inverseDistortionCoefficients = None
            pathToPsfFile                 = self["PSF/MappedFromFile/Filename"]
        else:
            includeFieldDistortion = self["Camera/IncludeFieldDistortion"]
            mappedDistortion       = False
            inverseDistortionCoefficients = self["Camera/FieldDistortion/ConstantInverseCoefficients"]
            pathToPsfFile          = None

        solarPanelOrientation = np.deg2rad(float(self["Platform/SolarPanelOrientation"]))

        # Convert the pixel coordinates to focal plane coordinates [mm]

        xFPmm, yFPmm = rf.pixelToFocalPlaneCoordinates(cols, rows, pixelSize,
                                                       ccdZeroPointX, ccdZeroPointY, CCDangle)

        # If distortion is required in the yaml input file, distort the focal plane coordinates [mm]
        if mappedDistortion:
            for i in range(len(xFPmm)):
                xFPmm[i], yFPmm[i] = rf.mappedDistortedToUndistortedFocalPlaneCoordinates(xFPmm[i], yFPmm[i], pathToPsfFile)

        elif (includeFieldDistortion == "yes" or
              includeFieldDistortion == "1"   or
              includeFieldDistortion == True):
            xFPmm, yFPmm = rf.distortedToUndistortedFocalPlaneCoordinates(xFPmm, yFPmm, inverseDistortionCoefficients, focalLength)

        # Convert the focal plane coordinates to equatorial sky coordinates [rad]

        ra, dec = rf.focalPlaneToSkyCoordinates(xFPmm, yFPmm, raPlatform, decPlatform, solarPanelOrientation, tiltAngle, azimuthAngle, focalPlaneAngle, focalLength)

        # Convert sky coordinates to degrees

        ra = np.rad2deg(ra)
        dec = np.rad2deg(dec)

        # Save the sky coordinates (in [deg]) to the star catalog file

        myFile = open(starCatalogPath, "w")
        myFile.write("# RA DEC Vmag starID\n")
        for n in range(len(ra)):
            myFile.write("{0}  {1}  {2}  {3}\n".format(ra[n], dec[n], magnitudes[n], starIDs[n]))
        myFile.close()

        # Set the "ObservingParameters/StarCatalogFile" tag in the yaml tree

        self["ObservingParameters/StarCatalogFile"] = starCatalogPath





    def createPhotometryTargetFile(self, starIDs, fileName):

        """Create a photometry file list in ascii format and sets it to the YAML input.

        Parameters
        ----------
        starIDs : ndarray
            Array with IDs of the star (integers)
        fileName : str
            Path of the photometry file that will be written.

        Return
        ------
        A file will be saved, containing the star IDs that photometry should be performed on.
        The "Photometry/TargetFileName" tag in the yaml tree will be changed to the given 
        starCatalogPath
        """

        # Check if starIDs are just a number

        if type(starIDs) is int:
            starIDs = [starIDs]

        # Create photometry list file

        np.savetxt(fileName, np.transpose(starIDs), delimiter=" ", fmt="%d")

        # Set this to simulation

        self["Photometry/TargetFileName"] = fileName






    def createDriftFile(self, quarter, fileName, model="poly", plot=False):

        """Create a photometry file list in ascii format and sets it to the YAML input.

        Parameters
        ----------
        starIDs : ndarray
            Array with IDs of the star (integers)
        fileName : str
            Path of the photometry file that will be written.

        Return
        ------
        A file will be saved, containing the star IDs that photometry should be performed on.
        The "Telescope/UseDriftFromFile" tag in the yaml tree will be changed to the given
        DriftFileName.
        """

        # Create TED file

        it.getTED(quarter=quarter, model=model, outfile=fileName, plot=plot)

        # Set this to simulation

        self["Telescope/UseDrift"]         = True
        self["Telescope/UseDriftFromFile"] = True
        self["Telescope/DriftFileName"]    = fileName





    def getReadoutTime(self):

        """Fetch the readout time.
        
        This function determine the duration of:
          - the readout that takes place before the next exposure starts,
          - and the readout that takes place during the next exposure,
            depending on the camera type (normal / fast) and the readout mode
            (nominal / partial readout).

        For the normal cameras the entire CCD is read out (with open shutter) after
        the exposure during a time interval called 'readoutTimeBeforeNextExposure'.
        Only after this readout, a new exposure is started.
        
        For the fast camera, half of the CCD is first quickly frame-transferred,
        after which it is read out slowly. In this case a new exposure is already
        started after the quick frame-transfer, and starts thus during the slow readout
        of the previous exposure. Hence the need for two parameters: 
        'readoutTimeBeforeNextExposure' and 'readoutTimeDuringNextExposure'.

        Return
        ------
        readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure
        """

        isFastCamera = self["Telescope/GroupID"] == "Fast"
        ccdPosition  = self["CCD/Position"]

        if ccdPosition == "Custom":
            numRows         = self["CCD/NumRows"]         # [pixels]
            numColumns      = self["CCD/NumColumns"]      # [pixels]
            firstRowExposed = self["CCD/FirstRowExposed"] # [pixels]
        else:
            # Positions start with 1, while the index starts at 0
            idx = ccdPosition - 1

            numRows    = self["CCDPositions/NumRows"][idx]     # [pixels]
            numColumns = self["CCDPositions/NumColumns"][idx]  # [pixels]

            if isFastCamera:
                firstRowExposed = self["CCDPositions/FirstRowForFastCamera"][idx]    # [pixels]
            else:
                firstRowExposed = self["CCDPositions/FirstRowForNormalCamera"][idx]  # [pixels]

        readoutMode = self["CCD/ReadoutMode/ReadoutMode"]

        if (readoutMode != "Nominal") and (readoutMode != "Partial"):
            raise ValueError("Simulation::getReadoutTime() Unknown readout mode " +
                             f"specification in configuration file: {readoutMode}")


        serialTransferTime       = self["CCD/SerialTransferTime"]       * 1e-9  # [ns] -> [s]
        parallelTransferTime     = self["CCD/ParallelTransferTime"]     * 1E-6  # [micro s] -> [s]
        parallelTransferTimeFast = self["CCD/ParallelTransferTimeFast"] * 1E-6  # [micro s] -> [s]

        numColumnsBiasMap =  self["SubField/NumBiasPrescanColumns"]    # [pixels]
        numRowsSmearingMap = self["SubField/NumSmearingOverscanRows"]  # [pixels]

        # Both detector halves are read out simultaneously
        # -> columns read out by the FEE:
        # 		- half of the CCD
        # 		- serial pre-scan
        # 		- (serial over-scan)

        numColumnsReadout = numColumns / 2 + numColumnsBiasMap # + numRowsSerialOverScan

        # How many rows will be actually read out by the FEE?
        # 	- nominal mode: image area + parallel over-scan
        #      normal camera: image area = whole CCD
        #      fast camera: image area = lower half of the CCD
        #	- partial readout: configurable
        # The rest of the image area will be dumped

        #--- Fast camera

        if isFastCamera:
            
            # Move the upper half of the CCD down to the lower half, row-by-row

            numRowsFrameTransfer          = numRows - firstRowExposed
            readoutTimeBeforeNextExposure = numRowsFrameTransfer * parallelTransferTimeFast

            # The actual readout of the lower half of the CCD (after frame transfer) is done
            # while the next exposure has already started

            # Nominal mode

            if readoutMode == "Nominal":
                numRowsReadout = firstRowExposed + numRowsSmearingMap
                numRowsDump = 0

            # Rows read out by the FEE: rows in the block (other rows in image area are dumped)
            # Note: no parallel over-scan

            elif readoutMode == "Partial":
                numRowsReadout = self["CCD/ReadoutMode/Partial/NumRowsReadout"]
                numRowsDump    = firstRowExposed - numRowsReadout

            readoutTimeDuringNextExposure = (numRowsDump * parallelTransferTimeFast +
                                             numRowsReadout * (parallelTransferTime +
                                                               numColumnsReadout *
                                                               serialTransferTime))

        #--- Normal camera

        else:
            # Nominal mode (full-frame readout)

            if readoutMode == "Nominal":
                
                # Rows read out by the FEE:
                #  - rows of image area
                #  - parallel over-scan

                numRowsReadout = numRows + numRowsSmearingMap

                # No rows dumped

                numRowsDump = 0;

            # Partial readout

            elif readoutMode == "Partial":
                
                # Rows read out by the FEE: rows in the block (other rows in image area are dumped)
                # Note: no parallel over-scan

                numRowsReadout = self["CCD/ReadoutMode/Partial/NumRowsReadout"]
                numRowsDump = numRows - numRowsReadout


            readoutTimeBeforeNextExposure = (numRowsDump * parallelTransferTimeFast +
                                             numRowsReadout * (numColumnsReadout *
                                                               serialTransferTime +
                                                               parallelTransferTime))
            readoutTimeDuringNextExposure = 0

        return readoutTimeBeforeNextExposure, readoutTimeDuringNextExposure





    def getStarsWithinCameraGroup(self, raPF, decPF, ra, dec, kappa=-8,
                                  camGroup=1, quarter=1, sizeSubfield=6):

        """Fetch all stars within a camera group.

        This function determines if any star from a catalogue is within the FOV of a specific
        PLATO camera group.

        TODO: Could be implemented in a better and faster manner!

        Parameters
        ----------
        camGroup : int [1, 2, 3, 4]
            N-CAM camera group ID
        raPF : float
            Right acsension of pointing field [deg]
        decPF : float
            Declination of pointing field [deg]
        ra : list, array
            Right ascension of stars to be checked against [deg]
        dec : list, tuple, array
            Declination of stars to be checked against [deg]
        kappa : float
            Rotation of the platform [deg]

        Return
        ------
        dex : list
           Indices booleans with True being all stars within camera group FOV
        distanceOA : list
           Distance of each from the camera's optical axis [mm]
        """

        # Telescope config

        raPlatformDeg  = self["ObservingParameters/RApointing"]  = raPF   # [deg]
        decPlatformDeg = self["ObservingParameters/DecPointing"] = decPF  # [deg]

        raPlatformRad  = np.deg2rad(raPlatformDeg)   # [rad]
        decPlatformRad = np.deg2rad(decPlatformDeg)  # [rad]

        focalLength      = float(self["Camera/FocalLength/ConstantValue"]) * 1000.0  # [m] -> [mm]
        focalPlaneAngle  = np.deg2rad(float(self["Camera/FocalPlaneOrientation/ConstantValue"]))

        solarPanelOrientation = self["Platform/SolarPanelOrientation"] = math.fmod(quarter * 90. - kappa, 360.)
        solarPanelOrientation = np.deg2rad(float(solarPanelOrientation))

        raTargetsRad  = np.deg2rad(ra)   # [rad]
        decTargetsRad = np.deg2rad(dec)  # [rad]

        # Loop over each star for this cam-group

        self["Telescope/GroupID"] = camGroup
        azimuthTelescope = np.deg2rad(self["CameraGroups/AzimuthAngle"][camGroup-1])
        tiltTelescope    = np.deg2rad(self["CameraGroups/TiltAngle"][camGroup-1])

        dexGroup   = np.zeros(len(ra), dtype=bool)
        distanceOA = np.zeros(len(ra))

        for i in range(len(ra)):

            subfieldIsOnCCD = self.setSubfieldAroundCoordinates(raTargetsRad[i], decTargetsRad[i],
                                                                sizeSubfield, sizeSubfield,
                                                                normal=True)
            if subfieldIsOnCCD:

                xFP, yFP = rf.skyToFocalPlaneCoordinates(raTargetsRad[i], decTargetsRad[i],
                                                         raPlatformRad, decPlatformRad,
                                                         solarPanelOrientation,
                                                         tiltTelescope, azimuthTelescope,
                                                         focalPlaneAngle, focalLength)

                distanceOA[i] = np.rad2deg(rf.gnomonicRadialDistanceFromOpticalAxis(xFP, yFP,
                                                                                    focalLength))

                if distanceOA[i] < self['CCD/RelativeTransmissivity/RadiusFOV']:
                    dexGroup[i] = True
                else:
                    dexGroup[i] = False

        # Return parameters

        return dexGroup, distanceOA





    def getSubPixelPositions(self, raPF, decPF, ra, dec, camGroup=1, quarter=1):

        """Fetch the subpixel postion for all stars.

        This function determines the subpixel postion of all stars for the first
        exposure in the simulation. The function also returns the CCD ID of each
        star in the focal plane.

        Parameters
        ----------
        raPF : float
            Right acsension of pointing field [deg]
        decPF : float
            Declination of pointing field [deg]
        ra : list, array
            Right ascension of stars to be checked against [deg]
        dec : list, tuple, array
            Declination of stars to be checked against [deg]

        Return
        ------
        ccdCode : int
           CCD ID of belonging to the star
        xCCD : float32
           Subpixel position of star in x (row)
        yCCD : float32
           Subpixel position of star in y (column)
        """

        # Telescope config

        raPlatformDeg  = self["ObservingParameters/RApointing"]  = raPF   # [deg]
        decPlatformDeg = self["ObservingParameters/DecPointing"] = decPF  # [deg]

        raPlatformRad  = np.deg2rad(raPlatformDeg)   # [rad]
        decPlatformRad = np.deg2rad(decPlatformDeg)  # [rad]

        focalLength      = float(self["Camera/FocalLength/ConstantValue"]) * 1000.0  # [m] -> [mm]
        focalPlaneAngle  = np.deg2rad(float(self["Camera/FocalPlaneOrientation/ConstantValue"]))

        solarPanelOrientation = self["Platform/SolarPanelOrientation"] = math.fmod(quarter*90.,360.)
        solarPanelOrientationRad = np.deg2rad(float(solarPanelOrientation))

        raTargetsRad  = np.deg2rad(ra)   # [rad]
        decTargetsRad = np.deg2rad(dec)  # [rad]

        # Loop over each star for this cam-group

        self["Telescope/GroupID"] = camGroup
        azimuthTelescopeRad = np.deg2rad(self["CameraGroups/AzimuthAngle"][camGroup-1])
        tiltTelescopeRad    = np.deg2rad(self["CameraGroups/TiltAngle"][camGroup-1])

        # CCD properties

        pixelSize = float(self["CCD/PixelSize"])

        ccdCode = np.zeros(len(ra))
        xCCD    = np.zeros(len(ra))
        yCCD    = np.zeros(len(ra))

        for i in range(len(ra)):

            subfieldIsOnCCD = self.setSubfieldAroundCoordinates(raTargetsRad[i], decTargetsRad[i],
                                                                6, 6, normal=True)
            if subfieldIsOnCCD:

                # Fetch CCD code and pixel coordinates (account for field distortion in included)

                if self["Camera/IncludeFieldDistortion"]:
                    includeFieldDistortion = self["Camera/IncludeFieldDistortion"]
                    if self["Camera/FieldDistortion/Type"] == 'FromFile':
                        mappedDistortion = True
                        distortionCoefficients = self["Camera/FieldDistortion/CoefficientsFromFile"]
                    else:
                        mappedDistortion = False
                        distortionCoefficients = self["Camera/FieldDistortion/ConstantCoefficients"]
                else:
                    includeFieldDistortion = False
                    mappedDistortion = False
                    distortionCoefficients = False

                out = rf.getCCDandPixelCoordinates(raTargetsRad[i],
                                                   decTargetsRad[i],
                                                   raPlatformRad,
                                                   decPlatformRad,
                                                   solarPanelOrientationRad,
                                                   tiltTelescopeRad,
                                                   azimuthTelescopeRad,
                                                   focalPlaneAngle,
                                                   focalLength,
                                                   pixelSize,
                                                   includeFieldDistortion,
                                                   normal=True,
                                                   mappedDistortion=mappedDistortion,
                                                   distortionCoefficients=distortionCoefficients)
                ccdCode[i], xCCD[i], yCCD[i] = out[0], out[1], out[2]
                
        # That's it!
        
        return ccdCode, xCCD, yCCD
