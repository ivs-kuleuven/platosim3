#!/usr/bin/env python3

"""
This script uses the PIC targets and their contaminants to simulate realistic imagettes
or lightcurves using PlatoSim. It can simulate any of the 24 normal cameras/telescopes,
for an arbitary number of mission quarters. Since a nescessary rotation (along the roll
axis) of the spacecraft platform is required in order to repoint the solar panels every
90 days, simulations cannot realistically extend beyond a quarter of a year. Given an
PIC input sample, the script will automatically simulation all targets being visible by
any camera falling on one of the 4 CCDs. The fact that the PLATO spacecraft will be
equipped with 4 camera groups consisting of 6 cameras each, constrains the efficient use
of nodes and cores on the HPC. In order to make the simulation as realistic as possible,
random seats of various intrinsic- and instrumental effects are included, meaning that
each camera within each camera group needs to be simulated independently (which indeed
are the raw imaging output of PLATO).
"""

# Built-in
import os
import sys
import glob
import math
import shutil
import argparse
import tracemalloc

# External packages
import h5py
import datetime
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d, CubicSpline
from pathlib import Path

# PlatoSim
import platosim.utilities       as ut
import platosim.referenceFrames as rf
from platosim.simulation   import Simulation
from platosim.simfile      import SimFile
from platosim.utilities    import errorcode, pdAddColumn, getPointingField
from platosim.plot         import drawStarInCCDfocalPlane, plotSubfieldAnimation
from platosim.matplotlibrc import setup
setup()

#==============================================================#
#                        PLATOnium CLASS                       #
#==============================================================#

class PLATOnium(object):

    """Class for running multi-camera and multi-quarter PlatoSim simulations.
    """
    
    def __init__(self, args):

        # PARSED ARGUMENTS
        
        self.targetNo     = args.starID
        self.group        = args.groupID
        self.camera       = args.cameraID
        self.quarter      = args.quarter

        self.inputFile     = args.infile
        self.outputDir     = args.outdir
        self.hpcDir        = args.hpcdir
        self.project       = args.project
        self.inputFileName = args.yaml
        self.simPrefix     = args.prefix
        self.starcatFile   = args.starcat
        self.varSourceFile = args.varfile
        self.varSourceList = args.varlist
        self.compress      = args.compress
        self.statistics    = args.statistics
        
        self.cadence      = args.cadence
        self.simTime      = args.tdur
        self.simExposures = args.nexp
        self.simBeginExp  = args.bexp
        self.picID        = args.pic
        self.maskUpdate   = args.mask
        self.seed         = args.seed
        self.mag          = args.mag
        self.performance  = args.performance
        
        self.noCon        = args.nocon
        self.reuseJitter  = args.reusejitter
        self.fullFrame    = args.fullframe
        
        self.sample          = args.sample
        self.conDeltaMag     = args.con_dmag
        self.conFluxError    = args.con_ferr
        self.tarFluxError    = args.tar_ferr
        self.tarAbsCenError  = args.tar_cerr
        self.jitterDriftOff  = args.jit_off

        # Make sure that parsing of mandatory arguments are physical

        # Normal cameras
        if self.group in [1, 2, 3, 4]:
            if self.camera in [1, 2, 3, 4, 5, 6]:
                self.groupID = self.group
            else:
                errorcode('error', 'Camera can only be [1, 2, 3, 4, 5, 6]')
        # Fast cameras
        elif self.group == 5:
            self.groupID = 'Fast'
            if self.camera == 1:
                self.cameraID = 'blue'
            elif self.camera == 2:
                self.cameraID = 'red'
            else:
                errorcode('error', 'Fast camera can only be [1, 2] = [blue, red]')
        else:
            errorcode('error', 'Camera-group can only be [1, 2, 3, 4, 5] (Fast = 5)')

        # Select full-frame CCD
        if self.fullFrame:
            self.ccdCode = self.targetNo
            if not self.ccdCode in [1, 2, 3, 4]:
                errorcode('error', 'CCD code can only be [1, 2, 3, 4]')
        
        # VERBOSITY (a.k.a log level) -> Identical to PlatoSim usage
        # verbose = 0: Cluster mode: Disabling print and warnings, and no log files are saved
        # verbose = 1: Default mode: Print details to bash but do not save log files
        # verbose = 3: Debug mode  : Print details to bash and saves all log files        
        if args.verbose == 0 or self.statistics:
            self.verbose = 0
            self.verbose_platosim = 0
            # Bash extention to write no output for the L1 pipeline
            self.devnull = '> /dev/null'
            # Turn off warnings
            import warnings
            warnings.filterwarnings("ignore")
        elif args.verbose is None or args.verbose == 1:
            self.verbose = 1
            self.verbose_platosim = 0
            self.devnull = ''
            # Turn off warnings
            import warnings
            warnings.filterwarnings("ignore")
        else:
            self.verbose = 3
            self.verbose_platosim = 3
            self.devnull = ''

        # Overwrite simulation
        if args.overwrite or self.statistics:
            self.overwrite = True
        else:
            self.overwrite = False

        # Save animation
        if args.animation:
            self.animation = True
        else:
            self.animation = False

        # Start software writing
        if self.verbose > 0:
            errorcode('software', '\nPLATOnium')

            
        # PATHS AND INPUT
        
        # Absolute pwd path
        self.path = Path(__file__).parent.resolve()

        # PlatoSim inputfiles
        self.platosimInputDir = self.path.joinpath(os.getenv('PLATO_PROJECT_HOME'), 'inputfiles')

        # Sort the workdir input paths and files
        if self.project and self.inputFileName:
            self.simDir    = self.path.joinpath(os.getenv('PLATO_WORKDIR'), self.project)
            self.inputDir  = self.simDir / 'input'
            self.inputFile = self.inputDir / self.inputFileName
        elif self.project:
            self.simDir    = self.path.joinpath(os.getenv('PLATO_WORKDIR'), self.project)
            self.inputDir  = self.simDir / 'input'
            self.inputFile = self.inputDir / 'inputfile.yaml'
        elif self.inputFile:
            self.inputFile = Path(self.inputFile).resolve()
            self.inputDir  = self.inputFile.parents[0]
            self.simDir    = self.inputFile.parents[1]
        else:
            errorcode('error', 'Either -i, --project or --yaml is needed to locate the inputfile!')
        
        # Check if the inputfile.yaml exists
        if not self.inputFile.is_file():
            errorcode('error', f'File inputfile.yaml do not exist! Alternamtively use "-i" or "--yaml"')        
            
        # L1 pipeline paths
        if self.sample:
            self.platoBin = self.path.joinpath(os.getenv('PLATO'), 'bin')
            self.platoLib = Path(os.path.split(sys.executable)[0]).resolve()

            
        # PHYSICAL PARAMETERS

        # Inclusion thresholds for contaminants [delta mag]
        if not self.conDeltaMag: self.conDeltaMag = 5
        
        # Defualt L1 pipeline parameters
        self.bsres           = 10   # [subpixel]
        self.prnuError       = 0.1  # [%]
        self.maskUpdateThres = 0.0  # [pixel]
        if not self.conFluxError:   self.conFluxError   = 10    # Flux error of contaminat stars [%]
        if not self.tarFluxError:   self.tarFluxError   = 1     # Flux error of target star [%]
        if not self.tarAbsCenError: self.tarAbsCenError = 0.02  # Absolute error of target centroid position [pixel]

        # Monitor script speed
        self.tic  = datetime.datetime.now()
        self.tic0 = datetime.datetime.now()

        


        
    def load_stars(self):
        """
        Module to load the stellar targets and contaminants.
        """

        if (self.verbose == 3) or (self.fullFrame and self.verbose > 0):
            print('\nLoading stellar catalogue..')

            
        # FULL FRAME
            
        if self.fullFrame:

            # Load stellar catalogue
            starcat = Path(glob.glob(f'{str(self.inputDir)}/starcat**group{self.group}.ftr')[0])
            if not starcat.is_file():
                errorcode('error', 'No star catalogue found in the project "input/" directory!')
            self.dx = pd.read_feather(starcat)
            
            # Save star catalogue
            self.ds = pd.DataFrame()
            self.ds['ra']  = self.dx.ra
            self.ds['dec'] = self.dx.dec
            self.ds['mag'] = self.dx.mag
            self.ds['ids'] = np.arange(0, len(self.ds.ra)).astype(int)            

            return


        # SUBFIELD

        # Fetch stars either from custum catalogue are the default PIC setup
        
        if self.starcatFile is not None:

            df = pd.read_csv(self.starcatFile, sep=' ',
                             names=['PIC', 'ra', 'dec', 'mag', 'dis'])
            
            # Change IDs all to be the target
            df.PIC = np.ones(len(df))
            df.PIC = df.PIC.astype('int')
            
            # Define data frames
            self.df = df.loc[0]
            dc = df.iloc[1:]

        else:

            # Add sample name if pipeline is activated
            if self.sample is not None:
                extra_str = f'_{self.sample}'
            else:
                extra_str = ''

            # Fetch PIC targets and contaminants
            try:
                picTarFile = glob.glob(str(self.inputDir) + f'/starcat{extra_str}**targets.ftr')[0]
                picConFile = glob.glob(str(self.inputDir) + f'/starcat{extra_str}**contaminants.ftr')[0]
                df = pd.read_feather(picTarFile)
                dc = pd.read_feather(picConFile)
            except IndexError:
                errorcode('error', f'Stellar catalogue for {self.sample} sample do not exist!')

            # Merge for full frame
            self.dx = pd.concat([df, dc])
            
            # Correct indicing and allow a specific star to be choosen
            if self.targetNo == 0:
                errorcode('error', 'Star ID indicing starts from 1 and not 0!')
            elif self.picID is not None:
                try:
                    self.targetNo = np.where(df['PIC'] == self.picID)[0][0]
                except IndexError:
                    errorcode('error', f'PIC {self.picID} star does not exist in catalogue:' +
                              f'\n{picTarFile}')
            else:
                self.targetNo -= 1

            # Select target star
            self.df = df.iloc[self.targetNo]

            
        # Additional info for subfield simulations
            
        if not self.fullFrame:
            
            # If requested select only the target or including contaminants
            if self.noCon:
                self.numCon = 0
                self.dc = dc[dc['PIC'] == self.numCon]
            else:
                self.dc = dc[dc['PIC'] == self.df['PIC']]
                self.numCon = len(self.dc['PIC'])
                # Sort contaminants after their distance
                self.dc = self.dc.sort_values(by=['dis'])

            # If requested overwrite magnitude of target star
            if self.mag:
                self.df.mag = self.mag
            
            # Save star catalogue
            self.ds = pd.DataFrame()
            self.ds['ra']  = np.append(self.df['ra'],  self.dc['ra'])
            self.ds['dec'] = np.append(self.df['dec'], self.dc['dec'])
            self.ds['mag'] = np.append(self.df['mag'], self.dc['mag'])
            self.ds['ids'] = np.arange(1, self.numCon+2)




            
    def configure_output(self):

        """Module to create, configure, and select the correct output folders.
        """

        # Final destination on the cluster
        if self.hpcDir:
            self.hpcDir = Path(self.hpcDir).resolve()
            # Create directory
            if not self.hpcDir.is_dir():
                self.hpcDir.mkdir(parents=True, exist_ok=True)
                os.system(f'chmod 755 {self.hpcDir}')
            
        # Default output path
        if self.outputDir is None:
            self.outputDir = self.simDir.joinpath('output')
        else:
            self.outputDir = Path(self.outputDir).resolve()
        # Create directory
        if not self.outputDir.is_dir():
            self.outputDir.mkdir(parents=True, exist_ok=True)
            os.system(f'chmod 755 {self.outputDir}')
            
        # Custom prefix
        if self.simPrefix is not None:
            self.simPrefix = f'{self.simPrefix}_'
        else:
            self.simPrefix = ''
        
        # Set general output filename
        self.outputFileName = f'{self.targetNo + 1}'.zfill(9)
            
        # Select suffix of observation
        if self.groupID == 'Fast':
            self.obsPrefix = f'Fcam{self.camera}_Q{self.quarter}'
        else:
            self.obsPrefix = f'Ncam{self.group}.{self.camera}_Q{self.quarter}'

        # Combine file name
        
        if self.fullFrame:
            # Full frame mode
            self.outputFileName = f'{self.obsPrefix}_ccd{self.ccdCode}'
        elif self.sample is None:
            # Standard mode
            self.outputFileName = f'{self.simPrefix}{self.outputFileName}_{self.obsPrefix}'
        else:
            # Pipline mode
            self.starID = self.outputFileName
            relativeDirStar = f'{self.sample}/Q{self.quarter}/Ncam{self.group}{self.camera}'

            # L1 pipeline paths
            self.pipelineDir = f'{self.outputDir}/{relativeDirStar}'
            self.outputDirStarIDsim = f'{self.pipelineDir}/{self.starID}'
            self.outputDirStarIDnew = f'{self.outputDir}/reduced/{self.sample}/{self.starID}'
                        
            # Microscan paths
            self.microscanDir = f'{self.outputDir}/microscan/{relativeDirStar}/{self.starID}'
            self.microscanDirStarID = f'{self.microscanDir}/{self.starID}'
            self.microscanDirInvers = f'{self.microscanDir}/inversion'
            
            # Create paths
            try: os.makedirs(self.outputDirStarIDsim)
            except: pass
            try : os.makedirs(self.outputDirStarIDnew)
            except: pass
            try: os.makedirs(self.microscanDirStarID)
            except: pass
            try: os.makedirs(self.microscanDirInvers)
            except: pass

            # Normal output path
            self.outputDir = Path(self.outputDirStarIDsim).resolve()

            # Path + Prefix of simulations
            self.microscanSimName = f'{self.microscanDirStarID}/{self.outputFileName}'

        # Store final file name
        self.outputSimName = f'{self.outputDir}/{self.outputFileName}'

        


        

    def track_statistics(self, flag):

        with open(self.inputDir / 'statistics.txt', 'a+') as f:
            # Move read cursor to the start of file.
            f.seek(0)
            # If file is not empty then append '\n'
            if len(f.read(100)) > 0: f.write("\n")
            # Append text at the end of file
            f.write(f'{int(self.targetNo+1)},{self.group},{self.camera},{self.quarter},{flag}')



            

    def init_sim(self):

        """Module to initialize the the PlatoSim simulation object.
        """
        
        # INITIALIZE SIMULATION

        # Print to bash
        if self.verbose > 0:
            errorcode('module', '\nInitialize and configure PlatoSim\n')

        # Setting up a test simulation environement
        sim = Simulation(self.outputFileName, self.inputFile)

        # Start time of simulation
        timeQuarter = ut.year()/86400/4  # [days]
        self.timeStart = round(timeQuarter * (self.quarter - 1) * 86400.)

        
        # CONFIGURE CAMERA

        # NOTE these function sets the correct CCD configuration and cadence
        #      and if requested also performance and time conditions
        
        if self.groupID == 'Fast':
            # Parameter "normal" used in subfield selection
            normal = False
            sim.useFastCamera(self.cameraID, self.performance, self.timeStart)
        else:
            normal = True
            sim.useNormalCamera(self.performance, self.timeStart)
            
        
        # CONFIGURE TIMING

        # NOTE: CCD offset is automatically set by setSubfieldAroundCoordinates()
        
        # Cadence of time series [s]
        if self.cadence:
            sim['ObservingParameters/CycleTime'] = self.cadence
        else:
            self.cadence = sim['ObservingParameters/CycleTime']
            
        # Check of begin exposure number is parsed
        if not self.simBeginExp:
            self.simBeginExp = 0
        
        # Apply start time relative mission BOL
        self.beginExposureNr = round(self.timeStart / self.cadence) + self.simBeginExp
        sim['ObservingParameters/BeginExposureNr'] = self.beginExposureNr

        # Duration of time series
        if self.simExposures:
            # Setting timeseries by N exposures given by user
            self.numExposures = self.simExposures
        elif self.simTime is not None:
            # Setting timeseries by time given by user
            self.numExposures = round(self.simTime * 86400. / self.cadence)
        else:
            # Setting time series to full quarter
            # NOTE Two days are lost due to platform roll, thermal stabilisation,
            # downlink, microscanning, etc.
            # TODO make configurable/flag option to change inter-quarter downtime
            # (normal distributed?)
            self.numExposures = round((timeQuarter - 2.) * 86400. / self.cadence)
            

        # PHOTOMETRY ALA MARCHIORI

        # The mask-update interval [days]
        # The mask-update rate [event after N number of exposures]
        if self.maskUpdate:
            sim['Photometry/MaskUpdateInterval'] = self.maskUpdate
            self.maskUpdateRate = round(self.maskUpdate * 86400. / self.cadence)
        else:
            self.maskUpdateRate = round(sim['Photometry/MaskUpdateInterval'] * 86400./self.cadence)

            
        # CONFIGURE PAYLOAD

        if sim["Platform/Orientation/Source"] == "Quaternion":
            print('Using Quaternion to setup platform!')
        
        # Select PLATO pointing field from inputfile
        self.pointingField = sim['ObservingParameters/StarCatalogFile']
        alpha, delta, kappa = getPointingField(self.pointingField)
        sim["Platform/Orientation/Angles/RAPointing"]  = alpha
        sim["Platform/Orientation/Angles/DecPointing"] = delta

        # Solar panel orientation [deg]: Q(N) = {0, 90, 180, 270} + kappa
        solarPanelOrientationDeg = ut.getSolarPanelOrientation(kappa, self.quarter)
        sim["Platform/Orientation/Angles/SolarPanelOrientation"] = solarPanelOrientationDeg

        # Set the Camera-group ID, Alt (tilt) [deg], and Az [deg]
        sim["Telescope/GroupID"]      = 'Custom'
        sim["Telescope/TiltAngle"]    = sim["CameraGroups/TiltAngle"][self.group-1]
        sim["Telescope/AzimuthAngle"] = sim["CameraGroups/AzimuthAngle"][self.group-1]


        # POINTING ERRORS
        
        # Include spacecraft Pointing Repeatability Error (PRE) between consecutive quarters
        # NOTE: Included if the file "instrumentPRE.txt" is available in the input folder
        inputFilePRE = self.inputDir.joinpath('instrumentPRE.txt')
        if inputFilePRE.is_file():
            PRE = np.loadtxt(inputFilePRE)
            # Catch one dimentional arrays
            try: PRE.shape[1]
            except: PRE = np.array([PRE])
            # Apply pointing errors
            dex = np.where(PRE[:,0] == self.quarter)[0]
            try: dex[0]
            except:
                errorcode('warning', 'Cannot apply pointing error model: ' +
                          'no matching quarters in instrumentPRE.txt!')
            else:
                if self.verbose > 0:
                    print('Applying pointing error      (PRE FromFile)')
                sim["Platform/Orientation/Angles/RAPointing"]            += PRE[dex, 1][0]
                sim["Platform/Orientation/Angles/DecPointing"]           += PRE[dex, 2][0]
                sim["Platform/Orientation/Angles/SolarPanelOrientation"] += PRE[dex, 3][0]

        # Absolute Pointing Error (APE) due to camera misalignments
        # NOTE: Included if "instrumentAPE.txt" is available in the input
        inputFileAPE = self.inputDir.joinpath('instrumentAPE.txt')
        if inputFileAPE.is_file():
            if self.verbose > 0 :
                print('Applying camera misalignment (APE FromFile)')
            APE = np.loadtxt(inputFileAPE)
            dex = (self.group - 1) * 6 + self.camera - 1
            sim["Telescope/TiltAngle"]    += APE[dex, 0]
            sim["Telescope/AzimuthAngle"] += APE[dex, 1]

        # Thermo-Elastic Distortion (TED)
        # The camera(s) drift due to the thermal gradient of
        # the interface between the camera and the optical bench.
        # NOTE: Included if "instrumentTED.txt" is available in input
        inputFileTED = self.inputDir.joinpath('instrumentTED.txt')
        if inputFileTED.is_file():
            sim["Telescope/UseDrift"]      = True
            sim["Telescope/DriftSource"]   = 'FromFile'
            sim["Telescope/DriftFileName"] = inputFileTED
        if sim["Telescope/UseDrift"] and self.verbose > 0:
            if sim["Telescope/DriftSource"] == 'FromFile':
                source = 'FromFile'
            else:
                source = 'RedNoise'
            print(f'Applying camera drift        (TED {source})')

        # Attitude Orbit Control System (AOCS) jitter
        # First check if file from payload is present
        # NOTE: Included if "instrumentACS.txt" is available in input
        inputFileAOCS = self.inputDir.joinpath('instrumentACS.txt')
        if inputFileAOCS.is_file():
            sim["Platform/UseJitter"]      = True
            sim["Platform/JitterSource"]   = 'FromFile'
            sim["Platform/JitterFileName"] = inputFileAOCS
            
        # Check if "AOCS_Q<quarterNo>.txt" is present to be reused for all quarters
        # NOTE: Not recommended but used for PLATO-KUL-PL-TN-0023
        elif self.reuseJitter:
            sim["Platform/UseJitter"]    = True
            sim["Platform/JitterSource"] = 'FromFile'
            inputFileAOCS   = self.inputDir.joinpath(sim["Platform/JitterFileName"])
            inputFileAOCS_Q = f'{self.inputDir}/AOCS_Q{self.quarter}.txt'
            # Check if the exists or else create new time column
            if Path(inputFileAOCS_Q).is_file():
                pass
            elif self.quarter == 1:
                os.system(f'cp {inputFileAOCS} {inputFileAOCS_Q}')
            else:
                # Load file given in YAML input
                data = np.loadtxt(inputFileAOCS)
                t, x, y, z = data[:,0], data[:,1], data[:,2], data[:,3]
                # Generate new time column
                cadence = np.diff(t)[0]  # [s -> 8 Hz]
                t = np.arange(len(t)) * cadence + self.timeStart
                # Save data to input folder
                np.savetxt(inputFileAOCS_Q, np.transpose([t,x,y,z]),
                           fmt=['%.3f', '%.9f', '%.9f', '%.9f'])
            # Set filepath to new file
            sim["Platform/JitterFileName"] = inputFileAOCS_Q
        # Print to bash that jitter is included
        if sim["Platform/UseJitter"] and self.verbose > 0:
            if sim["Platform/JitterSource"] == 'FromFile':
                source = 'FromFile'
            else:
                source = 'RedNoise'
            print(f'Applying platform jitter     (ACS {source})')

        # Thermal transients from data gaps
        # NOTE: Included if "instrumentCCD.txt" is available in input
        inputFileCCD = self.inputDir.joinpath('instrumentCCD.txt')
        if inputFileCCD.is_file():
            sim["CCD/Temperature"]         = "FromFile"
            sim["CCD/TemperatureFileName"] = inputFileCCD
            print(f'Applying gain transients     (CCD FromFile)')

            
        # FULL-FRAME SIMULATION

        if self.fullFrame:

            # Turn off photometry
            sim['Photometry/IncludePhotometry'] = False

            # Set CCD parameters
            self.isOnCCD         = True
            sim["CCD/Position"]  = str(self.ccdCode)

            if self.groupID == 'Fast':
                shieldRows = sim["CCDPositions/MetallicShield/ShieldRowCoordinates"]
                shieldCols = sim["CCDPositions/MetallicShield/ShieldColumnCoordinates"]
                sim["SubField/ZeroPointRow"]    = shieldRows[0]
                sim["SubField/ZeroPointColumn"] = shieldCols[0]
                sim["SubField/NumRows"]         = shieldRows[1] - shieldRows[0]
                sim["SubField/NumColumns"]      = shieldCols[1] - shieldCols[0]
            else:
                sim["SubField/NumRows"]         = sim["CCDPositions/NumRows"][0]
                sim["SubField/NumColumns"]      = sim["CCDPositions/NumColumns"][0]

            # Control output requirements
            sim["ControlHDF5Content/GroupByExposure"]    = True
            sim["ControlHDF5Content/WritePixelMaps"]     = True
            sim["ControlHDF5Content/WriteStarPositions"] = True

            return sim
        
        # SUBFIELD SIMULATION
        
        # Try to set a subfield around the coordinates on one of the 4 CCDs of the telescope.
        # This will fail (return = False) if the subfield is not visible by any of the 4 CCDs
        # or if the subfield is too large to entirely fit on a CCD. If successful, the function
        # sets the CCD and subfield parameters correctly in the 'sim' object.
        numColSubfield = sim["SubField/NumColumns"]
        numRowSubfield = sim["SubField/NumRows"]
        raTargetRad    = np.deg2rad(self.df['ra'])
        decTargetRad   = np.deg2rad(self.df['dec'])
        
        self.isOnCCD = sim.setSubfieldAroundSkyCoordinates(raTargetRad, decTargetRad,
                                                           numColSubfield, numRowSubfield,
                                                           normal=normal)
        if self.isOnCCD is False:
            if self.verbose > 0:
                message  = (f"PIC {self.df.PIC} (subfield {self.targetNo}) " +
                            'do not fall on any of the CCDs for ' +
                            f'N-CAM {self.group}.{self.camera} and Q{self.quarter}!')
                errorcode('warning', message)
            # If requested keep track on non-observability
            if self.statistics: self.track_statistics(flag=1)
            exit()

        # If the psf is MappedFromFile we need to include mapped field distortion
        if sim["PSF/Model"] == "MappedFromFile":
            includeFieldDistortion = True
            mappedDistortion       = True
            pathToPsfFile          = sim["PSF/MappedFromFile/Filename"]
            distortionCoefficients = None
        elif (sim["Camera/IncludeFieldDistortion"] == "yes" or
              sim["Camera/IncludeFieldDistortion"] == "1"   or
              sim["Camera/IncludeFieldDistortion"] == True):
            includeFieldDistortion = True
            mappedDistortion       = False
            pathToPsfFile          = None 
            distortionCoefficients = sim["Camera/FieldDistortion/ConstantCoefficients"]
        else:
            includeFieldDistortion = False
            mappedDistortion       = False
            pathToPsfFile          = None 
            distortionCoefficients = None

        # Fetch info from YAML file since setting the subfied figure the inputfile
        self.raPlatformDeg            = sim["Platform/Orientation/Angles/RAPointing"]
        self.decPlatformDeg           = sim["Platform/Orientation/Angles/DecPointing"]
        self.solarPanelOrientationDeg = sim["Platform/Orientation/Angles/SolarPanelOrientation"]
        self.tiltTelescopeDeg         = sim["Telescope/TiltAngle"]
        self.azimuthTelescopeDeg      = sim["Telescope/AzimuthAngle"]
        # Transform to radians
        raPlatformRad            = np.deg2rad(self.raPlatformDeg)
        decPlatformRad           = np.deg2rad(self.decPlatformDeg)
        solarPanelOrientationRad = np.deg2rad(self.solarPanelOrientationDeg)
        tiltTelescopeRad         = np.deg2rad(self.tiltTelescopeDeg)
        azimuthTelescopeRad      = np.deg2rad(self.azimuthTelescopeDeg)
        # Hardware parameters
        pixelSize       = float(sim["CCD/PixelSize"])
        focalLength     = float(sim["Camera/FocalLength/ConstantValue"]) * 1000.0  # [m] -> [mm]
        focalPlaneAngle = np.deg2rad(float(sim["Camera/FocalPlaneOrientation/ConstantValue"]))

        # Fetch focal plane coordinates [mm]
        self.xFP, self.yFP = rf.skyToFocalPlaneCoordinates(raTargetRad, decTargetRad,
                                                           raPlatformRad, decPlatformRad,
                                                           solarPanelOrientationRad,
                                                           tiltTelescopeRad, azimuthTelescopeRad,
                                                           focalPlaneAngle, focalLength)
        if includeFieldDistortion == True or includeFieldDistortion == "yes":
            if mappedDistortion:
                self.xFP, self.yFP = rf.mappedUndistortedToDistortedFocalPlaneCoordinates(self.xFP, self.yFP,
                                                                                          pathToPsfFile, focalLength)
            else: 
                self.xFP, self.yFP = rf.undistortedToDistortedFocalPlaneCoordinates(self.xFP, self.yFP,
                                                                                    distortionCoefficients,
                                                                                    focalLength)
        
                
        # Fetch CCD code and pixel coordinates (account for field distortion if included)
        infoCCD = rf.getCCDandPixelCoordinates(raTargetRad, decTargetRad,
                                               raPlatformRad, decPlatformRad,
                                               solarPanelOrientationRad,
                                               tiltTelescopeRad, azimuthTelescopeRad,
                                               focalPlaneAngle, focalLength, pixelSize,
                                               includeFieldDistortion, normal=normal,
                                               mappedDistortion=mappedDistortion,
                                               distortionCoefficients=distortionCoefficients,
                                               pathToPsfFile=pathToPsfFile)
        self.ccdCode, self.xCCD, self.yCCD = infoCCD[0], infoCCD[1], infoCCD[2]
        
        # Only continue if ccdCode is found
        
        if self.ccdCode:

            # Check if string is F-CAM
            if len(self.ccdCode) == 2:
                self.ccdCode = self.ccdCode[0]

            # Add CCD time-shift to time points
            self.timeStart += float(sim['CCDPositions/TimeShift'][int(self.ccdCode)-1])
            self.time = np.arange(self.numExposures) * self.cadence + self.timeStart
        else:
            if self.verbose > 0:
                errorcode('warning', 'Star falls within a CCD gap!')
            # If requested keep track on non-observability
            if self.statistics: self.track_statistics(flag=2)                
            exit()

        # Calculate radial distance of coordinate away from OA
        self.rOA = np.rad2deg(rf.gnomonicRadialDistanceFromOpticalAxis(self.xFP, self.yFP,
                                                                       focalLength));
        if self.rOA > 19.0:
            if self.verbose > 0:
                message  = (f"PIC {self.df.PIC} (subfield {self.targetNo}) " +
                            f'is outside camera FOV (d={self.rOA:.2f} deg) ' +
                            f'for N-CAM {self.group}.{self.camera} and Q{self.quarter}!')
                errorcode('warning', message)
            # If requested keep track on non-observability
            if self.statistics: self.track_statistics(flag=3)
            exit()

        # If requested keep track on non-observability
        if self.statistics:
            self.track_statistics(flag=0)
            exit()
            
        # Create data frame for printing and saving
        c = ['PIC', 'ra [deg]', 'dec [deg]', 'mag',
             'CCD', 'xCCD [pix]', 'yCCD [pix]',
             'rOA [deg]', 'xFP [mm]', 'yFP [mm]', 'Ncon']
        d = {'PIC': [self.df['PIC']], 'mag': [self.df['mag']],
             'ra [deg]': [self.df['ra']], 'dec [deg]': [self.df['dec']],
             'CCD': [self.ccdCode], 'xCCD [pix]': [self.xCCD], 'yCCD [pix]': [self.yCCD],
             'rOA [deg]': [self.rOA], 'xFP [mm]': [self.xFP], 'yFP [mm]': [self.yFP],
             'Ncon': self.numCon}
        self.df0 = pd.DataFrame(d, columns=c)

        # Print data frame
        if self.verbose > 0:
            print('\nInformation about stellar target')
            print(self.df0)

        # Finito!
        return sim
        
        
        


    
    def create_seeds(self, sim):
        """
        Module to select and load the random seeds.
        """

        # Initialise random number generator after user or clock
        
        if not self.seed:
            rng = np.random.default_rng()
            self.seed = 123456789
        else:
            rng = np.random.default_rng(self.seed)

        # Seeds needs to be available for L1 pipeline
        self.seedJitter = self.seed * self.quarter
        self.seedTarget = self.seed * self.quarter + 1000 * self.targetNo + self.group * self.camera
            
        # Jitter (only relevant for red noise) only depends on the quarter
        if sim["Platform/UseJitter"] == 'yes' and sim["Platoform/JitterSource"] == 'RedNoise':
            sim["RandomSeeds/JitterSeed"] += self.seedJitter

        # Drift (only relevant for red noise) depends on quarter, group, and camera
        if sim["Telescope/UseDrift"] and sim["Telescope/DriftSource"] == 'RedNoise':
            sim["RandomSeeds/DriftSeed"] += (self.seed * self.quarter + 4 * self.group +
                                             24 * self.camera + 100000)

        # Lastly, the remaining seeds needs to differ for all times and for all subfields
        if self.seed:
            sim["RandomSeeds/ReadOutNoiseSeed"] += self.seedTarget
            sim["RandomSeeds/PhotonNoiseSeed"]  += self.seedTarget
            sim["RandomSeeds/DarkSignalSeed"]   += self.seedTarget
            sim["RandomSeeds/FlatFieldSeed"]    += self.seedTarget
            sim["RandomSeeds/CosmicSeed"]       += self.seedTarget
        else:
            sim["RandomSeeds/ReadOutNoiseSeed"] += rng.integers(1e9, size=1)[0]
            sim["RandomSeeds/PhotonNoiseSeed"]  += rng.integers(1e9, size=1)[0]
            sim["RandomSeeds/DarkSignalSeed"]   += rng.integers(1e9, size=1)[0]
            sim["RandomSeeds/FlatFieldSeed"]    += rng.integers(1e9, size=1)[0]
            sim["RandomSeeds/CosmicSeed"]       += rng.integers(1e9, size=1)[0]



            



            
    def create_inputfiles(self, sim):
        """
        Function to create ascii input files for PlatoSim.
        """
        
        # SAVE STELLAR CATALOGS AND TARGET LISTS
        
        # Save catalog and load it into the inputfile
        self.starCatalogFile = f'{self.outputDir}/{self.outputFileName}.cat'        
        sim.createStarCatalogFile(self.ds.ra, self.ds.dec, self.ds.mag, self.ds.ids,
                                  self.starCatalogFile)
        
        # Print catalogue
        if self.verbose and not self.fullFrame:
            print('\nStar catalog used in simulation')
            df1 = pd.DataFrame({'RA [deg]': self.ds.ra,
                                'Dec [deg]': self.ds.dec,
                                'P [mag]': self.ds.mag,
                                'Dis [pix]': np.append(0, self.dc['dis'])/15.})
            print(df1)
        

        # SAVE LIST OF VARIABLE SOURCES

        # Automatically activate varsource (if the user forgets)
        if self.varSourceFile or self.varSourceList:
            sim['Sky/IncludeVariableSources'] = True

        # NOTE if a user-defined filename for the variable data is parsed
        # then a variable source list is created automatically
        # TODO binary files are not yet compatible with PlatoSim!
        if sim['Sky/IncludeVariableSources'] is True:

            # First priority to varSourceList -> Allows variability for all stars
            if self.varSourceList:

                # Make sure that wrong paths are detected
                self.varSourceList = Path(self.varSourceList).resolve()
                if self.varSourceList.is_file():
                    sim["Sky/VariableSourceList"] = self.varSourceList
                else:
                    errorcode('error', 'VariableSourceList do not exist, check file path!')
                    
            # Second priority to varSourceFile -> Allows variability for target only
            elif self.varSourceFile:

                # If varSourceFile is parsed
                self.varSourceFile = Path(self.varSourceFile).absolute()
                if self.varSourceFile.is_file():
                    # Automaticallt create and save varSourceList to file
                    self.varSourceList = f'{self.outputDir}/{self.outputFileName}.var'
                    sim.createVariableSourceList('1', str(self.varSourceFile), self.varSourceList)
                else:
                    errorcode('error', 'VariableSourceFile do not exist, check file path!')

        
        # SAVE PHOTOMETRY FILE
        
        # NOTE if a user defined file name for the photometry file is parsed
        # then a photometry file list is created automatically
        if sim['Photometry/IncludePhotometry'] is True:
            photometryList = self.inputDir.joinpath('photometry.txt')
            if os.path.exists(photometryList) is False:
                np.savetxt(photometryList, np.array([]), header='1', comments='')
            sim["Photometry/TargetFileName"] = photometryList


            


            
            
    def show_subfield(self, sim):
        """
        Function to show 
        (1) where the target star is situated in the CCD focal plane
        (2) the first subfield/imagette of the simulation
        This function exit the simulation.
        """

        # Print to bash
        if self.verbose > 0:
            errorcode('message', f'\n[PlatoSim]: Visualizing simulation for N-CAM ' +
                      f'{self.group}.{self.camera} Q{self.quarter}\n')
        
        # Control the content of the hdf5 output files
        sim.turnOffAllOutput()
        sim["ObservingParameters/NumExposures"]      = 1
        sim["ControlHDF5Content/WritePixelMaps"]     = True
        sim["ControlHDF5Content/WriteStarPositions"] = True
        
        # Add photometric mask to plot if available
        if sim['Photometry/IncludePhotometry']: mask = 1
        else: mask = None

        # Run simulation for first image cadence
        self.outputSimName = self.outputDir.joinpath(self.outputFileName)
        sim.outputDir = self.outputDir

        f = sim.run(removeOutputFile=self.overwrite)
            
        # Plot star in CCD focal plane
        if not self.fullFrame:
            fig = plt.figure(figsize=(12,10))
            drawStarInCCDfocalPlane(fig, sim,
                                    self.df0['xCCD [pix]'][0], self.df0['yCCD [pix]'][0],
                                    self.df0['CCD'][0], self.group,
                                    self.raPlatformDeg, self.decPlatformDeg,
                                    self.solarPanelOrientationDeg)

        # Show subfield for first cadence
        if self.fullFrame:
            figsize = (12,12)
            showStarPositions = False
            showMaskOfStarID  = False
            showGrid          = False
            title             = False
            imgScale          = "auto"
            cmap              = 'cubehelix'
            clipPercentile    = 1
            
        else:
            figsize = (6,6)
            showStarPositions = 'PIC'
            showMaskOfStarID  = '1'
            title = f'Imagette of PIC {int(self.df.PIC)} ({float(self.df.mag):.2f} mag)'
            clipPercentile    = 2
            imgScale          = "auto"
            cmap              = 'magma' #'gist_stern'
            showGrid          = True
            
        # Check that if any stars are detected
        try:
            f.getStarCoordinates(self.beginExposureNr)
        except:
            errorcode('warning', 'No stars detected within the subfield!')
            showStarPositions = False
            
        # Plot the subfield    
        f.showImage(self.beginExposureNr,
                    showStarPositions=showStarPositions,
                    clip=clipPercentile,
                    showMaskOfStarID=mask,
                    useTitle=title,
                    colorMap=cmap,
                    colorBar=True,
                    imgScale=imgScale,
                    showGrid=showGrid,
                    figsize=figsize)

        # Remove the output files
        if os.path.isfile(str(self.outputSimName) + '.hdf5'):
            os.system(f'rm {self.outputDir}/{self.outputFileName}*')

            

            

    def run_sim_normal(self, sim):
        """
        Module to run PlatoSim only.
        """

        # Print to bash
        if self.verbose > 0:
            tracemalloc.start()
            if self.fullFrame:
                ccdID = f' CCD {self.ccdCode}'
            else:
                ccdID = ''
            if self.groupID == 'Fast':
                errorcode('message', '\n[PlatoSim]: Simulating' +
                          f'{ccdID} F-CAM {self.camera} Q{self.quarter} ' + 
                          f'for {self.numExposures} exposures')
            else:
                errorcode('message', f'\n[PlatoSim]: Simulating' +
                          f'{ccdID} N-CAM {self.group}.{self.camera} Q{self.quarter} '
                          f'for {self.numExposures} exposures')

        # Select the number of exposures
        sim["ObservingParameters/NumExposures"] = self.numExposures

        # Save images if animation is requested
        if self.animation:
            sim["ControlHDF5Content/WritePixelMaps"]     = True
            sim["ControlHDF5Content/WriteStarPositions"] = True
            
        # Run simulation
        sim.outputDir = self.outputDir
        simFile = sim.run(removeOutputFile=self.overwrite, logLevel=self.verbose_platosim)
        
        # Common files to always remove (unless debug mode)
        if self.isOnCCD:
            if self.varSourceFile:
                os.remove(self.varSourceList)
            if self.sample is None and self.verbose < 3:
                if os.path.isfile(str(self.outputSimName) + '.log'):
                    os.remove(str(self.outputSimName) + '.log')
                if os.path.isfile(str(self.outputSimName) + '.yaml'):
                    os.remove(str(self.outputSimName) + '.yaml')
                if os.path.isfile(self.starCatalogFile):
                    os.remove(self.starCatalogFile)

        # Define output file name
        outputFile = f'{self.outputSimName}.hdf5'
        
        # Save full-frame catalogue for first exposure
        if self.fullFrame:
            f = SimFile(f'{self.outputSimName}.hdf5')
            PIC, row, col, xFP, yFP, flux = f.getStarCoordinates(self.beginExposureNr) 
            # Select detected stars
            df = self.dx.iloc[PIC]
            df = df.reset_index()
            # Add stellar positions.
            df['xCCD'] = col - 0.5
            df['yCCD'] = row - 0.5
            df['xFP']  = xFP
            df['yFP']  = yFP
            # Only keep stars within the rOA FOV
            focalLength = float(sim["Camera/FocalLength/ConstantValue"]) * 1000
            df['rOA'] = np.rad2deg(rf.gnomonicRadialDistanceFromOpticalAxis(df.xFP, df.yFP,
                                                                            focalLength))
            df = df[df.rOA <= 19.55]
            df = df.drop(columns=['index'])
            df = df.reset_index()
            # Sort index and starID
            df = df.rename(columns={'level_0': 'starID'})
            df = df.drop(columns=['index'])
            # Save to file
            df.to_feather(self.outputSimName + '.ftr')

        # Make a animation if requested
        if self.animation:
            
            # Adjust number of images to skip and frame rate
            if self.cadence == 25.0:
                fps, nskip = 50, 1000
            elif self.cadence == 50.0:
                fps, nskip = 25, 500
            elif self.cadence == 600.0:
                fps, nskip = 25, 50
            plotSubfieldAnimation(outputFile,
                                  outputFileName=str(self.outputSimName),
                                  cadence=self.cadence,
                                  frameRate=fps,
                                  skipNimages=nskip,
                                  numImages=False,
                                  colorMap="gist_stern",
                                  clipPercentile=8.0, 
                                  showStarPositions='PIC',
                                  showMaskOfStarID='1',
                                  useTitle=True,
                                  showGrid=True,
                                  figsize=(6,6))
            
        # Resources
        if self.verbose > 0:
            # Execution time of module
            self.tocPlatoSim = datetime.datetime.now() - self.tic
            self.tic = datetime.datetime.now()
            # Max RAM memory [Mb]
            self.memRamPlatoSim = round(tracemalloc.get_traced_memory()[1]/1e3)
            tracemalloc.stop()
            # Storage memory of HDF5 file [Mb]
            self.memDiskPlatoSim = round(Path(outputFile).stat().st_size/1e6)

        

            
    #--------------------------------------------------------------#
    #                    L1 PIPELINE MODULES                       #
    #--------------------------------------------------------------#

    def control_hdf5(self):
        """
        Module to control HDF5 content for L1 pipeline.
        """

        # Include HDF5 content
        sim["ControlHDF5Content/GroupByExposure"]             = True
        sim["ControlHDF5Content/WritePixelMaps"]              = True
        sim["ControlHDF5Content/WriteBiasMaps"]               = False
        sim["ControlHDF5Content/WriteSmearingMaps"]           = True
        sim["ControlHDF5Content/WriteFlatfieldMap"]           = True
        sim["ControlHDF5Content/WriteThroughputMaps"]         = True
        sim["ControlHDF5Content/WriteTransmissionEfficiency"] = True
        sim["ControlHDF5Content/WriteBackgroundMap"]          = False
        sim["ControlHDF5Content/WriteCTI"]                    = False        
        sim["ControlHDF5Content/WriteSubPixelImages"]         = False
        sim["ControlHDF5Content/WriteHighResolutionPSF"]      = True
        sim["ControlHDF5Content/WriteACS"]                    = True
        sim["ControlHDF5Content/WriteTelescopeACS"]           = False
        sim["ControlHDF5Content/WriteStarCatalog"]            = True
        sim["ControlHDF5Content/WriteStarPositions"]          = True
        sim["ControlHDF5Content/WriteGhostPositions"]         = False
        sim["ControlHDF5Content/WriteCosmics"]                = True
        
        # Check for high res mapped PSF
        if sim["PSF/Model"] == 'MappedFromFile':
            sim["ControlHDF5Content/WriteDiffusedPSF"] = True
        else:
            sim["ControlHDF5Content/WriteDiffusedPSF"] = False




            
    def run_microscan(self, sim):
        """
        Module to run a microscan sequence with PlatoSim.
        """

        # Print to bash
        if self.verbose > 0:
            errorcode('module', '\nMicroscanning & PSF inversion')

        # Check if mapped PSFs are used and apply correct resolution
        if sim["PSF/Model"] == 'MappedFromFile':
            sim["SubField/SubPixels"] = 64
        else:
            sim["SubField/SubPixels"] = 128
            
        # Prepare for the Archimedean jitter spiral pattern
        nimages = 428
        sim["Platform/UseJitter"]               = 'yes'
        sim["Platform/JitterSource"]            = 'FromFile'
        sim["ObservingParameters/NumExposures"] = nimages
        sim['ObservingParameters/CycleTime']    = 25
        
        # Time of microscan simulation needs to match quarter for PlatoSim to run successfully
        # NOTE Here a new file is created with a appropriate time column to match the observation
        # NOTE This is done once per quarter and CCD time-shift and the file is saved to the
        # simulation folder
        spiralFileName     = 'microscan_spiral_8Hz_3h_BC.txt'
        spiralFileBase     = f'{self.platosimInputDir}/{spiralFileName}'
        spiralFileQuarterN = f'{self.microscanDirStarID}/{spiralFileName}'

        # Download microscan file if first run of L1 pipeline
        if not Path(spiralFileBase).is_file():
            print(f'Downloading miscroscanning file..')
            ut.downloadFromFTP(filename=spiralFileName, outputDir=self.platosimInputDir, server='plato')
            
        if self.quarter == 1 and int(self.df0['CCD']) == 1:
            os.system(f'cp {spiralFileBase} {spiralFileQuarterN}')
        else:
            # Load file
            data = np.loadtxt(spiralFileBase)
            t, x, y, z = data[:,0], data[:,1], data[:,2], data[:,3]
            # Generate new time column
            cadence = np.diff(t)[0]  # [s -> 8 Hz]
            t = np.arange(len(t)) * cadence + self.timeStart
            # Save data to input folder
            np.savetxt(spiralFileQuarterN, np.transpose([t,x,y,z]), fmt=['%.3f', '%.9f', '%.9f', '%.9f'])

        # Load Archimedean spiral
        sim["Platform/JitterFileName"] = spiralFileQuarterN

        # HDF5 content to always exclude
        sim["Telescope/UseDrift"]                     = False
        sim["Sky/IncludeVariableSources"]             = False
        sim["Photometry/IncludePhotometry"]           = False
        sim["ControlHDF5Content/WriteSubPixelImages"] = False
        
        # HDF5 content to always include
        sim["ControlHDF5Content/WritePixelMaps"]              = True
        sim["ControlHDF5Content/WriteBiasMaps"]               = True
        sim["ControlHDF5Content/WriteSmearingMaps"]           = True
        sim["ControlHDF5Content/WriteThroughputMaps"]         = True
        sim["ControlHDF5Content/WriteFlatfieldMap"]           = True
        sim["ControlHDF5Content/WriteHighResolutionPSF"]      = True
        sim["ControlHDF5Content/WriteTransmissionEfficiency"] = True
        sim["ControlHDF5Content/WriteStarPositions"]          = True
        sim["ControlHDF5Content/WriteACS"]                    = True
        sim["ControlHDF5Content/WriteCosmics"]                = True
        sim["ControlHDF5Content/WriteStarCatalog"]            = True
        sim["ControlHDF5Content/WriteTelescopeACS"]           = True
        sim["ControlHDF5Content/WriteCTI"]                    = True

        # If mapped PSF is used the diffused PSFs need to be saved
        if sim["PSF/Model"] == 'MappedFromFile':
            sim["ControlHDF5Content/WriteDiffusedPSF"] = 'yes'
        
        # Save catalog and load it into the inputfile
        # NOTE we need to replace the target name for correct handling by the LESIA pipeline
        numStar = self.numCon + 1
        self.ds.ids = np.arange(self.targetNo, self.targetNo + numStar, 1) + 1
        starCatalogFile = f'{self.microscanSimName}.coo'
        np.savetxt(starCatalogFile, self.ds, fmt=['%11.6f', '%11.6f', '%8.4f', '%i'])
        sim["ObservingParameters/StarCatalogFile"] = starCatalogFile

        # MICROSCANNING SIMULATION

        if self.verbose > 0:
            errorcode('message', f'\n[PlatoSim]: Simulating {nimages} imagettes along Archimedean spiral')
        sim.outputDir = self.microscanDirStarID
        simFile = sim.run(removeOutputFile=self.overwrite, logLevel=self.verbose_platosim)

        # Execution time of module
        if self.verbose > 0:
            self.tocMicroscan = datetime.datetime.now() - self.tic
            self.tic = datetime.datetime.now()
            
        # PRE-PROCESSING

        # Change directory needed to execute scripts
        os.chdir(self.microscanDir)
        
        # Run pre-processing
        if self.verbose > 0:
            errorcode('message', '\n[pproc]: Pre-processing imagettes')
        cmd = os.system(f'{self.platoLib}/pproc.py ' +
                        f'--platosim --auto-bg -f {self.starID} {self.devnull}')
        if cmd != 0: errorcode('error', 'pproc.py failed due to the above error!')

        # EXTRACT CONTAMINANTS

        # Run contaminant extraction
        if self.verbose > 0:
            errorcode('message', f'\n[extract_contaminants]: Model contaminant stars')
            print(f'Include contaminats with dmag < {self.conDeltaMag} from target')
        cmd = os.system(f'{self.platoLib}/extract_contaminant.py ' +
                        f'-D {self.conDeltaMag} -e {self.conFluxError} ' + 
                        f'-s {self.seedTarget} {self.starID} {self.starID} {self.devnull}')
        if cmd != 0: errorcode('error', 'extract_contaminant.py failed due to the above error!')
        if self.verbose > 0:
            print('Modelling of contaminants done')
            
        # PSF INVERSION
        
        # Find Regularization parameter for each star
        PV   = -0.34  # P-V magnitude offset
        Vmag = self.df['mag'] - PV
        regs = np.format_float_scientific(10.**(0.51 * Vmag - 14.61))
        
        # Make sure to only use one thread since we will use the HPC
        os.system('export OMP_NUM_THREADS=1')
        
        # Run the inversion module
        # NOTE Type of microscanning (-t)
        # NOTE Sub-pixel resolution of invertion (-p)
        # NOTE Sub-pixel resolution of inverted PSF (-l)
        # NOTE Sub-pixel resolution of original PSF (-r)
        # NOTE Regularisation parameter for the wPRLS method (-u)
        # NOTE Number of elementary steps over which to calculate averaged positions (-N)
        if self.verbose > 0:
            errorcode('message', '\n[invert_parabolic_multi]: Run the PSF inversion')
        cmd = os.system(f'{self.platoLib}/invert_parabolic1_multi ' +
                        f'-Q -t continuous -m PRLS ' +
                        f'-N 1 -r 128 -l 128 -p {self.bsres} -u {regs} ' +
                        f'-d . -i {self.starID} -q {self.starID}/{self.starID}_offsets.txt ' + 
                        f'-o inversion {self.devnull}')
        if cmd != 0: errorcode('error', 'invert_parabolic_multi failed due to the above error!')

        # TODO check how well the inversion is!
        
        # Execution time of module
        if self.verbose > 0:
            self.tocInversion = datetime.datetime.now() - self.tic
            self.tic = datetime.datetime.now()


            

        

        
    def run_L1_onground(self):
        """
        Module to for the on-ground L1 pipeline processing chain.
        """

        # Print to bash
        if self.verbose > 0:
            errorcode('module', '\nOn-ground L1 pipeline')

        # Change directory needed to execute scripts
        os.chdir(self.pipelineDir)

        # PRE-PROCESSING

        if self.verbose > 0:
             errorcode('message', '\n[pproc]: Pre-processing imagettes')
        cmd = os.system(f'{self.platoLib}/pproc.py ' +
                        f'--platosim --auto-bg -f {self.starID} {self.devnull}')
        if cmd != 0: errorcode('error', 'pproc.py failed due to the above error!')
        
        # PSF FIITING

        # NOTE using Dierckx's knot distribution (-K 1)
        # NOTE using Levenberg-Marquardt minimization method (default: -M 0)
        # NOTE using B-spline resolution of 10 (-b) matching the PSF resolution!
        # NOTE PRNU knowledge errorin % (-p)
        if self.verbose > 0:
            errorcode('message', '\n[psffit]: PSF fitting for light curve generation')
        # Change directory needed to execute scripts

        #path = self.path / '../pipeline'
        cmd = os.system(f'{self.platoLib}/psffit.py ' +
                        f'-K 1 -b {self.bsres} --seed {self.seedTarget} ' +
                        f'-F {self.tarFluxError} -s {self.tarAbsCenError} -p {self.prnuError} ' + 
                        f'-f {self.microscanDirInvers}/{self.starID}_PRLS.vec ' +
                        f'-o {self.starID} {self.starID} {self.devnull}')
        if cmd != 0: errorcode('error', 'psffit.py failed due to the above error!')
        
        # PROLOGUE
                
        # Execution time of module
        if self.verbose > 0:
            self.tocOnground = datetime.datetime.now() - self.tic
            self.tic = datetime.datetime.now()





        

    def run_L1_onboard(self):
        """
        Module to for the on-board L1 pipeline processing chain.
        """

        # Print to bash
        if self.verbose > 0:
            errorcode('module', '\nOn-board L1 pipeline')
            
        # Change directory needed to execute scripts
        os.chdir(self.pipelineDir)

        # PRE-PROCESSING

        if self.verbose > 0:
            errorcode('message', '\n[pproc]: Pre-processing imagettes')
        cmd = os.system(f'{self.platoLib}/pproc.py ' +
                        f'--platosim --auto-bg -f {self.starID} {self.devnull}')
        if cmd != 0: errorcode('error', 'pproc.py failed due to the above error!')
        
        # NOTE for the input parameters:
        # - Binary mask is default (-T 1)
        # - Write mask FITS files (-M)
        # - Calculating SPR vs. time (--spr_tot) -> only binary masks!
        # - Using B-spline PSF resolution (--bsres)
        # - Using Dierckx's knot distribution (-K 1)
        # - Update mask every 14 days (--update-period 48384)
        # - Update mask given 0.1 pixel threshold (--update-thres 0.1)
        # Conditions for a mask update are:
        # - the displacement since the last update exceeds the threshold
        # - the current exposure is = last update time + update period
        # - the exposure number must be a multiple of 24 such that the mask update
        #   always occurs at the beginning of a 600s cycle
        if self.verbose > 0:
            errorcode('message', '\n[lightcurve.py]: Aperture photometry ala Marchiori+2019')
        cmd = os.system(f'{self.platoLib}/lightcurve.py ' +
                        f'-M --input-hdf5 --spr_tot ' + 
                        f'--include-contaminants --add_chromatic_abberation ' +
                        f'--update-period {self.maskUpdateRate} ' +
                        f'--update-thres {self.maskUpdateThres} ' +
                        f'--bsres {self.bsres} ' +
                        f'-I {self.microscanDirInvers}/{self.starID}_PRLS.vec ' +
                        f'-B {self.starID} -o {self.starID} {self.starID} {self.devnull}')
        if cmd != 0: errorcode('error', 'lightcurve.py failed due to the above error!')

        # JITTER AND DRIFT CORRECTION

        if not self.jitterDriftOff:
            if self.verbose > 0:
                errorcode('message', '\n[jittercorrection.py]: Jitter & Drift Correction')
            cmd = os.system(f'{self.platoLib}/jittercorrection.py ' +
                            f'--add_chromatic_abberation ' +
                            f'-W 128 -r 128 ' +
                            f'--bsres {self.bsres} --seed {self.seedJitter} ' + 
                            f'-f {self.prnuError} -b {self.tarAbsCenError} -a {self.conDeltaMag} ' +
                            f'-I {self.microscanDirInvers}/{self.starID}_PRLS.vec ' +
                            f'-o {self.starID} {self.starID} {self.starID} {self.devnull}')
            if cmd != 0: errorcode('error', 'psffit.py failed due to the above error!')
            
        # PROLOGUE
                
        # Execution time of module
        if self.verbose > 0:
            self.tocOnboard = datetime.datetime.now() - self.tic
            self.tic = datetime.datetime.now()


        


    #--------------------------------------------------------------#
    #                            OUTPUTS                           #
    #--------------------------------------------------------------#

    def create_sim_table(self, odir):

        """Module to create a overview table of the simulation.
        """

        # Write PlatoSim info to a table
        filename = f'{odir}/{self.outputFileName}_table.ftr'
        data = {"ID":      self.targetNo+1,
                "PIC":     self.df.PIC,
                "ra":      self.df.ra,
                "dec":     self.df.dec,
                "mag":     self.df.mag,
                "group":   self.group,
                "camera":  self.camera,
                "quarter": self.quarter,
                "rOA":     self.rOA,
                "xFP":     self.xFP,
                "yFP":     self.yFP,
                "ccd":     self.ccdCode,
                "xCCD":    self.xCCD,
                "yCCD":    self.yCCD,
                "ncon":    self.numCon}
        df1 = pd.DataFrame(data, index=[0])
        df1.to_feather(filename)

            
            
    

    def sort_output_normal(self):

        # Create a info table of simulation
        if not self.fullFrame:
            self.create_sim_table(self.outputDir)

        # Give full read and write access to output files
        os.system(f'chmod 755 {self.outputSimName}*')
                    
        # Compress files
        if self.compress and os.path.isfile(self.outputSimName + '.hdf5') and not self.sample:
            if self.verbose > 0:
                errorcode('module', '\nRestructuring data output\n')
                print('Compressing files')
            os.system(f'zip -j {self.outputSimName}.zip {self.outputSimName}* ' +
                      f'{self.devnull}')
            # Give read and write access to file
            os.system(f'chmod 755 {self.outputSimName}.zip')
            # Remove non-compressed files
            os.remove(self.outputSimName + '.hdf5')
            if not self.fullFrame:
                os.remove(self.outputSimName + '_table.ftr')

        # If requested move file to final output directory (for cluster)
        if self.hpcDir:
            os.system(f'mv {self.outputSimName}.* {self.hpcDir}')

            
        # Execution time of module
        if self.verbose > 0:
            self.tocPrologue = datetime.datetime.now() - self.tic
            self.tic = datetime.datetime.now()




            
    def sort_output_pipeline(self):

        """Module to move output files and delete redundant files.
        """

        if self.verbose > 0:
            errorcode('module', '\nPrologue')
        
        if self.verbose > 0:
            errorcode('message', '\nRestructuring data output')
            print(f'L1 light curve is saved to {self.outputDirStarIDnew}')

        # Select prefix-files
        self.outputFileName = f'/{self.starID}_{self.obsPrefix}'
        prefixInversion = self.microscanDirInvers + f'/{self.starID}'
        prefixStarIDtar = self.outputDirStarIDsim + f'/000000001'
        prefixStarIDsim = self.outputDirStarIDsim + f'/{self.starID}'
        prefixStarIDnew = self.outputDirStarIDnew + self.outputFileName

        # Create a info table of simulation
        self.create_sim_table(self.outputDirStarIDnew)

        # Rewrite data time series
        if args.sample == 'P1':
            cols = ['flux', 'cx', 'cy', 'bg', 'flux_err', 'cx_err', 'cy_err', 'bg_err',
                    'chi2', 'iter', 'lamb']
            try:
                df = pd.read_csv(prefixStarIDtar + '.dat', delimiter=' ', comment='#',
                                 names=cols, usecols=np.arange(0, len(cols), 1))
            except:
                # Open the file in append & read mode ('a+')
                with open(self.inputDir / 'failed.txt', 'a+') as f:
                    # Move read cursor to the start of file.
                    f.seek(0)
                    # If file is not empty then append '\n'
                    if len(f.read(100)) > 0: f.write("\n")
                    # Append text at the end of file
                    f.write(f'{int(self.starID)} {self.group} {self.camera} {self.quarter}')
                df = None
                
        if args.sample == 'P5':
            cols = ['flux', 'xc', 'yc', 'flux_cor']
            # If jitter/drift correction is not applied no JC correction neither
            if self.jitterDriftOff:
                df = pd.read_csv(prefixStarIDtar + '.dat', delimiter=' ', comment='#',
                                 names=cols[:3], usecols=cols[:3])
            else:
                df = pd.read_csv(prefixStarIDtar + '-jc.dat', delimiter=' ', comment='#',
                                 names=cols, usecols=cols)
            # Move the SPR file
            shutil.move(prefixStarIDtar + '-sprtot.dat', prefixStarIDnew + '.spr')
            # Move mask files number after mask update exposure
            maskfits = glob.glob(f'{prefixStarIDtar}**mask.fits')
            for i in range(len(maskfits)):
                shutil.move(maskfits[i], prefixStarIDnew + maskfits[i][-17:])
        
        # Add a proper time column
        if df is not None:
            df = pdAddColumn(df, self.time, 'time')
            if args.sample == 'P1': df = df.astype({'time':np.float64, 'iter':int})
            if args.sample == 'P5': df = df.astype({'time':np.float64})
            # Feather format needs to be indiced!
            df = df.reset_index()
            # Save new data frame
            df.to_feather(prefixStarIDnew + '.ftr')
            # Move files to new data directory
            shutil.move(prefixInversion + f'_PRLS_invert.log', prefixStarIDnew + '.invert')

            
        # Remove microscan-starID and simulation folder (and all its content)
        if self.verbose != 3:
            shutil.rmtree(self.microscanDirStarID)
            shutil.rmtree(self.microscanDirInvers)
            shutil.rmtree(self.outputDirStarIDsim)

        # Give full read/write access
        os.system(f'chmod 755 {prefixStarIDnew}*')
        
        # Compress files
        if self.compress:
            if self.verbose > 0:
                print('Compressing files')
            os.system(f'zip -j {prefixStarIDnew}.zip {prefixStarIDnew}* {self.devnull}')
            os.system(f"find {self.outputDirStarIDnew} -type f -not -name '*.zip' -delete {self.devnull}")
            os.system(f'chmod 755 {prefixStarIDnew}.zip')

        # If requested move file to final output directory (for cluster)
        if self.hpcDir:
            os.system(f'mv {prefixStarIDnew}.* {self.hpcDir}')
            
        # Execution time of module
        if self.verbose > 0:
            self.tocPrologue = datetime.datetime.now() - self.tic
            self.tic = datetime.datetime.now()

            
            


    def resources(self):

        """Module to print resources used by PLATOnium.
        """
        if not args.plot:
            errorcode('message', '\nSimulation statistics')
            print(f'Max RAM memory for PlatoSim      : {self.memRamPlatoSim} MB')
            print(f'Storage memory for PlatoSim      : {self.memDiskPlatoSim} MB')
            print(f'Execution time for PlatoSim      : {self.tocPlatoSim} [hh:mm:ss]')
            if self.sample:
                print(f'Execution time for Microscanning : {self.tocMicroscan} [hh:mm:ss]')
                print(f'Execution time for PSF inversion : {self.tocInversion} [hh:mm:ss]')
                if self.sample == 'P1':
                    print(f'Execution time for L1 On-ground  : {self.tocOnground} [hh:mm:ss]')
                if self.sample == 'P5':
                    print(f'Execution time for L1 On-board   : {self.tocOnboard} [hh:mm:ss]')
                print(f'Execution time for Prologue      : {self.tocPrologue} [hh:mm:ss]')
            print('------------------------------------------------------------')
            print(f'Total execution time             : {datetime.datetime.now() - self.tic0} [hh:mm:ss]\n')



            
#==============================================================#
#               PARSING COMMAND-LINE ARGUMENTS                 #
#==============================================================#

parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
                                 #description=errorcode('software', '\nPLATOnium'))

parser.add_argument('-p', '--plot',      action='store_true',   help='Flag to plot 1st subfield (this will not save the simulation)')
parser.add_argument('-w', '--overwrite', action='store_true',   help='Flag to overwrite an existing HDF5 output file that may exist')
parser.add_argument('-a', '--animation', action='store_true',   help='Flag to generate a animation of the subfield simulation')
parser.add_argument('-v', '--verbose', metavar='INT', type=int, help='Verbosity level [0, 1, 3] (Default: 1)')

man_group = parser.add_argument_group('MANDATORY PARAMETERS')
man_group.add_argument('starID',   type=int, help='Star ID in target list (or for --fullframe CCD in [1,2,3,4])')
man_group.add_argument('groupID',  type=int, help='Camera group ID: [1, 2, 3, 4, Fast]')
man_group.add_argument('cameraID', type=int, help='N-CAM: [1, 2, 3, 4, 5, 6], F-CAM: [1, 2]')
man_group.add_argument('quarter',  type=int, help='Mission quarter: [1, 2, 3, 4, ..]')

out_group = parser.add_argument_group('I/O PARAMETERS')
out_group.add_argument('-i', '--infile', metavar='FILE', type=str, help='Path to YAML input file')
out_group.add_argument('-o', '--outdir', metavar='PATH', type=str, help='Output directory')
out_group.add_argument('-d', '--hpcdir', metavar='PATH', type=str, help='Final output directory for storage on cluster')
out_group.add_argument('--project',      metavar='NAME', type=str, help='Name project folder within $PLATO_WORKDIR')
out_group.add_argument('--yaml',         metavar='NAME', type=str, help='Name of yaml file within $PLATO_WORKDIR/input')
out_group.add_argument('--prefix',       metavar='NAME', type=str, help='Prefix filename extention')
out_group.add_argument('--starcat',      metavar='FILE', type=str, help='Path to stellar catalogue    -> see PlatoSim docs')
out_group.add_argument('--varfile',      metavar='FILE', type=str, help='Path to variable source file -> see PlatoSim docs')
out_group.add_argument('--varlist',      metavar='FILE', type=str, help='Path to variable source list -> see PlatoSim docs')
out_group.add_argument('--compress',     action='store_true',      help='Flag to compress output files')
out_group.add_argument('--statistics',   action='store_true',      help='Flag to compute the camera observability')

sim_group = parser.add_argument_group('SIM PARAMETERS')
sim_group.add_argument('--cadence', metavar='SEC',  type=float, help='Cadence or cycle time step for each exposure [seconds] (default: 25 s)')
sim_group.add_argument('--tdur',    metavar='DAY',  type=float, help='Total lenght of shortened quarter time series [days]')
sim_group.add_argument('--nexp',    metavar='NO.',  type=int,   help='Number of exposures of shortened quarter time series')
sim_group.add_argument('--bexp',    metavar='NO.',  type=int,   help='Number of exposure from beginning of quarter')
sim_group.add_argument('--pic',     metavar='ID',   type=int,   help='ID from the PIC (overwrites "starID")')
sim_group.add_argument('--seed',    metavar='INT',  type=int,   help='Option to bootstrap seeds in order to reproduce results')
sim_group.add_argument('--mask',    metavar='DAY',  type=float, help='Option to alter/overwrite the mask-update in inputfile [days]')
sim_group.add_argument('--mag',     metavar='PMAG', type=float, help='Option to change target magnitude (overwrites catalogue magnitude)')
sim_group.add_argument('--performance', metavar='MODE', type=str, help='Option to set basic input parameters (required, expected)')
sim_group.add_argument('--nocon',       action='store_true',    help='Flag ignore all stellar contaminants')
sim_group.add_argument('--reusejitter', action='store_true',    help='Flag to re-use AOCS jitter file across all quarters')
sim_group.add_argument('--fullframe',   action='store_true',    help='Flag to simulate full-frame CCD -> CCDcode = starID')


pip_group = parser.add_argument_group('L1 PIPELINE PARAMETERS')
pip_group.add_argument('--sample',   metavar='PIC',     type=str,   help='PIC sample to activate the correct L0-L1 pipeline chain [P1, P5]')
pip_group.add_argument('--con_dmag', metavar='MAG',     type=float, help='Relative magnitude inclusion threshold for contaminant(s) (Default: 5 mag)')
pip_group.add_argument('--con_ferr', metavar='PERCENT', type=float, help='Flux error assumption for Target')
pip_group.add_argument('--tar_ferr', metavar='PERCENT', type=float, help='Flux error assumption for Contaminant(s)')
pip_group.add_argument('--tar_cerr', metavar='PIXEL',   type=float, help='Centroid error assumption for Target')
pip_group.add_argument('--jit_off',  action='store_true', help='Turn off the jitter/drift correction')

args = parser.parse_args()

# Load and run modules
p = PLATOnium(args)
p.load_stars()
p.configure_output()
sim = p.init_sim()
p.create_seeds(sim)
p.create_inputfiles(sim)

if args.plot:
    # Only show imagette
    p.show_subfield(sim)
    
elif args.sample == 'P1':
    # Run the entire on-ground L0-L1 pipeline chain
    p.control_hdf5()
    p.run_sim_normal(sim)
    p.run_microscan(sim)
    p.run_L1_onground()
    p.sort_output_pipeline()
    
elif args.sample == 'P5':
    # Run the entire on-board L0-L1 pipeline chain
    p.control_hdf5()
    p.run_sim_normal(sim)
    p.run_microscan(sim)
    p.run_L1_onboard()
    p.sort_output_pipeline()
    
else:
    # Only run PlatoSim time series
    p.run_sim_normal(sim)
    p.sort_output_normal()
    
# Finito!
if args.verbose != 0:
    p.resources()
