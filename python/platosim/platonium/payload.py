#!/usr/bin/env python3

"""
This script is an integrated part of PlatoSim's toolkit PLATOnium.
Call this script prior to running your simulations with "platonium"
to generate the following files:

  - instrumentGaps.txt : Realistic distribution of spacecraft down times
  - instrumentPRE.txt  : Pointing errors for each mission quarter pointing
  - instrumentAPE.txt  : Alignment errors between cameras and the optical bench
  - instrumentTED.txt  : Long-term camera drift due to Thermo-Elastic Distortion
  - run.pbs & data.pbs : Job and parameterisation file for SLURM parallelisation
  - inputfile.yaml     : Copy and adjust a YAML file (only if it doesn't exist)

While using the "--project" output path, the files will be generated
directly into your working directory and will immediately be known
and used when running "platonium". 
"""
# TODO instrumentGain.txt : Gain differences for the each CCD (F and E side) and FEE

# Python standard
import os
import math
import shutil
import datetime
import argparse
from textwrap import dedent

# PlatoSim standard
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from pathlib import Path

# PlatoSim functions
import platosim.noise as ns
from platosim.utilities import errorcode, convertQuarterRange, getPointingField

#==============================================================#
#                         BEGIN CLASS                          #
#==============================================================#

class Payload(object):

    """Class to generate a realistic setup for multi-camera simulations.
    
    This class makes it straight forward to generate a realistic random
    model for various systematic instrumental effects. It likewise add
    all the necessary files needed to launch the simulations on a high
    performance computing cluster.
    """

    def __init__(self, filename, mode="single", ncam=False):

        # Global parameters
        self.day2sec = 86400.
        self.prefix  = 'cluster'

        # Flags
        self.plot = args.plot
        self.fcam = args.fcam
        
        # Output directory
        if args.outdir:
            self.odir = Path(args.outdir).resolve()
        elif args.project:
            self.project = args.project
            self.odir = Path(os.getenv('PLATO_WORKDIR')) / self.project / "input"
        else:
            self.odir = False

        # File names
        if self.odir:
            self.fileNameGap = f"{self.odir}/instrumentGap.txt"
            self.fileNameCCD = f"{self.odir}/instrumentCCD.txt"
            self.fileNamePRE = f"{self.odir}/instrumentPRE.txt"
            self.fileNameAPE = f"{self.odir}/instrumentAPE.txt"
            self.fileNameTED = f"{self.odir}/instrumentTED.txt"
            self.fileNameACS = f"{self.odir}/instrumentACS.txt"
        else:
            self.fileNameGap = self.fileNameCCD = self.fileNamePRE = self.fileNameAPE = self.fileNameACS = self.fileNameTED = False

        # Number of images in a quarter
        self.nimg = (60 * 60 * 24 * 90) / 25.

        # Select pointng field
        self.field = args.field
        self.ra, self.dec, self.kappa = getPointingField(self.field)
        
        # Short-hand definitions 
        N = args.stars
        G = args.group
        C = args.camera
        Q = args.quarter

        # Default ranges
        N = range(1, int(N)+1)
        if G is None: G = range(1, 4+1)
        if C is None: C = range(1, 6+1)
        if Q is None: Q = range(1, 8+1)

        # When number is given (also check for obvious mistakes)
        if len(G) == 1:
            if int(G) > 4:
                errorcode('warning', 'Maximum number of camera groups is 4!')
                parser.print_help()
            else:
                G = [int(G)]

        if len(C) == 1:
            if int(C) > 6:
                errorcode('warning', 'Maximum number of cameras in a group is 6!')
                parser.print_help()
            C = [int(C)]

        if len(Q) == 1:
            Q = [int(Q)]

        # When a range is given (also check for obvious mistakes)
        if len(G) == 3:
            if int(G[2]) > 4:
                errorcode('warning', 'Maximum number of camera groups is 4!')
                parser.print_help()
            else:
                G = range(int(G[0]), int(G[2])+1)

        if len(C) == 3:
            if int(C[2]) > 6:
                errorcode('warning', 'Maximum number of cameras in a group is 6!')
                parser.print_help()
            else:
                C = range(int(C[0]), int(C[2])+1)

        if isinstance(Q, str): 
            Q = convertQuarterRange(Q)
            Q = range(int(Q[0]), int(Q[1])+1)

        # Store parameters
        self.N = N
        self.G = G
        self.C = C
        self.Q = Q

        # Time column
        t0 = round(90. * (self.Q[0]  - 1) * self.day2sec)
        t1 = round(90. * (self.Q[-1])     * self.day2sec)
        self.time = np.arange(t0, t1, 25)




                        
    def createInputYAML(self):

        """Function to copy and adjust a yaml ready to launch.

        Parameters
        ----------
        field : str
            Observational PLATO field (e.g. SPF, NPF, LOPS2, LOPN1)
        odir : str, pathlib object
            Absolute output directory
        """

        # Get files names of YAML files
        yaml_old = Path(os.getenv("PLATO_PROJECT_HOME") + "/inputfiles/inputfile.yaml")
        yaml_new = self.odir / "inputfile.yaml"

        # Copy YAML if it doesn't exist already
        if not yaml_new.is_file():

            #print(f"Copying YAML configuration file : {yaml_new}")
            shutil.copy(yaml_old, yaml_new)

            # Find and replace a few strings:
            with open(yaml_new, 'r') as file:
                filedata = file.read()
                filedata = filedata.replace('inputfiles/starcatalog.txt', self.field)
                # Photon flux of a P=0 G2V-star [phot/s/m^2/nm]
                filedata = filedata.replace('1.00179e8       #', '0.73244782244e8 #')
                filedata = filedata.replace( 'NumColumns:                      100',
                                            f'NumColumns:                      7  ')
                filedata = filedata.replace( 'NumRows:                         100',
                                            f'NumRows:                         7  ')
                filedata = filedata.replace('IncludePhotometry:               no ',
                                            'IncludePhotometry:               yes')
                filedata = filedata.replace('MaskUpdateInterval:              14.0',
                                            'MaskUpdateInterval:              30.0')
                filedata = filedata.replace('GroupByExposure:                 yes',
                                            'GroupByExposure:                 no ')
                filedata = filedata.replace('WriteBiasMaps:                   yes',
                                            'WriteBiasMaps:                   no ')
                filedata = filedata.replace('WriteSmearingMaps:               yes',
                                            'WriteSmearingMaps:               no ')
                filedata = filedata.replace('WriteFlatfieldMap:               yes',
                                            'WriteFlatfieldMap:               no ')
                filedata = filedata.replace('WriteThroughputMaps:             yes',
                                            'WriteThroughputMaps:             no ')
                filedata = filedata.replace('WriteTransmissionEfficiency:     yes',
                                            'WriteTransmissionEfficiency:     no ')
                filedata = filedata.replace('WriteBackgroundMap:              yes',
                                            'WriteBackgroundMap:              no ')
                filedata = filedata.replace('WriteCTI:                        yes',
                                            'WriteCTI:                        no ')
                filedata = filedata.replace('WriteACS:                        yes',
                                            'WriteACS:                        no ')
                filedata = filedata.replace('WriteTelescopeACS:               yes',
                                            'WriteTelescopeACS:               no ')
                filedata = filedata.replace('WriteStarCatalog:                yes',
                                            'WriteStarCatalog:                no ')
                filedata = filedata.replace('WriteStarPositions:              yes',
                                            'WriteStarPositions:              no ')
                filedata = filedata.replace('WriteGhostPositions:             yes',
                                            'WriteGhostPositions:             no ')
                filedata = filedata.replace('WriteCosmics:                    yes',
                                            'WriteCosmics:                    no ')
                # Write the file out again
                if self.odir:
                    with open(yaml_new, 'w') as file:
                        file.write(filedata)




        
    def createJobScript(self):

        """Function to create a job script to be used on the VSC.
        """

        errorcode('module', f'\nJob script for parallisation\n')

        # Max nodes used
        maxNodes = 1

        # Number of node-cores for cluster
        coresCluster = 36

        # Memory per node for cluster [gb]
        pmemCluster = 192 - 8

        # Memory needed for 1 image [kb]
        memoryPerImage = 4.25

        # Run time to produce 1 image [s]
        timePerImage = 1.1 / 1000.

        # Overestimation of calculated memory and walltime [%]
        limitSafe = 5
        
        # Number of idividual simulations
        numSim = self.N[-1] * self.G[-1] * self.C[-1] * self.Q[-1]

        # Estimate number of nodes
        nodes = int(math.ceil(numSim/float(coresCluster)))

        # Estimate number of node-cores
        if nodes == 1:
            cores = coresCluster
        elif nodes > maxNodes:
            nodes = maxNodes
            cores = int(coresCluster * nodes)
        else:
            cores = int(math.ceil(numSim/nodes))

        # Estimate memory per core-node [mb]
        pmemSim = memoryPerImage * self.nimg / (nodes * cores) * 1e-3
        pmem    = int(math.ceil(pmemSim + limitSafe/100 * pmemSim))
        pmem = 50
        
        # Check if simulation quarter is larger than 5gb per node-core
        if pmemSim > 5000:
            errorcode('warning', 'Memory per node-core exceedes 5 GB!')

        # Check if if the total memory exceedes that of the cluster [gb]
        memPerNode = coresCluster * pmemSim*1e-3
        if memPerNode > pmemCluster-8:
            errorcode('error', 'Simulation node-memory of {int(memPerNode)} GB ' +
                      'exceedes the cluster node-memory of {pmemCluster}-8 GB!')

        # Estimate wall time given the resources above:
        walltimeQuarter = timePerImage * self.nimg
        walltimeSafe    = walltimeQuarter + limitSafe/100 * walltimeQuarter
        walltime        = str(datetime.timedelta(seconds=int(walltimeSafe)))
            
        # Print if simulation will run and if so, with what resources
        print(f'Preparing a SLURM job script for {numSim} individual simulations:')
        # print(f'Preparing {numSim} simulations with resources: ' +
        #       f'nodes={nodes}, cores={cores}, pmem={pmem}mb, walltime={walltime}')

        # Generate job script

        if self.odir:
            print(f"Creating parameterization job-script: {self.odir}/{self.prefix}.slurm")

        script = \
        f"""
        #!/bin/bash

        #SBATCH -A <account>       # Account name of cluster (mandatory)
        #SBATCH -o <stdout>        # Name of standard output
        #SBATCH -e <stderr>        # Name of standard errors
        #SBATCH -M <cluster>       # Cluster name (VSC: genius or wice)
        #SBATCH -p <software>      # Cluster software (VSC: batch_debug, batch, bigmen, etc).
        #SBATCH -N 1               # Number of nodes
        #SBATCH -n 36              # Number of CPU per node
        #SBATCH --mem-per-cpu=1G   # Amount of RAM memory per CPU
        #SBATCH -t 00:15:00        # Totoal execution of slurm job 

        cd $SLURM_SUBMIT_DIR
        module purge
        module restore <name>      # Name of cluster env with modules 
        
        # User defined parameters
        WORKDIR=<workdir_name>     # PLATOnium working directory
        PROJECT=<project_name>     # PLATOnium project name

        # Export paths
        export PLATO=$VSC_DATA/plato
        export PLATO_PROJECT_HOME=$VSC_DATA/PlatoSim3
        export PLATO_WORKDIR=$VSC_DATA/$WORKDIR
        export TEMDIR=$VSC_SCRATCH_NODE
        export OUTDIR=$VSC_SCRATCH/platosim/$PROJECT
        export PYTHONPATH=$PLATO:$PLATO_PROJECT_HOME/python
        export PLATONIUM=$PLATO_PROJECT_HOME/python/platosim/platonium/platonium
        export CONDA=/data/leuven/341/vsc34166/miniconda3/etc/profile.d/conda.sh

        # Activate environment
        source $CONDA
        conda activate platonium

        # Run PLATOnium
        python $PLATONIUM $star $group $camera $quarter --project $PROJECT -o $TEMDIR -d $OUTDIR --compress
        """

        # Save textfile for worker
        if self.odir:
            filename = f"{self.odir}/{self.prefix}.slurm"
            with open(filename, "w") as ofile:
                ofile.write(dedent(script).strip())




            
    def createParamFile(self):

        """Function to create a job script to be used on the VSC.
        """
        
        if self.odir:
            print(f"Creating HPC parameterization file  : {self.odir}/{self.prefix}.data")

        # Check if F-CAM is requested -> Group 5
        if self.fcam:
            self.G = range(5,6)
            self.C = range(1,3)
            
        # Add rows in a loop
        textfile = "star,group,camera,quarter\n"
        for starNo in self.N:
            for groupNo in self.G:
                for cameraNo in self.C:
                    for quarterNo in self.Q:
                        textfile += f"{starNo},{groupNo},{cameraNo},{quarterNo}\n"

        # Save textfile for worker
        if self.odir:
            filename = f"{self.odir}/{self.prefix}.data"
            with open(filename, "w") as ofile:
                ofile.write(dedent(textfile).strip())




                        
    def createDataGapsAndTransients(self):

        """Function to create a time columns including Data Gaps.
        Used to include data gaps relavant for space mission orbiting in L2.
        """

        # Generate a data-gap file
        errorcode('module', '\nDate gaps & Thermal transients\n')
        print('Downtime due to quater interuptions')
        print('Downtime due to montly data downlinks')
        print('Downtime due to loss of fine guidance')
        print('Downtime due to safe-mode events')
        _, self.t0, self.td = ns.getDataGaps(self.time, self.Q,
                                             outfile=self.fileNameGap, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameGap}")

        # Generate thermal transient file
        ns.temperatureTransients(self.time, self.t0, self.td, tempCCD=203.5,
                                 outfile=self.fileNameCCD, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameCCD}")





        
    def createPRE(self):

        """Function to create a Pointing Repeatability Error (PRE) file.
        Used to realistically include errors for each pointing.
        """

        # Generete PRE file
        errorcode('module', '\nPointing repeatability error (PRE)')
        ns.getPRE(self.ra, self.dec, self.kappa, self.Q, sigma=3,
                  outfile=self.fileNamePRE, show_table=True, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNamePRE}")

        


        
    def createAPE(self):

        """Function to create a Absolute Pointing Error (APE) file.
        Used to realistically include camera misalignments errors.
        """

        # Generete APE file
        errorcode('module', '\nAbsolute Pointing Error (APE)')
        ns.getAPE(self.ra, self.dec, self.kappa, sigma=3, fcam=self.fcam,
                  outfile=self.fileNameAPE, show_table=True, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameAPE}")




        
    def createTED(self):

        """Function to create a Thermo-Elastic Distortion (TED) file.
        """
        
        # Generate TED file
        errorcode('module', '\nThermo-Elastic Distortion (TED)')
        ns.getTED(self.Q, outfile=self.fileNameTED, show_table=True, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameTED}")




        
    def createACS(self):

        """Function to create a Attitude and orbit Control System (ACS) jitter file. 
        """
        
        # Generate TED file
        errorcode('module', '\nAttitude Control System (ACS)\n')
        ns.getACS(self.time, outfile=self.fileNameACS, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameACS}")



        
#--------------------------------------------------------------#
#                PARSING COMMAND-LINE ARGUMENTS                #
#--------------------------------------------------------------#

software = '\nInstrumental Noise Simulator'
parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', software))

man_group = parser.add_argument_group('MANDATORY PARAMETERS')
man_group.add_argument('stars',   type=str, help='Number of stars')
man_group.add_argument('field',   type=str, help='LOP (SPF, NPF)')

obs_group = parser.add_argument_group('OBSERVATION PARAMETERS')
obs_group.add_argument('--group',   metavar='NO.', type=str, help='Group   no.: 1, 2, .. (Default: 1-4 = all)')
obs_group.add_argument('--camera',  metavar='NO.', type=str, help='Camera  no.: 1, 2, .. (Default: 1-6 = all)')
obs_group.add_argument('--quarter', metavar='NO.', type=str, help='Quarter no.: 1, 2, .. (Default: 1-8 = 2yr)')
obs_group.add_argument('--fcam',    action='store_true', help='Setup for F-CAMs')

out_group = parser.add_argument_group('I/O PARAMETERS')
out_group.add_argument('-p', '--plot',   action='store_true',      help='Flag to plot each action')
out_group.add_argument('-o', '--outdir', metavar='PATH', type=str, help='Output directory to save')
out_group.add_argument('--project',      metavar='NAME', type=str, help='Name of PLATOnium project')

# Initialize instance of class
args = parser.parse_args()
x = Payload(args)

# Run each module
x.createInputYAML()
x.createJobScript()
x.createParamFile()
x.createDataGapsAndTransients()
x.createPRE()
x.createAPE()
x.createTED()
x.createACS()
print('')
