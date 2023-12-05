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
import shutil
import argparse
import datetime

# PlatoSim standard
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from pathlib import Path

# PlatoSim imports
import platosim.noise     as ns
import platosim.slurm     as sm
import platosim.utilities as ut
from platosim.utilities import errorcode


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
        self.aocs = args.aocs

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
            self.fileNamePRE = f"{self.odir}/instrumentPRE.txt"
            self.fileNameAPE = f"{self.odir}/instrumentAPE.txt"
            self.fileNameTED = f"{self.odir}/instrumentTED.txt"
            self.fileNameACS = f"{self.odir}/instrumentACS.txt"
            self.fileNameGap = f"{self.odir}/instrumentGap.ftr"
            self.fileNameCCD = f"{self.odir}/instrumentCCD.txt"
        else:
            self.fileNameGap = self.fileNameCCD = self.fileNamePRE = self.fileNameAPE = self.fileNameACS = self.fileNameTED = False

        # Number of images in a quarter
        self.nimg = round(ut.year() / 4 / 25)

        # Select pointng field
        self.field = args.field
        self.alpha, self.delta, self.kappa = ut.getPointingField(self.field)

        # Short-hand definitions 
        N = args.ids
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
            Q = ut.convertQuarterRange(Q)
            Q = range(int(Q[0]), int(Q[1])+1)

        # Store parameters
        self.N = N
        self.G = G
        self.C = C
        self.Q = Q

        # Time column
        tQ = ut.year() / self.day2sec / 4
        t0 = round(tQ * (self.Q[0]  - 1) * self.day2sec)
        t1 = round(tQ * (self.Q[-1])     * self.day2sec)
        self.time = np.arange(t0, t1, 25)




                        
    def createInputYAML(self):

        """Function to copy and adjust a yaml ready to launch.
        """

        if self.odir:
            ut.copyInputYAML(self.field, self.odir)


    


    def createParamFile(self):

        """Function to create a job script to be used on the VSC.
        """
        
        if self.odir:
            filename = f"{self.odir}/{self.prefix}_ncams.data"
            print(f"Creating HPC parameterization file  : {filename}")
            sm.getParamFile(self.N, self.G, self.C, self.Q,
                            fcam=False, ofile=filename)
            
            if self.fcam:
                filename = f"{self.odir}/{self.prefix}_fcams.data"
                print(f"Creating HPC parameterization file  : {filename}")
                sm.getParamFile(self.N, range(5,6), range(1,3), self.Q,
                                fcam=True, ofile=filename)




    def createJobScript(self):

        """Function to create a job script to be used on the VSC.
        """

        if self.odir:
            errorcode('module', f'\nJob script for parallisation\n')
            filename = f"{self.odir}/{self.prefix}.slurm"
            print(f"Creating parameterization job-script: {filename}")
            sm.getJobScript(self.N, self.G, self.C, self.Q,
                            ofile=filename)




                            
    def createPRE(self):

        """Function to create a Pointing Repeatability Error (PRE) file.
        Used to realistically include errors for each pointing.
        """

        # Generete PRE file
        errorcode('module', '\nPointing repeatability error (PRE)')
        ns.getPRE(self.alpha, self.delta, self.kappa, self.Q, sigma=3,
                  ofile=self.fileNamePRE, table=True, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNamePRE}")

        


        
    def createAPE(self):

        """Function to create a Absolute Pointing Error (APE) file.
        Used to realistically include camera misalignments errors.
        """

        # Generete APE file
        errorcode('module', '\nAbsolute Pointing Error (APE)')
        ns.getAPE(self.alpha, self.delta, self.kappa, sigma=3,
                  ofile=self.fileNameAPE, table=True, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameAPE}")


        


    def createGain(self):

        """Function to create a gain variation file.
        Used to realistically include gain variations for CCD and FEE.
        TODO under struction!
        """

        # Generete APE file
        errorcode('module', '\nCCD and FEE gain variations')
        #ns.getGain(gain0CCD=, gain0FEE=, sigma=3,
        #           ofile=self.fileNameAPE, table=True, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameAPE}")




        
    def createGap(self):

        """Function to create a time columns including Data Gaps.
        Used to include data gaps relavant for space mission orbiting in L2.
        """

        # Generate a data-gap file
        errorcode('module', '\nDate gaps & Thermal transients\n')
        print('Downtime due to quater interuptions')
        print('Downtime due to loss of fine guidance')
        print('Downtime due to safe-mode events')
        _, self.t0, self.td = ns.getDataGaps(self.time, self.Q,
                                             ofile=self.fileNameGap, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameGap}")

        # Generate thermal transient file
        ns.temperatureTransients(self.time, self.t0, self.td, tempCCD=203.5,
                                 ofile=self.fileNameCCD, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameCCD}")

        


        
    def createTED(self):

        """Function to create a Thermo-Elastic Distortion (TED) file.
        """
        
        # Generate TED file
        errorcode('module', '\nThermo-Elastic Distortion (TED)\n')
        ns.getTED(self.Q, ofile=self.fileNameTED, table=True, plot=self.plot)
        if self.odir: print(f"File saved: {self.fileNameTED}")




        
    def createACS(self):

        """Function to create a Attitude and orbit Control System (ACS) jitter file. 
        """
        
        # Generate ACS file
        if self.aocs and self.odir:
            errorcode('module', '\nAttitude Control System (ACS)\n')
            ns.getACS(self.time, ofile=self.fileNameACS, plot=self.plot)
            if self.odir: print(f"File saved: {self.fileNameACS}")




            
#--------------------------------------------------------------#
#                PARSING COMMAND-LINE ARGUMENTS                #
#--------------------------------------------------------------#

software = '\nInstrumental Noise Simulator'
parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', software))

man_group = parser.add_argument_group('MANDATORY PARAMETERS')
man_group.add_argument('ids',   type=str, help='Number of IDs (stars or CCDs=4)')
man_group.add_argument('field', type=str, help='PLATO LOP field [LOPS2, LOPN1, SPF, NPF]')

out_group = parser.add_argument_group('I/O PARAMETERS')
out_group.add_argument('-p', '--plot',    action='store_true',      help='Flag to plot each action')
out_group.add_argument('-v', '--verbose', metavar='INT',  type=int, help='Verbosity level {0, 1, 3} (Default: 1)')
out_group.add_argument('-o', '--outdir',  metavar='PATH', type=str, help='Output directory to save')
out_group.add_argument('--project',       metavar='NAME', type=str, help='Name of PLATOnium project')

obs_group = parser.add_argument_group('OBSERVATION PARAMETERS')
obs_group.add_argument('--group',   metavar='INT', type=str, help='Group   number: 1, 2, .... (Default: 1-4 = all)')
obs_group.add_argument('--camera',  metavar='INT', type=str, help='Camera  number: 1, 2, ... (Default: 1-6 = all)')
obs_group.add_argument('--quarter', metavar='INT', type=str, help='Quarter number: 1, 2, .. (Default: 1-8 = 2yr)')
obs_group.add_argument('--fcam',    action='store_true', help='Flag to generate files for the F-CAMs')
obs_group.add_argument('--aocs',    action='store_true', help='Flag to generate a red noise AOCS jitter file')

args = parser.parse_args()

#--------------------------------------------------------------#
#                            WORKFLOW                          #
#--------------------------------------------------------------#

# Start time tracking
tic = datetime.datetime.now()

# Initialize instance of class
x = Payload(args)

# Run each module
x.createInputYAML()
x.createParamFile()
x.createJobScript()
x.createPRE()
x.createAPE()
x.createGap()
x.createTED()
x.createACS()

# Finish with output
if (args.verbose is None) or (args.verbose > 0):
    toc = datetime.datetime.now()
    print(f'\nTotal execution time: {toc-tic} [hh:mm:ss]\n')
