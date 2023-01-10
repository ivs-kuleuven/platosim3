#!/usr/bin/env python3

"""
Copy this script to your prefered working directory and change the
'User defined parameters' given below. By running this script a new
folder, with the name specified using 'simPrefix', will be created
along with a job script and textfile for each simulation configuration.
It is important that the 'memoryPerImage' and timePerIMage' parameters
are given sufficiently accurately as an underestimation may terminate
the job and a over estimation might leave the job in the que forever.
Retrieve these two parameters from a smaller simulation with the same
simulation setup on your local machine, and normalize to get the
parameters per image.

For each input argument select either
- A number  : e.g. "1"
- A range   : e.g. "1-2"
- Full range: i.e. "0"     (default quarterNo: 8Q = 2 years)
"""

import os
import math
import datetime
import argparse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt 
from pathlib import Path
from textwrap import dedent

import platosim.utilities  as ut
import platosim.instrument as it

#==============================================================#
#                         BEGIN CLASS                          #
#==============================================================#

class InstrumentSetup(object):

    """Class to generate a realistic setup for multi-camera simulations.
    
    This class makes it straight forward to generate a realistic random
    model for various systematic instrumental effects. It likewise add
    all the necessary files needed to launch the simulations on a high
    performance computing cluster

    Usage example for a single feather file:

    """

    def __init__(self, filename, mode="single", ncam=False):

        # Constants
        self.day2sec = 86400.

        # Output directory
        if args.outdir:
            self.odir = Path(args.outdir)
        elif args.project:
            self.project = args.project
            self.odir = Path(os.getenv('PLATO_WORKDIR')) / self.project / "input"
        else:
            self.odir = False
            
        # Number of images in a quarter
        self.nimg = (60 * 60 * 24 * 90) / 25

        # Provisional Long-duration Observation Phases (LOP):
        PF = {'SPF': [ 86.79870508, -46.39594703, 0.0],  # Galactic [253.0, -30.0]
              'NPF': [265.08002279,  39.5836954,  0,0]}  # Galactic [65.0,  +30.0]
        # Select pointng field
        self.ra, self.dec, self.kappa = PF[args.field]
        
        # Short-hand definitions 
        N = args.stars
        G = args.group
        C = args.camera
        Q = args.quarter

        # Default ranges
        N = range(1, int(N)+1)
        if G == '0': G = range(1, 4+1)
        if C == '0': C = range(1, 6+1)
        if Q == '0': Q = range(1, 8+1)

        # When number is given (also check for obvious mistakes)
        if len(G) == 1:
            if int(G) > 4:
                ut.errorcode('warning', 'Maximum number of camera groups is 4!')
                parser.print_help()
            else:
                G = [int(G)]

        if len(C) == 1:
            if int(C) > 6:
                ut.errorcode('warning', 'Maximum number of cameras in a group is 6!')
                parser.print_help()
            C = [int(C)]

        if len(Q) == 1:
            Q = [int(Q)]

        # When a range is given (also check for obvious mistakes)
        if len(G) == 3:
            if int(G[2]) > 4:
                ut.errorcode('warning', 'Maximum number of camera groups is 4!')
                parser.print_help()
            else:
                G = range(int(G[0]), int(G[2])+1)

        if len(C) == 3:
            if int(C[2]) > 6:
                ut.errorcode('warning', 'Maximum number of cameras in a group is 6!')
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
        t0 = round(90. * (self.Q[0]  - 1) * self.day2sec)
        t1 = round(90. * (self.Q[-1] - 1) * self.day2sec)
        self.time = np.arange(t0, t1, 25)
        
        

    def create_job_script(self):
        """
        Function to create a job script to be used on the VSC.
        """

        # prefix
        prefix = "run"
        
        # Name of yaml inputfile
        inputfile = "inputfile.yaml"

        # VSC project account details
        account = "default_project"

        # Email to recieve job notifications (use "None" to deactivate notifications)
        email = "nicholas.jannsen@kuleuven.be"

        # VSC module environment for PlatoSim3
        environment = "plato"

        # Max nodes used
        maxNodes = 10

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

        # Check if simulation quarter is larger than 5gb per node-core
        if pmemSim > 5000:
            ut.errorcode('warning', 'Memory per node-core exceedes 5 GB!')

        # Check if if the total memory exceedes that of the cluster [gb]
        memPerNode = coresCluster * pmemSim*1e-3
        if memPerNode > pmemCluster-8:
            ut.errorcode('error', "Simulation node-memory of {} GB exceedes the cluster node-memory of {}-8 GB!".format(int(memPerNode), pmemCluster))

        # Estimate wall time given the resources above:
        walltimeQuarter = timePerImage * self.nimg
        walltimeSafe    = walltimeQuarter + limitSafe/100 * walltimeQuarter
        walltime        = str(datetime.timedelta(seconds=int(walltimeSafe)))

        # Print if simulation will run and if so, with what resources
        ut.errorcode('message', f"\nPrepared {numSim} simulations with resources: nodes={nodes}, cores={cores}, pmem={pmem}mb, walltime={walltime}")

        # Check email argument
        if email is None: email = ''

        # Generate job script
        
        print(f"Creating HPC parameterization job-script: {self.odir}/{prefix}.pbs")

        script = \
        """
        #!/bin/bash

        #PBS -A lp_mesa_modeling
        #PBS -N output
        #PBS -l nodes=10:ppn=36
        #PBS -l pmem=5gb
        #PBS -l walltime=100:00:00

        cd $PBS_O_WORKDIR
        module purge
        module restore plato

        # Define sample
        PX=P1
        PROJECT=kul21
        OUTDIR=/scratch/leuven/341/vsc34166/platonium/$PROJECT

        # Load paths
        POETRY=$VSC_DATA/poetry/virtualenvs/platonium-KLYW-dHg-py3.8/bin/activate
        PLATO=$VSC_DATA/plato
        PLATO_PROJECT_HOME=$VSC_DATA/PlatoSim3
        PLATO_WORKDIR=$VSC_DATA/workdir
        SIMDIR=$PLATO_PROJECT_HOME/python/platosim/platonium
        TEMDIR=$VSC_SCRATCH_NODE
        PYTHONPATH=$POETRY:$PLATO:$PLATO_PROJECT_HOME/python:$PLATO_WORKDIR:$SIMDIR:$TEMDIR:$OUTDIR
        export POETRY
        export PLATO
        export PLATO_PROJECT_HOME
        export PLATO_WORKDIR
        export SIMDIR
        export TEMDIR
        export OUTDIR
        export PYTHONPATH

        # Load python environment
        source $POETRY

        # Run PLATOnium
        starID=$(printf "%09d" $star)
        python $SIMDIR/platonium $star $group $camera $quarter --project $PROJECT --sample $PX -o $TEMDIR --compress

        # Move data from the scrath node to output
        mkdir -p $OUTDIR/$starID
        zipfile=$TEMDIR/reduced/$PX/$starID/${starID}_Ncam${group}.${camera}_Q${quarter}.zip
        if [ -f $zipfile ]; then
            mv $zipfile $OUTDIR/$starID
        fi

        # Check if directory is not removed
        olddir=$TEMDIR/$PX/$starID/Ncam${group}${camera}
        if [ -d $olddir ]; then
            rm -rf $$olddir
        fi
        """

        # Save textfile for worker
        if self.odir:
            filename = f"{self.odir}/{prefix}.pbs"
            with open(filename, "w") as ofile:
                ofile.write(dedent(script).strip())




            
    def create_param_file(self):
        """
        Function to create a job script to be used on the VSC.
        """

        # Select prefix
        prefix = "data"
        
        print(f"Creating HPC parameterization file : {self.odir}/{prefix}.txt")

        # Add rows in a loop
        textfile = "star,group,camera,quarter\n"
        for starNo in self.N:
            for groupNo in self.G:
                for cameraNo in self.C:
                    for quarterNo in self.Q:
                        textfile += f"{starNo},{groupNo},{cameraNo},{quarterNo}\n"

        # Save textfile for worker
        if self.odir:
            filename = f"{self.odir}/{prefix}.txt"
            with open(filename, "w") as ofile:
                ofile.write(dedent(textfile).strip())




    def create_time_gaps(self):
        """
        Function to create a job script to be used on the VSC.
        """

        # Initialise random generator
        rng = np.random.default_rng()

        # Create pandas data frame for different flags
        flags = ["time", "roll", "downlink", "safemode", "fineguidance",
                 "attitude", "ccdanomaly", "argabrightning"]
        df = pd.DataFrame(columns=flags)

        # Start time of simulation
        t0 = round(90. * (self.Q[0]  - 1) * self.day2sec)
        t1 = round(90. * (self.Q[-1] - 1) * self.day2sec)
        time = np.arange(t0, t1, 25)
        # Create continous time array
        df["time"] = time

        # Quarterly rolls
        for Q in self.Q[1:]:                                    # [d]
            roll_gap    = 2 + np.random.uniform(-0.5, 0.5)      # [d]
            roll_event0 = (90. * (Q - 1) - roll_gap) * self.day2sec  # [s]
            roll_event1 = 90. * (Q - 1) * self.day2sec               # [s]
            df["roll"]  = np.logical_and(time>=roll_event0, time<=roll_event1)

        # Downlink gaps
        downlink_events = np.linspace(self.Q[1], self.Q[-1], 3*(self.Q[-1]-self.Q[1]))
        for M in downlink_events:
            downlink_gap    = (5/24. + np.random.uniform(-0.5/24., 0.5/24))
            downlink_event0 = (M - downlink_gap) * self.day2sec
            downlink_event1 = M * self.day2sec
            df["downlink"]  = np.logical_and(time>=downlink_event0, time<=downlink_event1)


        # Create new time array with all gaps
        df["tsim"] = ~df.roll * ~df.downlink
        df["tsim"] = df["tsim"].where(df["tsim"]!='T', 1) * df["time"]
        df["tsim"].replace({0: np.nan}, inplace=True) 
        df["tsim"].iloc[0] = 0

        df["tgap"] = ~df.roll * ~df.downlink
        df["tgap"] = ~df["tgap"]
        df["tgap"] = df["tgap"].where(df["tgap"]!='T', 1) * df["time"]

        df["tgap"].replace({0: np.nan}, inplace=True) 

        plt.plot(df["tsim"] / self.day2sec, np.ones(df.tsim.shape[0]), "b.")
        plt.plot(df["tgap"] / self.day2sec, np.ones(df.tsim.shape[0]), "r.")
        plt.show()
        
        print(df); exit()

        # Apply start time relative mission BOL
        self.beginExposureNr = round(self.timeStart / self.cadence)
        
        # Select prefix
        prefix = "instrument"
        
        print(f"Creating file with instrumental time gaps : {self.odir}/{prefix}.txt")
        




    def createPRE(self):
        """
        Function to create a Pointing Repeatability Error (PRE) file.
        Used to realistically include errors for each pointing.
        """

        # Generete PRE file
        if self.odir: print("\nSaving PRE file")
        it.getPRE(self.ra, self.dec, self.kappa, self.Q, sigma=3,
                  outdir=self.odir + "/instrumentPRE.txt", show_table=True)

        

    def createAPE(self):
        """
        Function to create a Absolute Pointing Error (APE) file.
        Used to realistically include camera misalignments errors.
        """

        # Generete APE file
        if self.odir: print("\nSaving APE file")
        it.getAPE(self.ra, self.dec, self.kappa, sigma=3,
                  outdir=self.odir + "/instrumentAPE.txt", show_table=True)



    def createTED(self):
        """
        Function to create a Absolute Pointing Error (APE) file
        """
        
        # Generate TED file
        if self.odir: print("\nSaving TED file and plot")
        it.getTED(self.Q, outdir=self.odir + "/instrumentTED.txt", plot=True)
        
        
#--------------------------------------------------------------#
#                PARSING COMMAND-LINE ARGUMENTS                #
#--------------------------------------------------------------#

software = '\nInstrument file generator'
parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=ut.errorcode('software', software))

man_group = parser.add_argument_group('MANDATORY PARAMETERS')
man_group.add_argument('field',   type=str, help='Number of stars')
man_group.add_argument('stars',   type=str, help='Number of stars')
man_group.add_argument('group',   type=str, help='Either: 1, 2, 3, 4')
man_group.add_argument('camera',  type=str, help='Either: 1, 2, 3, 4, 5, 6')
man_group.add_argument('quarter', type=str, help='Either: 1, 2, .. (Default Q1-Q8 -> 2 years)')

out_group = parser.add_argument_group('I/O PARAMETERS')
out_group.add_argument('-o', '--outdir',  metavar='PATH', type=str, help='Output directory to save')
out_group.add_argument('-p', '--project', metavar='NAME', type=str, help='Name of PLATOnium project')

# Initialize instance of class
args = parser.parse_args()
x = InstrumentSetup(args)

# Run each module
x.create_job_script()
x.create_param_file()
#x.create_time_gaps()
x.createPRE()
x.createAPE()
x.createTED()
