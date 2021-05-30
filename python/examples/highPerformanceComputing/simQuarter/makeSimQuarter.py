#!/usr/bin/env python3

"""
Copy this script to your prefered working directory and change the 'User defined parameters' given below. By running this script a new folder, with the name specified using 'simPrefix', will be created along with a job script and textfile for each simulation configuration. It is important that the 'memoryPerImage' and timePerIMage' parameters are given sufficiently accurately as an underestimation may terminate the job and a over estimation might leave the job in the que forever. Retrieve these two parameters from a smaller simulation with the same simulation setup on your local machine, and normalize to get the parameters per image.
"""

import os
import sys
import math
import datetime
from textwrap import dedent
from colorama import Fore, Style, Back

#--------------------------------------------------------------#
#                    USER DEFINED PARAMEERS                    #
#--------------------------------------------------------------#

# Name of simulation
simPrefix = "test0"

# Name of yaml inputfile
inputfile = "inputfile.yaml"

# VSC project account details
account = "default_project"

# Email to recieve job notifications (use "None" to deactivate notifications)
email = "nicholas.jannsen@kuleuven.be"

# VSC module environment for PlatoSim3
environment = "env-PlatoSim3"

# Number of node-cores for cluster
coresCluster = 36

# Memory per node for cluster [gb]
pmemCluster  = 192

# Memory needed for 1 image [kb]
memoryPerImage = 50

# Run time to produce 1 image [s]
timePerImage = 1

# Overestimation of calculated memory and walltime [%]
limitSafe = 2

#--------------------------------------------------------------#
#                PARSING COMMAND-LINE ARGUMENTS                #
#--------------------------------------------------------------#

# Help usage function:
usage = \
"""
Usage: {0} <cameraGroupNo> <cameraNo> <quarterNo>

cameraGroupNo : [1, 2, 3, 4]        (Maximum 4 camera groups)
cameraNo      : [1, 2, 3, 4, 5, 6]  (Maximum 6 cameras per group)
quarterNo     : [1, 2, 3, 4, .. n]  (No maximum, hence use logic)

For each input argument select either
- A number  : e.g. "1"
- A range   : e.g. "1-2"
- Full range: i.e. "0"         (default quarterNo is 8 = 2 years)

Example: python {0} 1 0 1-4

Before running this script all the 'user defined parameters' needs
to be modified to ensure that the automatic estimation of needed
nodes, cores, and memory are rigthfully calculated.
""".format(sys.argv[0])

def error(): print(Fore.RED + Style.BRIGHT + '[ERROR]: Wrong input!' + Style.RESET_ALL)
def help(): print(Fore.BLUE + Style.BRIGHT + usage + Style.RESET_ALL); exit()

# Print usage if no arguments are given:
if len(sys.argv) == 1: help()

# Check for opvious mistakes
if len(sys.argv) < 4 or len(sys.argv) > 4:
    error(); help(); exit(1)

# Simplify expression
G = sys.argv[1]
C = sys.argv[2]
Q = sys.argv[3]

# Default ranges
if G=='0': G = range(1, 4+1)
if C=='0': C = range(1, 6+1)
if Q=='0': Q = range(1, 8+1)

# When number is given (also check for obvious mistakes)
if len(G)==1:
    if int(G) > 4: error(); help()
    else: G = [int(G)]
if len(C)==1:
    if int(C) > 6: error(); help()
    C = [int(C)]
if len(Q)==1: Q = [int(Q)]

# When a range is given (also check for obvious mistakes)
if len(G)==3:
    if int(G[2]) > 4: error(); help()
    else: G = range(int(G[0]), int(G[2])+1)
if len(C)==3:
    if int(C[2]) > 6: error(); help()
    else: C = range(int(C[0]), int(C[2])+1)
if len(Q)==3: Q = range(int(Q[0]), int(Q[2])+1)

#--------------------------------------------------------------#
#                     RESOURCES AND PARAMETER                  #
#--------------------------------------------------------------#

# Number of images in a quarter
numImgQuarter = 60 * 60 * 24 * 90 / 25

# Number of idividual simulations
numSim = G[-1] * C[-1] * Q[-1]

# Estimate number of nodes
nodes = int(math.ceil(numSim/float(coresCluster)))

# Estimate number of node-cores
if nodes == 1: cores = coresCluster
if nodes > 1:  cores = int(math.ceil(numSim/nodes))

# Estimate memory per core-node [mb]
pmemSim = memoryPerImage * numImgQuarter / (nodes * cores) * 1e-3
pmem    = int(round(pmemSim + limitSafe/100 * pmemSim))

# Check if simulation quarter is larger than 5gb per node-core
if pmemSim > 5000:
    print(Fore.YELLOW + Style.BRIGHT + "[WARNING]: Memory per node-core exceedes 5 GB!" + Style.RESET_ALL)

# Check if if the total memory exceedes that of the cluster [gb]
memPerNode = coresCluster * pmemSim*1e-3
if memPerNode > pmemCluster-8:
    print(Fore.RED + Style.BRIGHT + "[ERROR]: Simulation node-memory of {} GB exceedes the cluster node-memory of {}-8 GB!".format(int(memPerNode), pmemCluster) + Style.RESET_ALL); exit()

# Estimate wall time given the resources above:
walltimeQuarter = timePerImage * numImgQuarter / (nodes * cores)
walltimeSafe    = walltimeQuarter + limitSafe/100 * walltimeQuarter
walltime        = str(datetime.timedelta(seconds=int(walltimeSafe)))

# Print if simulation will run and if so, with what resources
print(Style.BRIGHT + "Preparing {} simulations with resources: nodes={}, cores={}, pmem={}mb, walltime={}".format(numSim, nodes, cores, pmem, walltime) + Style.RESET_ALL)

# Check email argument
if email is None: email = ''

# Create output directory for data
outputDir = os.getcwd() + "/" + simPrefix
createDir = os.mkdir(outputDir)

#--------------------------------------------------------------#
#            CREATE TEXTFILE OF PARAMETER INSTANCES            #
#--------------------------------------------------------------#

print("Creating HPC parameterization text-file : " + outputDir + "/" + simPrefix + ".txt")

textfile = "group camera quarter \n"
for groupNo in G:
    for cameraNo in C:
        for quarterNo in Q:
            textfile += "{0} {1} {2} \n".format(groupNo, cameraNo, quarterNo)

# Save textfile for worker
filenameTextfile = outputDir + "/" + simPrefix + ".txt"
with open(filenameTextfile, "w") as outputFile:
    outputFile.write(dedent(textfile).strip())

#--------------------------------------------------------------#
#                       CREATE JOB SCRIPT                      #
#--------------------------------------------------------------#

print("Creating HPC parameterization job-script: " + outputDir + "/" + simPrefix + ".pbs")

script = \
"""
#!/bin/bash

#PBS -N {0}
#PBS -o stdout_$group_$camera$_quarter
#PBS -e stderr_$group_$camera_$quarter
#PBS -A {1}
#PBS -m a -M {2}

#PBS -l nodes={3}:ppn={4}
#PBS -l pmem={5}mb
#PBS -walltime {6}

cd $PBS_O_WORKDIR

module purge
module restore {7}

PLATO_PROJECT_HOME=$VSC_DATA/PlatoSim3
export PLATO_PROJECT_HOME
PYTHONPATH=$PLATO_PROJECT_HOME/python:$VSC_DATA/python_lib/lib/python3.7/site-packages/
export PYTHONPATH

python $PLATO_PROJECT_HOME/python/examples/runSimQuarter.py {8} $group $camera $quarter
""".format(simPrefix, account, email, nodes, cores, pmem, walltime, environment, inputfile)

# Save textfile for worker
filenameScript = outputDir + "/" + simPrefix + ".pbs"
with open(filenameScript, "w") as outputFile:
    outputFile.write(dedent(script).strip())
