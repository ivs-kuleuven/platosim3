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
- Full range: i.e. "0"         (default quarterNo is 8 = 2 years)

User examples:
   python makeScriptVSC.py 1 0 1-4

Before running this script all the 'user defined parameters' needs
to be modified to ensure that the automatic estimation of needed
nodes, cores, and memory are rigthfully calculated.

Author: Nicholas Jannsen (nicholas.jannsen@kuleuven.be)
"""

import os
import math
import datetime
import argparse
from textwrap import dedent
from platosim.utilities import errorcode, convertQuarterRange

#--------------------------------------------------------------#
#                    USER DEFINED PARAMEERS                    #
#--------------------------------------------------------------#

# Name of simulation
simPrefix = "run"

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
pmemCluster  = 192

# Memory needed for 1 image [kb]
memoryPerImage = 4.25

# Run time to produce 1 image [s]
timePerImage = 1.1 / 1000.

# Overestimation of calculated memory and walltime [%]
limitSafe = 5

#--------------------------------------------------------------#
#                PARSING COMMAND-LINE ARGUMENTS                #
#--------------------------------------------------------------#

software = '\nHPC script generator'
parser = argparse.ArgumentParser(epilog=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter,
                                 description=errorcode('software', software))

parser.add_argument('starNo',        type=str, help='Number of stars')
parser.add_argument('cameraGroupNo', type=str, help='Either: 1, 2, 3, 4')
parser.add_argument('cameraNo',      type=str, help='Either: 1, 2, 3, 4, 5, 6')
parser.add_argument('quarterNo',     type=str, help='Either: 1, 2, .. (Default Q1-Q8 -> 2 years)')
parser.add_argument('outDir',        type=str, help='Path to where the data shall be stored')

args = parser.parse_args()
N = args.starNo
G = args.cameraGroupNo
C = args.cameraNo
Q = args.quarterNo
outputdir = args.outDir

# Default ranges
N = range(1, int(N)+1)
if G == '0': G = range(1, 4+1)
if C == '0': C = range(1, 6+1)
if Q == '0': Q = range(1, 8+1)

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

#--------------------------------------------------------------#
#                     RESOURCES AND PARAMETER                  #
#--------------------------------------------------------------#

# Number of images in a quarter (2 days are lost due to rotation)
numImgQuarter = (60 * 60 * 24 * 90 - 2) / 25

# Number of idividual simulations
numSim = N[-1] * G[-1] * C[-1] * Q[-1]

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
pmemSim = memoryPerImage * numImgQuarter / (nodes * cores) * 1e-3
pmem    = int(math.ceil(pmemSim + limitSafe/100 * pmemSim))

# Check if simulation quarter is larger than 5gb per node-core
if pmemSim > 5000:
    errorcode('warning', 'Memory per node-core exceedes 5 GB!')

# Check if if the total memory exceedes that of the cluster [gb]
memPerNode = coresCluster * pmemSim*1e-3
if memPerNode > pmemCluster-8:
    errorcode('error', "Simulation node-memory of {} GB exceedes the cluster node-memory of {}-8 GB!".format(int(memPerNode), pmemCluster))

# Estimate wall time given the resources above:
walltimeQuarter = timePerImage * numImgQuarter
walltimeSafe    = walltimeQuarter + limitSafe/100 * walltimeQuarter
walltime        = str(datetime.timedelta(seconds=int(walltimeSafe)))

# Print if simulation will run and if so, with what resources
errorcode('message', f"Prepared {numSim} simulations with resources: nodes={nodes}, cores={cores}, pmem={pmem}mb, walltime={walltime}")

# Check email argument
if email is None: email = ''

# Create output directory for data
outputDir = outputdir

#--------------------------------------------------------------#
#            CREATE TEXTFILE OF PARAMETER INSTANCES            #
#--------------------------------------------------------------#

print("Creating HPC parameterization text-file : " + outputDir + "/" + simPrefix + ".txt")

textfile = "star,group,camera,quarter \n"
for starNo in N:
    for groupNo in G:
        for cameraNo in C:
            for quarterNo in Q:
                textfile += f"{starNo},{groupNo},{cameraNo},{quarterNo}\n"

# Save textfile for worker
filenameTextfile = outputDir + "/" + simPrefix + ".txt"
with open(filenameTextfile, "w") as outputFile:
    outputFile.write(dedent(textfile).strip())

#--------------------------------------------------------------#
#                       CREATE JOB SCRIPT                      # TODO!
#--------------------------------------------------------------#
exit()
print("Creating HPC parameterization job-script: " + outputDir + "/" + simPrefix + ".pbs")

script = \
"""
#!/bin/bash

#PBS -N {simPrefix}
#PBS -o stdout_$group_$camera_$quarter
#PBS -e stderr_$group_$camera_$quarter
#PBS -A {aaccount}
#PBS -m a -M {email}

#PBS -l nodes={nodes}:ppn={cores}
#PBS -l pmem={pmem}mb
#PBS -l walltime={walltime}

cd $PBS_O_WORKDIR

module purge
module restore {environment}

PLATO_PROJECT_HOME=$VSC_DATA/PlatoSim3
export PLATO_PROJECT_HOME
PYTHONPATH=$PLATO_PROJECT_HOME/python:$VSC_DATA/python_lib/lib/python3.7/site-packages/
export PYTHONPATH

python $PLATONIUM/platonium/platonium.py {simName} $star $group $camera $quarter
"""

# Save textfile for worker
filenameScript = outputDir + "/" + simPrefix + ".pbs"
with open(filenameScript, "w") as outputFile:
    outputFile.write(dedent(script).strip())

#--------------------------------------------------------------#
#                   CREATE COMPRESSION SCRIPT                  #
#--------------------------------------------------------------#

# print("Creating HPC tar-compression job-script: " + outputDir + "/" + simPrefix + ".pbs")

# script = \
# """
# #!/bin/bash -l

# #PBS -N FileCompression
# #PBS -o stdout_$group_$camera_$quarter
# #PBS -e stderr_$group_$camera_$quarter
# #PBS -A {1}
# #PBS -m a -M {2}

# #PBS -l nodes={3}:ppn={4}
# #PBS -l pmem={5}mb
# #PBS -l walltime={6}

# cd $PBS_O_WORKDIR

# module purge
# module restore {7}

# source $VSC_HOME/.bashrc

# PYTHONPATH=$PLATO_PROJECT_HOME/python:$VSC_DATA/python_lib/lib/python3.7/site-packages/
# export PYTHONPATH

# python createCompression.py {8}
# """.format(simPrefix, account, email, nodes, cores, pmem, walltime, environment, inputfile)

# # Save textfile for worker
# filenameScript = outputDir + '/' + simPrefix + 'Compression' + '.pbs'
# with open(filenameScript, 'w') as outputFile:
#     outputFile.write(dedent(script).strip())
