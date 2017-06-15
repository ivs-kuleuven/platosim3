
# Example script on how to use slurm. 
#
# Don't run this script directly with Python, but launch the script demo_slurmscript.sh
# which will in turn launch this python script. For example if you run on the prompt:
#      $ sbatch --array=0-9 demo_slurmscript.sh
# the time series will be into 10 different parts that are computed simultaneously,
# and written to 10 different HDF5 output files.


import os
import sys
import numpy as np
from math import ceil

from simfile import SimFile
from simulation import Simulation

# Slurm sets some environment variables for each job, which we can use. 
# The slurm script specifies an array of jobs, for which we can derive 
# the total number of jobs, and the sequential number of the current job.
# argv[1]: the task ID of the first task in the slurm array
# argv[2]: the task ID of the last task in the slurm array
# argv[3]: the task ID of the current task in the slurm array

Njobs = int(sys.argv[2]) - int(sys.argv[1]) + 1
jobNr = int(sys.argv[3])

# Specify the absolute paths of the input files and the output folder.

inputDir    = os.getenv("PLATO_PROJECT_HOME") + "/inputfiles"
inputFile   = inputDir + "/inputfile.yaml"

outputDir   = os.getcwd()
outputFile  = "SlurmDemoOutput{0:04d}".format(jobNr)

# Configure a Simulation object

sim = Simulation(outputFile, inputFile)
sim.outputDir = outputDir

sim["PSF/Model"] = "MappedGaussian"
sim["ObservingParameters/RApointing"]  = 180.0 
sim["ObservingParameters/DecPointing"] = -70.0
sim["Telescope/AzimuthAngle"]          =   0.0
sim["Telescope/TiltAngle"]             =   0.0

# Specify the number of exposures and the nr of the begin exposure for this job

totalNumberOfExposures = 105
numExposuresPerJob = int(ceil(float(totalNumberOfExposures) / Njobs))
beginExposureNrOfThisJob = jobNr * numExposuresPerJob

sim["ObservingParameters/BeginExposureNr"] = beginExposureNrOfThisJob

# Only for the very last job, we need to adapt the number of exposures because
# what's left after (Njobs-1) jobs may not be the full 'numExposuresPerJob'.

if jobNr == Njobs-1:
    numExposuresOfTheLastJob = totalNumberOfExposures - (Njobs-1) * numExposuresPerJob
    sim["ObservingParameters/NumExposures"] = numExposuresOfTheLastJob
else:
    sim["ObservingParameters/NumExposures"] = numExposuresPerJob

# Run the simulation

simFile = sim.run()
