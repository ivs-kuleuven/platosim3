
# Usage:
#        $ python3 makeSimQuarterSlurmScripts.py inputfile.yaml
#        $ cd folder_with_slurm_files
#        $ for f in `ls slurm*.sh`; do sbatch $f; done
#        $ sacct --format=jobid,jobname,account,user,partition,ntasks,alloccpus,elapsed,state


import os
import sys
from textwrap import dedent

inputfile = sys.argv[1]

simulationPrefix = "Run1"
slurmScriptOutputFolder = "/home/joris/slurm"

print("Creating slurm scripts " + slurmScriptOutputFolder + "/slurm_" + simulationPrefix + "_group*_camera*_Q*.sh")

platoDir = os.getenv("PLATO_PROJECT_HOME")

for groupNr in [1,2,3,4]:
    for cameraNr in [1,2,3,4,5,6]:
        for quarterNr in [1,2,3,4]:

            script = \
            """
            #!/bin/sh

            #SBATCH --job-name=platosimslurm
            #SBATCH --account=ivsusers
            #SBATCH --time 30
            #         (estimated run time in minutes for each of the jobs)
            #SBATCH --output=stdout_{0}_{1}_{2}.txt
            #SBATCH --error=stderr_{0}_{1}_{2}.txt
            #SBATCH --cpus-per-task=1
            #         (only change if multithreaded)
            #         (if OpenMP set export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK)
            #SBATCH --ntasks=1
            #         (default value)
            #         (if openmpi, use this value to execute a job on multiple nodes)
            #SBATCH --mem=4096
            #         (memory in MB)
            #SBATCH --partition=normal
            #         (partitions: high, normal, low, desktops, longjobs)

            python3 {4}/python/Examples/simQuarter.py {0} {1} {2} {3}
            """.format(inputfile, groupNr, cameraNr, quarterNr, platoDir)

            path = slurmScriptOutputFolder + "/slurm_" + simulationPrefix + "_group{0}_camera{1}_Q{2}.sh".format(groupNr, cameraNr, quarterNr)
            with open(path, "w") as outputFile:
               outputFile.write(dedent(script).strip())

