#!/bin/sh

# Run this script on the bash prompt with with e.g. 
#      $ sbatch --array=0-9 demo_slurmscript.sh
# to split the time series in 10 different parts that are computed simultaneously.

#SBATCH --job-name=platosimslurm
#SBATCH --account=ivsusers
#SBATCH --time 60
#         (estimated run time in minutes for each of the jobs)
#SBATCH --output=stdout_%A_%a.txt
#SBATCH --error=stderr_%A_%a.txt
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

python demo_slurm.py $SLURM_ARRAY_TASK_MIN $SLURM_ARRAY_TASK_MAX $SLURM_ARRAY_TASK_ID
