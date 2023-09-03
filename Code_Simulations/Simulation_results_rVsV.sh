#!/usr/bin/env bash
## name of jobb
#SBATCH -J PREVAC_2000Simu_7500Sub_rVsV

## Resources: (nodes, procs, tasks, walltime, ... etc)
#SBATCH -N 1 -C zonda|bora|miriel

#SBATCH -t 5:00:00
# #  standard output message
#SBATCH -o batch_run_files/batch%A_%a.out
# # output error message
#SBATCH -e batch_run_files/batch%A_%a.err
#SBATCH --exclusive
#SBATCH --array=1-200%15

module purge

i="${SLURM_ARRAY_TASK_ID}"
echo $i

module load tools/R/3.6.2

# Model estimation ---
Rscript --vanilla Rcode/Simulation_rVsV_2.R "7500" $i