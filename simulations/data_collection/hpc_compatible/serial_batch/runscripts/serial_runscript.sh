#!/bin/bash

#SBATCH --job-name=Ph
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=phys030424
#SBATCH --array=1-10
#SBATCH --time=01:00:00
#SBATCH --mem=2G

module add languages/julia
JULIA_NUM_THREADS=1 julia /user/work/hb21877/Quasicrystal_Majorana_project/simulations/data_collection/hpc_compatible/serial_batch/main_phason.jl $SLURM_ARRAY_TASK_ID