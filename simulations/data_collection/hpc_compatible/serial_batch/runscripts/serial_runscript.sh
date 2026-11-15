#!/bin/bash

#SBATCH --job-name=MP_QC-SC_rerun
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --account=phys030424
#SBATCH --array=1-121
#SBATCH --time=03:00:00
#SBATCH --mem=5G

module add languages/julia
export JULIA_DEPOT_PATH="/user/work/hb21877/.julia"
JULIA_NUM_THREADS=1 julia /user/work/hb21877/Quasicrystal_Majorana_project/simulations/data_collection/hpc_compatible/serial_batch/main.jl $SLURM_ARRAY_TASK_ID