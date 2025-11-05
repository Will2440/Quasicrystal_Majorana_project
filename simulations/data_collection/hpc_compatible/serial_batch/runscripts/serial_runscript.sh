#!/bin/bash
  
#SBATCH --job-name=QC-PT
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --ntasks-per-node=1
#SBATCH --account=phys033184
#SBATCH --array=1-1608
#SBATCH --time=00:40:00
#SBATCH --mem=500M


module add languages/julia
julia /user/home/hb21877/Quasicrystal_Majorana_project_clone/Quasicrystal_Majorana_project/simulations/data_collection/hpc_compatible/serial_batch/main.jl $SLURM_ARRAY_TASK_ID