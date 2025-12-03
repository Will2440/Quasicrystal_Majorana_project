#!/bin/bash
  
#SBATCH --job-name=QC-PT
#SBATCH --nodes=1
#SBATCH --partition=compute
#SBATCH --ntasks-per-node=1
#SBATCH --account=phys035904
#SBATCH --array=1-2
#SBATCH --time=00:10:00
#SBATCH --mem=500M

module add languages/julia
julia /user/home/hb21877/Quasicrystal_Majorana_project/simulations/data_collection/hpc_compatible/serial_batch/main.jl $SLURM_ARRAY_TASK_ID