#!/bin/bash
  
#SBATCH --job-name=TMQC_mu_loop
#SBATCH --nodes=1
#SBATCH --partition=cpu
#SBATCH --ntasks-per-node=1
#SBATCH --account=phys033184
#SBATCH --array=1-1608
#SBATCH --time=00:40:00
#SBATCH --mem=500M


module add languages/julia
julia /user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/serial_run_main_mu_loop_TMQC.jl $SLURM_ARRAY_TASK_ID
