#!/bin/bash
  
#SBATCH --job-name=GQC_np_trial
#SBATCH --nodes=1
#SBATCH --partition=test
#SBATCH --ntasks-per-node=10
#SBATCH --account=phys033184
#SBATCH --array=1-1
#SBATCH --time=00:10:00
#SBATCH --mem=100M


module add languages/julia
JULIA_NUM_THREADS=10 julia /user/home/hb21877/Majorana_solver_for_BlueCrystal/serial_compute/serial_run_main_N_loop_np_trial.jl $SLURM_ARRAY_TASK_ID
