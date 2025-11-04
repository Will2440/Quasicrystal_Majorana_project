# HPC-Compatible Serial Batch Jobs

## Overview
This is the unedited work environment from the final calculations of the masters thesis project using BlueCrystal. This is less user-oriented than the scripts in <local_machine/> as there is not the same UserOptions architecture. However, much of the code is comparible and this readme will explain how to execute tasks.

## Usage Notes

1) Decide on the parameter ranges to be solved for and compelete one of the two param_prep[...].jl scripts. You can check the .txt file containing the parameters you just created under <param_sets/> These two scripts are essentially the same but tailroed towards looping over mu with fixed N (param_prep_mu_loop.jl), looping over different N (param_prep_N_loop.jl) or looping over special ranges of N defined by sequence inflation lengths (param_prep_N_fractality.jl).
2) Enter one of the serial_run_main[...].jl scripts. Check 
   1) the sequnece generation is long enough for your chosen N_range;
   2) the param_filepath is correct;
   3) there are the correct number of t_columns for your chosen sequence;
   4) check the sequence is defined correctly to the one you want;
   5) check the precision, if relevant;
   6) check the data_save_filepath is correct and has been created by you beforehand;
   7) check the solving function you intend to use is correctly being called from hp_serial_solver.jl or np_serial_solver.jl.
3) You are ready to run the code. Go to <runscripts/> and open the correct bash file. Check all of the info, especially that the array matches the length of the aprticular parameter set file being used, that the correct serial_run_main[...].jl is being called, and the correct nubmer of threads is being assigned to julia to match the above info.
4) The slurm out files from the HPC will dump in the <runscripts/> folder so they don't clutter the entire directory. Delete these as needed to avoid pileup.