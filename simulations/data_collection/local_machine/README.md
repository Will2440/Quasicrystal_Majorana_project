# Local Machine Simulations

## Overview
This folder <simulations/data_collection/local_machine/> is designed for generating simulation data for the quasicrystal Kitaev chain model on a local machine. It has been written for Mac and not tested on other platforms. User-interaction can be limited to __init__.jl and main.jl to run pre-existing calculations. There is a UserOptions architecture within main.jl which is designed to allow the user to choose exactly what calculations are made, what data is saved and what existing solving function should be used for their task.

## Folder Contents
This folder contains the following scripts
- __init__.jl
  - This script should be run once per Julia session to load all of the modules to be used in the simulation. To make changes or    reset this, restart the Julia REPL.
- main.jl
  - This script is the primary point of user-interaction. It contains the code for parameter range choice, data saving location and options for tailoring the calculation to the users specific needs. Follow the usage instructions in the lead comments of this script for more details.
- compute.jl
  - This script contains the auxilliary function needed to read UserOptions defined in Sec 3 of main.jl and execute the corresponding functions in solvers.jl
- solvers.jl
  - This script contains all of the physics. Functions pertaining to the creation and diagonalisation of the Kitaev Hamiltonian, as well as the further calculations on the solved system are here.
  - Secondly, Sec 3 contains the different possible 'solving functions' which loop over parameters in different ways. These are optimised for different kinds of jobs. See the examples in the lead comments of this script for more details.

## Usage Notes

1) Run __init__.jl (once per Julia REPL session)
2) Copmplete all 3 sections of main.jl
3) Run main.jl

## Editing Notes

- To add or change the possible calculations made on the system edit functions in Sec 2 of solvers.jl
- To add or change the UserOptions architecture, you must check the following:
  - The struct defininton in solver.jl
  - The get_user_options() function in Sec 3 of main.jl
  - The run_selected_solver() function arguments and internal logic. Find the function in compute.jl and the function call in Sec 3 of main.jl
- Once functions within solver.jl are changed the Julia REPL must be restarted and __init__.jl run again to refresh the module
- To achieve smarter parameter space search for large computation mu vs. rho plots, create an auxilliary function to generate parameter ranges and parse those into main.jl, or rewrite a smarter solver function which automated this during computation


## Work in Progress

[25/09/2025]
- Only the hp and np_generic_solver functions are tested fully integrated with the UserOptions architecture. The mu_loop and N_loop functions need to be updated.
  - Want to change hp and np_generic_solver functions so that the logic for what is calculated and saved happens once outside the for loop, not each time within it.
- Need to add the Arpack routine versions of the solving function.
