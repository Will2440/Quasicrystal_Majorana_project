# Quasicrystal_Majorana_project


## Original Outlook

In this work we will study the existence of Majorana modes in quasicrystalline systems. To do so we will start by studying and analysing the paradigmatic example of the 1D periodic Kitaev chain. In doing so, we will familiarise ourselves with the BdG formalism, the characterisation of a topological phase via topological invariants and the LDOS. 

Once we have fully characterised the 1D periodic Kitaev chain we will study the effect of aperiodicity in the 1D chain (e.g., a Fibonacci distribution of hoppings along the 1D chain, a Fibonacci distribution on hoppings and pairings, and other aperiodic distributions) to see if the Majorana modes are present in the absence of translational symmetry and if the topological phase transitions happen at the same points in the phase diagram and if they are characterised by the same topological invariants.

Finally, if we have time, we will extend our study to 2D systems, such as the AB tiling, Penrose tiling or, more interesting, the Hat/Spectre.


## Working Summary of Completed Content

The bulk of the exploratory work carried out in this project is contained in Jupyter Notebooks within <simulations/working_notebooks> and <simulations/simple_conceptual_models>. Here the 1D Kitaev chain model is generated generically for any desired sequence of site hopping potential $t$. The physically equivalent $mu$ varying system is analysed separately.

These notebooks work through the following: sequence generation; BdG formalism for solving superconducting Hamiltonians; the Pfaffian topological invariant for periodic chains; plotting the energy structure from the eigenstates; plotting the local density of states (LDOS); computing and plotting the spectral function; inverse participation ratio (IPR) and the Majorana Polarisation (MP) as a topological indicator for the aperiodic chain. The abundance of topological phase in parameter space is explored. Finally the gap labelling theorem is used to demonstrate fractality in the energy structure and attempts at calculating the fractal dimension ar emade.

This further work on the fractal dimension focused on both the topological phase space and the energy structure of a single system. The latter was already established in the literature on Kitaev chains, however, the former pointed towards a fractal topological phase transition which is not analysed in the literature, and is thus a focus of the continued project.

Various working scripts were created to optimise solving our model both on local machines (<simulations/system_eigensolver>) and on HPC (Bristol's BlueCrystal was used) (<simulations/Majorana_solver_for_BlueCrystal>).

Results and figures from all of these scripts are collated in <simulations/hp_results> (standing for high precision -- this relates to the precision of the float variables used in the calculation needed to avoid erroneous reults at near-zero values), <simulations/images_and_graphs> with final versions collecetd in <simulations/report_plots>.


## Ongoing and Planned work

Edited 19/09/25
The current focus is on the calculation of the fractal dimanesions of the topological phase space and the energy structure. There is a working folder for optimised data cllection for the fristof these two tasks in <simulations/phase_transition_fractality>. This workflow is incomplete.


## Usage

- To run your own simulations go to <simulations/data_collection>. If you are running these on the 
- For access to existing clean data from the final masters thesis submission look in <final_data_masters_thesis>. N.B. The data structure of these data will not always match the currently active bson_unpacker.jl functions. This function expects a full compliment of values recorded, whereas the old thesis data will vary. To process this data the same functions can be used, but care must be taken adjusting the DataFrame creation in process_bson_files() so that it matches the contents of the .bsons exactly.

