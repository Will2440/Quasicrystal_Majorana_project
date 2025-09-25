# Data Processing

## Overview
This folder <simulations/data_processing> contains dedicated scripts for making the data generated in <data_collection/> readable (conversion from folder of .bson files to usable DataFrames) and doing the necessary manipulation to make plots of interest or further calculations.

## Folder Contents
This folder contains the following scripts
- bson_unpacker.jl
  - This script contains the main function -- unpack_bason_standard() -- and its auxilliaries needed to unpack the structure of .bson files produced by the <data_collection> scripts.
- standard_plotting.jl
  - This script contains generic plotting functions for 1D lines and 2D heatmaps where the variables can be chosen on the function call using Symbols == the DataFrame Key for that variable.
  - There are further functions which execute standard applications of these generic functions in one go for a quick view of the results.
- abundance.jl
  - This script contains the essential code to calculate and save the adundance of MBS phase.
- ipr_fitting.jl
  - This script contains the essential code to generate N vs. ipr plots with different fits (exponential, 1/N or none) to datasets corresponding to different values of mu. This is useful to demonstrate state localisation.


## Usage Notes

- Run bson_unpacker.jl by itself to display the DataFrame in the terminal for a quick check.
- Use the unpack_bason_standard() function from bson_unpacker.jl in any other plotting or processing function to obtain the standard form of the DataFrame for that dataset in the present working memory of the script.