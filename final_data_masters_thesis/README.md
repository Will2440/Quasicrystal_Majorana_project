This folder contains the final data used to produce the heatmap and IPR/MP line plots in the masters thesis report. This data was collected using BlueCrystal and by Miguel using BluePebble. 

Contents:
- Folders beginning with the crystal name (GQC, SQC, TMQC, PQC) contain data for ranges defined in the directory names. There are many .bson files in each fodler due to the nature of the data collect, these are collated by the code file described below.
- Ignore folder named 'N_loops'
- File names <bson_unpacker_sharable.jl> contains the code used to unpack the data and generate the relevant plots. See file header for more detail on usage
  - N.B. the updated version of this unpacking file, which is fully compatible with the current data_collection methods in <simulations/data_collection>, can be found in <simulations/data_processing>.

